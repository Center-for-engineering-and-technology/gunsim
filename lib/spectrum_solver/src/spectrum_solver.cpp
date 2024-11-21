#pragma once

#include <spectrum_solver.h>

namespace spectrum_solver {

using namespace spectrum_solver_structs;

    std::variant<SpectrumResult, std::string> SpectrumSolver::solve() {
        // To Do: проверка заполненности массива

        std::vector<double> totalSignal(static_cast<size_t>(specSolverOptions.sigParams.sampleNum));
        // минимальное время задержки сигнала до точки наблюдения, для первоначального сравнения равно максимальному значению
        double minTimeShift = specSolverOptions.sigParams.sampleMax;

        for (size_t i = 0; i < specSolverParams.gunArray.size(); i++) {
            auto& gun = specSolverParams.gunArray[i];
            // To Do: считывание из файла сигнала единичной пушки
            auto result = getSignalFromFile(gun);
            if (std::holds_alternative<std::string>(result)) {
                return std::get<std::string>(result);
            }
            auto& fileData = *std::get<std::unique_ptr<GunSignalData>>(result);
            auto& gunSignal = fileData.gunSignal;
            auto& gunSigParams = fileData.gunSigParams;
            if (specSolverOptions.sigParams.sampleNum != gunSigParams.sampleNum) {
                // To Do: интерполяция считанных значений на заданную в интерфейсе сетку
                return std::string("Error: Different size of vectors. Please, change parameters count in signal.");
            }
            // добавляем влияние расстояния до точки наблюдения на сигнал
            {
                std::string s = addObservationPointInfluence(gun, gunSignal);
                if (!s.empty()) { return s; }
            }
            // добавление отражения сигнала
            if (checkReflection()) {
                std::string s = addReflection(gun, gunSignal);
                if (!s.empty()) { return s; }
            }
            // добавление сигнала от текущей пушки к итоговому сигналу
            addNewGunSignal(minTimeShift, totalSignal, gun, gunSignal);
        }

        spectrum_solver_structs::SpectrumResult result;
        result.signal = totalSignal;
        result.sigParams = specSolverOptions.sigParams;

        // создаем Фурье солвер для расчетов спектра
        fourier_solver::ConcreteFourierSolver solver(specSolverOptions.sigParams);
        // расчет спектра
        solver.solve(totalSignal);
        result.specParams = solver.getSpecParams();

        // To Do: добавить std::variant для нескольких типов фильтра
        if (specSolverOptions.filter.bandpass.mode == gund_structs::BandpassFilter::OFF) {
            // получение результатов для нефильтрованного сигнала (амплитуда и фаза)
            result.ampSpec = solver.getAmp();
            // To Do: result.phaseSpec = solver.getPhase();
        }
        else if (specSolverOptions.filter.bandpass.mode == gund_structs::BandpassFilter::EXTERNAL) {
            addFilter(solver, result);
        }

        // To Do: вычисление параметров модели (peak-to-peak, rms и тд)
        computeToPeakValues(result);

        return result;
    }

    void SpectrumSolver::setDataPath() {
        try {
            dataPath = gund_utility::getExecutableDir() + "/data/";
        }
        catch(...) {}
    }

    std::variant<std::unique_ptr<GunSignalData>, std::string> SpectrumSolver::getSignalFromFile(const gund_structs::Gun& gun) {
        try {
            // вид строки: 1500C_6m_V100_P2000.sig
            // название-пушки_глубина-в-метрах_объем_давление.формат-для-сигнала
            std::string signalFilename = gund_json_parser::reconvertGunType(gun.type) + "_";
            signalFilename += (std::to_string((int)gun.z) + "m_");
            signalFilename += ("V" + std::to_string((int)gun.volume) + "_");
            signalFilename += ("P" + std::to_string((int)gun.pressure));
            signalFilename += ".sig";
            const auto file = dataPath + signalFilename;
            
            // чтение сигнала
            gund_format_parser::GundalfOutputParser parser(file);
            const std::vector<double>& signal = parser.getData();
            const gund_structs::SignatureParameters& sigParams = parser.getSigParams();
            return std::make_unique<GunSignalData>(signal, sigParams);
        }
        catch (const std::exception& e) {
            return std::string(e.what());
        }
    }

    bool SpectrumSolver::checkFilter() {
        if (specSolverOptions.filter.bandpass.mode == gund_structs::BandpassFilter::OFF) {
            return false;
        }
        else {
            return true;
        }
    }

    void SpectrumSolver::addFilter(fourier_solver::ConcreteFourierSolver& solver, spectrum_solver_structs::SpectrumResult& result) {
        const std::string fltFile = specSolverOptions.filter.bandpass.filename;
        // зачитывание характеристики фильтра
        gund_format_parser::GundalfOutputParser fltParser(fltFile);
        std::vector<double> filter = fltParser.getData();

        fourier_solver::FilterFourierSolver filteredSolver(solver, filter);
        filteredSolver.applyFilter();

        result.filteredSignal = filteredSolver.getSignal();
        result.ampSpec = filteredSolver.getAmp();
        // To Do: result.phaseSpec = filteredSolver.getPhase();
    }

    bool SpectrumSolver::checkObservationPoint() {
        if (specSolverOptions.obsPoint.observationType == gund_structs::ObservationPoint::ObservationPointType::INFINITE) {
            return true;
        }
        return false;
    }

    std::string SpectrumSolver::addObservationPointInfluence(const gund_structs::Gun& gun, std::vector<double>& gunSignal) {
        std::string err;
        if (checkObservationPoint()) {}
        // To Do: проверить, чем заполняются поля в json, есдли значение пусто
        else {
            double distanceToObs = computeDistanceFromGunToObsPoint(gun);
            if (distanceToObs == 0) {
                err = "Detector coordinates are equil with gun coordinates";
                return err;
            }
            double scale = 1. / distanceToObs;
            std::transform(gunSignal.begin(), gunSignal.end(), gunSignal.begin(),
                [&scale](double element) { return element *= scale; });
        }
        return err;
    }

    bool SpectrumSolver::checkReflection() {
        // To Do: проверка заполненности структуры (при наличии большего числа параметров для отражения)
        if (specSolverOptions.reflection.refCoef == 0) {
            return false;
        }
        return true;
    }

    std::string SpectrumSolver::addReflection(const gund_structs::Gun& gun, std::vector<double>& gunSignal) {     
        try {
            // смещение по времени
            double timeShift = computeRefTimeShift(gun);
            // соответствующее смещение по элементам вектора
            size_t valueShift = static_cast<size_t>(timeShift / specSolverOptions.sigParams.sampleInterval);
            std::vector<double> refSignal(gunSignal.size() - valueShift);
            std::copy(gunSignal.begin(), gunSignal.end() - valueShift, refSignal.begin());
            if (!checkObservationPoint()) {
                double distanceToObs = computeDistanceFromGunToObsPoint(gun);
                double ghostDistanceToObs = std::sqrt(std::pow(gun.x - specSolverOptions.obsPoint.x, 2)
                    + std::pow(gun.y - specSolverOptions.obsPoint.y, 2)
                    + std::pow(specSolverOptions.obsPoint.z + gun.z, 2));

                double scale = distanceToObs / ghostDistanceToObs;
                std::transform(refSignal.begin(), refSignal.end(), refSignal.begin(),
                    [&scale](double element) { return element *= scale; });
            }
            for (size_t i = valueShift; i < gunSignal.size(); i++) {
                gunSignal[i] += refSignal[i - valueShift] * specSolverOptions.reflection.refCoef;
            }
            return std::string();
        }
        catch (const std::exception& e) {
            return std::string(e.what());
        }
    }

    double SpectrumSolver::computeRefTimeShift(const gund_structs::Gun& gun) {
    
        if (checkObservationPoint()) {
            // т. наблюдения находится в 1 м от пушки
            // To Do: проверить т. наблюдения - бесконечность для нескольких пушек
            return 2 * gun.z / specSolverParams.physParams.soundVelocity;
        }
        else {
            /*
            // оценочное решение из статьи от создателей Gundalf - учитывает диаграмму направленности
            // из-за отсутствия учета диаграммы направленности выбрано направление в сторону детектора
            double y = computeDistanceFromGunToObsPoint(gun);
            double x = std::abs(gun.z - specSolverOptions.obsPoint.z);

            if (y == 0) {
                throw std::runtime_error("Detector coordinates are equil with gun coordinates");
            }

            return 2 * gun.z * (x / y) / specSolverParams.physParams.soundVelocity;*/

            // решение из прямого вычисления длины траектории сигнала
            double distanceToObs = computeDistanceFromGunToObsPoint(gun);
            double ghostDistanceToObs = std::sqrt(std::pow(gun.x - specSolverOptions.obsPoint.x, 2)
                + std::pow(gun.y - specSolverOptions.obsPoint.y, 2)
                + std::pow(specSolverOptions.obsPoint.z + gun.z, 2));
            return (ghostDistanceToObs - distanceToObs) / specSolverParams.physParams.soundVelocity;
        }
    }

    double SpectrumSolver::computeTimeShiftForSum(const gund_structs::Gun& gun) {
        // для точки наблюдения на бесконечности не важна зависимость от расстояния
        if (checkObservationPoint()) {
            return 0.;
        }
        else {
            // прямое вычисление длины траектории сигнала
            double distanceToObs = computeDistanceFromGunToObsPoint(gun);
            return distanceToObs / specSolverParams.physParams.soundVelocity;
        }
    }

    void SpectrumSolver::addNewGunSignal(double& minTimeShift, std::vector<double>& totalSignal, const gund_structs::Gun& gun, const std::vector<double>& gunSignal) {
        double gunTimeShift = computeTimeShiftForSum(gun);
        // To Do: проверка размера массива
        if (gunTimeShift == 0 || gunTimeShift == minTimeShift) {
            std::transform(totalSignal.begin(), totalSignal.end(), gunSignal.begin(), totalSignal.begin(), std::plus<double>());
        }
        else if (gunTimeShift > minTimeShift) {
            // соответствующее относительное смещение по элементам вектора
            size_t valueShift = static_cast<size_t>((gunTimeShift - minTimeShift) / specSolverOptions.sigParams.sampleInterval);
            for (size_t i = valueShift; i < totalSignal.size(); i++) {
                totalSignal[i] += gunSignal[i - valueShift];
            }
        }
        else {
            // соответствующее относительное смещение по элементам вектора
            size_t valueShift = static_cast<size_t>((minTimeShift - gunTimeShift) / specSolverOptions.sigParams.sampleInterval);
            minTimeShift = gunTimeShift;
            std::vector<double> tmpTotalSignal(totalSignal.size() - valueShift);
            std::copy(totalSignal.begin(), totalSignal.end() - valueShift, tmpTotalSignal.begin());
            totalSignal = gunSignal;
            for (size_t i = valueShift; i < totalSignal.size(); i++) {
                totalSignal[i] += tmpTotalSignal[i - valueShift];
            }
        }
    }

    double SpectrumSolver::computeDistanceFromGunToObsPoint(const gund_structs::Gun& gun) {
        double distance = std::sqrt(std::pow(gun.x - specSolverOptions.obsPoint.x, 2)
            + std::pow(gun.y - specSolverOptions.obsPoint.y, 2)
            + std::pow(gun.z - specSolverOptions.obsPoint.z, 2));
        return distance;
    }

    void SpectrumSolver::computeToPeakValues(spectrum_solver_structs::SpectrumResult& result) {
        std::vector<double> tmp_vec;
        if (checkFilter()) {
            if (result.filteredSignal.empty()) {
                throw std::runtime_error("Filtered signal data is empty");
            }
            tmp_vec = result.filteredSignal;
        }
        else {
            if (result.signal.empty()) {
                throw std::runtime_error("Signal data is empty");
            }
            tmp_vec = result.signal;
        }
        auto [min, max] = std::minmax_element(tmp_vec.begin(), tmp_vec.end());
        // значение между максимальным и минимальным значением экстремумов (в Бар-м)
        result.rParams.peakToPeakBar = *max - *min;
        // для нахождения в dB = 20 * log10(val) + 220
        result.rParams.peakToPeakdB = 20 * std::log10(result.rParams.peakToPeakBar) + 220;
        // значение между максимальным экстремумом и нулевым значением
        result.rParams.zeroToPeakBar = *max;
        result.rParams.zeroToPeakdB = 20 * std::log10(result.rParams.zeroToPeakBar) + 220;
        // To Do: погрешности (зависят от конкретной пушки и можно хранить в БД) - не зависят от наличия фильтра
    }


} // namespace spectrum_solver
