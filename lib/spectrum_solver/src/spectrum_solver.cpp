#include <algorithm>
#define _USE_MATH_DEFINES

#include <signature_structs.h>

#include <filter.h>
#include <spectrum_solver.h>

#include <random>

namespace spectrum_solver
{

using namespace spectrum_solver_structs;

std::variant<SpectrumResult, std::string> SpectrumSolver::solve ()
{
    // To Do: проверка заполненности массива
    // To Do: добавить и сюда, если у фильтра проблемы с шагом - считаем без фильтра
    if (!checkObservationCoordinate())
    {
        // координаты точки наблюдения совпадают в координатами одной из пушек
        errors.messages.push_back("Coordinates of observation point coincide with coordinates of one of the guns");
        errors.messages.push_back("Continue simulation with observation point at infinity");
        specSolverOptions.obsPoint.observationType = gund_structs::ObservationPoint::ObservationPointType::INFINITE;
    }
    if (specSolverOptions.reflection.firstCableDepth != 0 && specSolverOptions.reflection.firstCableDepth == specSolverOptions.reflection.secondCableDepth)
    {
        // в Gundalf при совпадении глубины кабелей, второй игнорируется
        errors.messages.push_back("The cable depths are equal.");
        errors.messages.push_back("Continue simulation ignoring cable 2");
        specSolverOptions.reflection.secondCableDepth = 0.;
    }
    // заполняем соседей всех пушек (для drop-out и первоначальной генерации сигнала)
    fillNeighborGuns();
    // заполняем максимальную глубину среди всех пушек (для т. наблюдения на бесконечности)
    fillMaxDepth();

    if (specSolverOptions.diffModel)
    {
        gun_model::GunArraySolver<> solver(gunMap, specSolverParams, specSolverOptions, gun_model::GunArraySolver<>::InteractionModel::DEFAULT_INTERACTION);
        auto res = solver.solve();
        if (res != "Solving done well")
        {
            gun_model::GunArraySolver<> solver_close(gunMap, specSolverParams, specSolverOptions, gun_model::GunArraySolver<>::InteractionModel::CLOSE_INTERACTION);
            auto res_close = solver_close.solve();
            if (res_close != "Solving done well")
                return "Error from gun_model: " + res_close;
            else
                gun_model_result = std::make_unique<gun_model::GunModelResult>(solver_close.getResult());
        }
        else
        {
            gun_model_result = std::make_unique<gun_model::GunModelResult>(solver.getResult());
        }
    }

    std::vector<double> totalSignal(static_cast<size_t>(specSolverOptions.sigParams.sampleNum));
    // минимальное время задержки сигнала до точки наблюдения, для первоначального сравнения равно максимальному значению
    specSolverTmp.minTimeShift = specSolverOptions.sigParams.sampleMax;

    for (size_t i = 0; i < specSolverParams.gunArray.size(); i++)
    {
        auto& gun = specSolverParams.gunArray[i];
        auto itDublicate(std::find(specSolverTmp.dublicateGuns.begin(), specSolverTmp.dublicateGuns.end(), i));
        if (itDublicate != specSolverTmp.dublicateGuns.end())
        {
            // текущая пушка находится в той же точке, что и уже учтенная в расчете
            errors.messages.push_back("In position (" + std::to_string(gun.x) + ", " + std::to_string(gun.y) + ", " + std::to_string(gun.z) + ") gun already exists. Gun " + gun.name + " ignored");
            // сохранение нулевого вектора сигнала, для того, чтобы не сбивать нумерацию
            std::vector<double> nullGunSignal = {};
            fillIndividualSignals(0, nullGunSignal);
            continue;
        }
        // моделируем сигнал одиночного пневмоисточника с отражением и влиянием т. наблюдения
        std::variant<std::unique_ptr<std::vector<double>>, std::string> modeledSignal;
        if (specSolverOptions.diffModel)
            modeledSignal = modelIndividualGunSignal_diffModel(gun, i, *gun_model_result);
        else
            modeledSignal = modelIndividualGunSignal(gun);

        if (std::holds_alternative<std::string>(modeledSignal))
        {
            return std::get<std::string>(modeledSignal);
        }
        auto& gunSignal = *std::get<std::unique_ptr<std::vector<double>>>(modeledSignal);

        // обновление суммарных параметров при прохождении по циклу
        updateSumValues(gun, gunSignal);
        // добавление сигнала от текущей пушки к итоговому сигналу
        addNewGunSignal(specSolverTmp.minTimeShift, totalSignal, gun, gunSignal);
    }

    if (!std::isnan(specSolverOptions.sampleTimeShift))
    {
        // TODO: Очень костыльный сдвиг сигнатуры по методике от МАГЭ
        // TODO: Сдвигать через набег фазы + сохранять размер
        auto find_peak_idx = [&] () -> long long
        {
            if (totalSignal.empty())
            {
                return 0;
            }

            const auto max_elem = std::max_element(
                std::cbegin(totalSignal),
                std::cend(totalSignal)
            );

            if (std::cend(totalSignal) == max_elem)
            {
                return 0;
            }

            const auto peak_threshold = 0.8 * (*max_elem);
            const auto idx_elem = std::find_if(
                std::cbegin(totalSignal),
                std::cend(totalSignal),
                [peak_threshold] (double value) -> bool
                {
                    return value >= peak_threshold;
                }
            );

            if (std::cend(totalSignal) == idx_elem)
            {
                return 0;
            }

            return static_cast<long long>(idx_elem - totalSignal.begin());
        };

        auto time_shift_size = static_cast<long long>(specSolverOptions.sampleTimeShift / specSolverOptions.output_sigParams.sampleInterval) - find_peak_idx();
        if (time_shift_size > 0)
        {
            auto orig_size = totalSignal.size();
            std::vector<double> time_shift_zeros_padding(time_shift_size);
            totalSignal.insert(totalSignal.begin(), time_shift_zeros_padding.begin(), time_shift_zeros_padding.end());
            totalSignal.resize(orig_size);
        }
    }

    SpectrumResult result;
    result.sigParams = specSolverOptions.output_sigParams;
    result.signal = totalSignal;
    // создаем Фурье солвер для расчетов спектра
    fourier_solver::ConcreteFourierSolver solver(result.sigParams);
    solver.solve(result.signal);
    result.specParams = solver.getSpecParams();

    addFilter(solver, result);
    // To Do: сохранить сдвиг и сигналы в tmp

    // вычисление вспомогательных параметров модели (peak-to-peak, rms и тд)
    computeToPeakValues(result);
    computeContribAndErrorValues(result);
    computeCenterValues(result);
    computeEnergyParams(result);

    return result;
}

std::variant<std::unique_ptr<std::vector<double>>, std::string> SpectrumSolver::modelIndividualGunSignal (const gund_structs::Gun& gun)
{
    // считывание из файла сигнала единичной пушки
    auto result = getSignalFromFile(gun);
    if (std::holds_alternative<std::string>(result))
    {
        return std::get<std::string>(result);
    }
    auto& fileData = *std::get<std::unique_ptr<GunSignalData>>(result);
    auto gunSignal = std::make_unique<std::vector<double>>(fileData.gunSignal);
    auto& gunSigParams = fileData.gunSigParams;
    if (specSolverOptions.sigParams.sampleNum != gunSigParams.sampleNum)
    {
        // To Do: интерполяция считанных значений на заданную в интерфейсе сетку
        errors.criticalErrors = ("Error: Different size of vectors.Please, change parameters count in signal\n");
        return errors.criticalErrors;
    }
    // добавляем влияние расстояния до точки наблюдения на сигнал
    {
        std::string s = addObservationPointInfluence(gun, *gunSignal);
        if (!s.empty())
        {
            return s;
        }
    }
    // добавление отражения сигнала
    if (checkReflection())
    {
        std::string s = addReflection(gun, *gunSignal);
        if (!s.empty())
        {
            return s;
        }
    }
    return gunSignal;
}

std::variant<std::unique_ptr<std::vector<double>>, std::string> SpectrumSolver::modelIndividualGunSignal_diffModel (const gund_structs::Gun& gun, size_t i, const gun_model::GunModelResult& gmr)
{

    if (i >= gmr.signatures.size())
    {
        errors.criticalErrors = ("Error: Diff model didn't work correctly\n");
        return errors.criticalErrors;
    }
    auto gunSignal = std::make_unique<std::vector<double>>(gmr.signatures[i]);
    if (specSolverOptions.sigParams.sampleNum != gunSignal->size())
    {
        // To Do: интерполяция считанных значений на заданную в интерфейсе сетку
        errors.criticalErrors = ("Error: Different size of vectors.Please, change parameters count in signal\n");
        return errors.criticalErrors;
    }
    // добавляем влияние расстояния до точки наблюдения на сигнал
    {
        std::string s = addObservationPointInfluence(gun, *gunSignal);
        if (!s.empty())
        {
            return s;
        }
    }
    // добавление отражения сигнала
    if (checkReflection())
    {
        std::string s = addReflection(gun, *gunSignal);
        if (!s.empty())
        {
            return s;
        }
    }
    return gunSignal;
}

void SpectrumSolver::setDataPath ()
{
    try
    {
        dataPath = gund_utility::getExecutableDir() + "/data/";
    }
    catch (...)
    {
    }
}

std::variant<std::unique_ptr<GunSignalData>, std::string> SpectrumSolver::getSignalFromFile (const gund_structs::Gun& gun)
{
    try
    {
        // вид строки: 1500C_6m_V100_P2000.sig
        // название-пушки_глубина-в-метрах_объем_давление.формат-для-сигнала
        std::string signalFilename = gund_json_parser::reconvertGunType(gun.type) + "_";
        signalFilename += (std::to_string(( int )gun.z) + "m_");
        signalFilename += ("V" + std::to_string(( int )gun.volume) + "_");
        signalFilename += ("P" + std::to_string(( int )gun.pressure));
        signalFilename += ".sig";
        const auto file = dataPath + signalFilename;

        // чтение сигнала
        gund_format_parser::GundalfOutputParser parser(file);
        const std::vector<double>& signal = parser.getData();
        const gund_structs::SignatureParameters& sigParams = parser.getSigParams();
        return std::make_unique<GunSignalData>(signal, sigParams);
    }
    catch (const std::exception& e)
    {
        return std::string(e.what());
    }
}

void SpectrumSolver::fillRelErrors ()
{
    try
    {
        std::string relErrorsFile = "relative_errors.json";
        const auto file = dataPath + relErrorsFile;
        hashRelErrors = gund_json_parser::parseErrorsJson(file);
    }
    catch (...)
    {
    }
}

bool SpectrumSolver::checkFilter ()
{
    if (specSolverOptions.filter.bandpass.mode == gund_structs::BandpassFilter::OFF)
    {
        return false;
    }
    else
    {
        return true;
    }
}

std::string SpectrumSolver::addFilter (fourier_solver::ConcreteFourierSolver& solver, spectrum_solver_structs::SpectrumResult& result)
{
    try
    {
        std::vector<double> filtered_signal;
        const auto& sig_params = specSolverOptions.sigParams;
        const auto& out_sig_params = specSolverOptions.output_sigParams;

        if (specSolverOptions.pre_filter.bandpass.mode == gund_structs::BandpassFilter::EXTERNAL)
        {
            const std::string& fltFile = specSolverOptions.pre_filter.bandpass.filename;
            if (!fltFile.empty())
            {
                // зачитывание характеристики фильтра
                gund_format_parser::GundalfOutputParser fltParser(dataPath + fltFile);
                specSolverOptions.pre_filter.bandpass.filterData = fltParser.getData();
                specSolverOptions.pre_filter.bandpass.filterSigParams = fltParser.getSigParams();
            }
            if (specSolverOptions.pre_filter.bandpass.filterSigParams.sampleInterval
                != sig_params.sampleInterval)
            {
                // предварительно, пока мы не делаем интерполяцию значений (Gundalf выдает ту же ошибку)
                throw std::runtime_error(
                    "Pre Filter has diferent sample interval than the modelling sample interval."
                    "Please select another filter."
                );
            }

            fourier_solver::FilterFourierSolver filteredSolver(
                solver,
                specSolverOptions.pre_filter.bandpass.filterSigParams,
                specSolverOptions.pre_filter.bandpass.filterData
            );
            filteredSolver.applyFilter();
            filtered_signal = std::move(filteredSolver.getSignal());
        }
        else
        {
            filtered_signal = solver.getSignal();
        }

        if (specSolverOptions.antialias_filter.has_value()
            || (sig_params.sampleInterval != out_sig_params.sampleInterval))
        {
            size_t modelling_to_out_scale = static_cast<size_t>(
                std::ceil(out_sig_params.sampleInterval / sig_params.sampleInterval)
            );
            if (specSolverOptions.antialias_filter.has_value()
                && specSolverOptions.antialias_filter->getHighNyqPrcntFreqCut()
                < 99)
            {
                const auto& config = specSolverOptions.antialias_filter.value();

                gund_structs::NewSignatureParameters new_sig_params = gund_structs::ConvertOldToNewSignatureParameters(out_sig_params);

                gund_structs::RealSignature signal(std::move(new_sig_params), 1);

                signal.data = std::move(filtered_signal);

                signal = filters_common::ApplyFilter(
                    signal,
                    filters_common::constructFilter(config, result.sigParams.sampleNum, new_sig_params)
                );

                filtered_signal.resize(result.sigParams.sampleNum);
                for (size_t i = 0; i < result.sigParams.sampleNum; ++i)
                {
                    filtered_signal[i] = signal.data[modelling_to_out_scale * i];
                }
            }
            else
            {
                std::vector<double> out_signal(result.sigParams.sampleNum);
                for (size_t i = 0; i < result.sigParams.sampleNum; ++i)
                {
                    out_signal[i] = filtered_signal[modelling_to_out_scale * i];
                }
                filtered_signal = std::move(out_signal);
            }

            solver.solve(filtered_signal);
        }

        if (specSolverOptions.filter.bandpass.mode == gund_structs::BandpassFilter::EXTERNAL)
        {
            const std::string& fltFile = specSolverOptions.filter.bandpass.filename;
            if (!fltFile.empty())
            {
                // зачитывание характеристики фильтра
                gund_format_parser::GundalfOutputParser fltParser(dataPath + fltFile);
                specSolverOptions.filter.bandpass.filterData = fltParser.getData();
                specSolverOptions.filter.bandpass.filterSigParams = fltParser.getSigParams();
            }
            if (specSolverOptions.filter.bandpass.filterSigParams.sampleInterval != out_sig_params.sampleInterval)
            {
                // предварительно, пока мы не делаем интерполяцию значений (Gundalf выдает ту же ошибку)
                throw std::runtime_error("Filter has diferent sample interval than the signature. Please select another filter.");
            }

            fourier_solver::FilterFourierSolver filteredSolver(
                solver,
                specSolverOptions.filter.bandpass.filterSigParams,
                specSolverOptions.filter.bandpass.filterData
            );
            filteredSolver.applyFilter();
            filtered_signal = std::move(filteredSolver.getSignal());

            auto solver_cp = solver;
            solver_cp.solve(specSolverOptions.filter.bandpass.filterData);
            auto filter_spectrum = solver_cp.getSpectrum();

            result.bandpassfilterPhaseSpec = fourier_solver::computePhase(filter_spectrum);

            auto max_spectrum_val = std::abs(
                *std::max_element(filter_spectrum.begin(), filter_spectrum.end(), [] (auto a, auto b)
                                  {
                                      return std::abs(a) < std::abs(b);
                                  })
            );
            std::transform(filter_spectrum.begin(), filter_spectrum.end(), filter_spectrum.begin(), [&] (auto spec_val)
                           {
                               return std::abs(spec_val) / max_spectrum_val;
                           });
            result.bandpassfilterAmpSpec = fourier_solver::computeAmp(filter_spectrum);
        }

        if (specSolverOptions.q_filter.has_value())
        {
            const auto& config = specSolverOptions.q_filter.value();

            gund_structs::NewSignatureParameters new_sig_params = gund_structs::ConvertOldToNewSignatureParameters(out_sig_params);

            gund_structs::RealSignature signal(std::move(new_sig_params), 1);
            signal.data = std::move(filtered_signal);

            signal = filters_common::ApplyFilter(
                signal,
                filters_common::constructFilter(config, new_sig_params)
            );
            filtered_signal = std::move(signal.data);
        }

        if (specSolverOptions.wiener_filter.has_value())
        {
            const auto& config = specSolverOptions.wiener_filter.value();

            gund_structs::NewSignatureParameters new_sig_params = gund_structs::ConvertOldToNewSignatureParameters(out_sig_params);

            gund_structs::RealSignature signal(new_sig_params, 1);
            signal.data = std::move(filtered_signal);

            auto wiener_filter = filters_common::constructFilter(config, signal);

            auto wiener_fitlered = filters_common::ApplyFilter(signal, wiener_filter);
            filtered_signal = std::move(wiener_fitlered.data);

            if (config.scale_mode)
            {
                const auto& in_raw_data = signal.data;
                auto& out_raw_data = filtered_signal;
                auto in_norm = std::accumulate(in_raw_data.begin(), in_raw_data.end(), 0., [] (auto a, auto b)
                                               {
                                                   return a + b * b;
                                               });
                auto out_norm = std::accumulate(out_raw_data.begin(), out_raw_data.end(), 0., [] (auto a, auto b)
                                                {
                                                    return a + b * b;
                                                });
                if (out_norm > 0)
                {
                    auto scale = std::sqrt(in_norm / out_norm);

                    for (auto& out_v: out_raw_data)
                        out_v *= scale;
                }
            }
        }

        solver.solve(filtered_signal);
        result.filteredSignal = std::move(filtered_signal);
        result.ampSpec = std::move(solver.getAmp());
        result.phaseSpec = std::move(solver.getPhase());

        return std::string();
    }
    catch (const std::exception& e)
    {
        return std::string(e.what());
    }
}

bool SpectrumSolver::checkObservationPoint ()
{
    if (specSolverOptions.obsPoint.observationType == gund_structs::ObservationPoint::ObservationPointType::INFINITE)
    {
        return true;
    }
    return false;
}

std::string SpectrumSolver::addObservationPointInfluence (const gund_structs::Gun& gun, std::vector<double>& gunSignal)
{
    std::string err;
    if (checkObservationPoint())
    {
    }
    else
    {
        double distanceToObs = computeDistanceFromGunToObsPoint(gun);
        if (distanceToObs == 0)
        {
            err = "Detector coordinates are equil with gun coordinates";
            return err;
        }
        double scale = 1. / distanceToObs;
        std::transform(gunSignal.begin(), gunSignal.end(), gunSignal.begin(), [&scale] (double element)
                       {
                           return element *= scale;
                       });
    }
    return err;
}

bool SpectrumSolver::checkReflection ()
{
    // To Do: проверка заполненности структуры (при наличии большего числа параметров для отражения)
    if (specSolverOptions.reflection.refCoef == 0)
    {
        return false;
    }
    return true;
}

std::string SpectrumSolver::addReflection (const gund_structs::Gun& gun, std::vector<double>& gunSignal)
{
    try
    {
        // смещение по времени
        double timeShift = computeRefTimeShift(gun);
        // соответствующее смещение по элементам вектора
        size_t valueShift = static_cast<size_t>(timeShift / specSolverOptions.sigParams.sampleInterval);
        std::vector<double> refSignal(gunSignal.size() - valueShift);
        std::copy(gunSignal.begin(), gunSignal.end() - valueShift, refSignal.begin());
        if (!checkObservationPoint())
        {
            double distanceToObs = computeDistanceFromGunToObsPoint(gun);
            double ghostDistanceToObs = std::sqrt(
                std::pow(gun.x - specSolverOptions.obsPoint.x, 2)
                + std::pow(gun.y - specSolverOptions.obsPoint.y, 2)
                + std::pow(specSolverOptions.obsPoint.z + gun.z, 2)
            );

            double scale = distanceToObs / ghostDistanceToObs;
            std::transform(refSignal.begin(), refSignal.end(), refSignal.begin(), [&scale] (double element)
                           {
                               return element *= scale;
                           });
        }
        else
        {
            // отражение от кабеля только в случае точки наблюдения на бесконечности
            std::vector<double> cableRefSignal(gunSignal.size(), 0);
            if (specSolverOptions.reflection.firstCableDepth != 0)
            {
                // сценарий source->cable1->surface->reciever (ref)
                // в случае расположения кабеля над ПИ Gundalf все равно считает пройденное расстояние как 2 глубины кабеля. Поступаем аналогично
                addReflectionLineFromCable(cableRefSignal, gunSignal, specSolverOptions.reflection.firstCableDepth, 1);
                // сценарий source->surface->cable1->surface->reciever (ref^2)
                addReflectionLineFromCable(cableRefSignal, gunSignal, specSolverOptions.reflection.firstCableDepth + gun.z, 2);
                // остальными сценариями распространения сигнала Gundalf пренебрегает
            }
            if (specSolverOptions.reflection.secondCableDepth != 0)
            {
                // сценарий source->cable2->surface->reciever (ref)
                // в случае расположения кабеля над ПИ Gundalf все равно считает пройденное расстояние как 2 глубины кабеля. Поступаем аналогично
                addReflectionLineFromCable(cableRefSignal, gunSignal, specSolverOptions.reflection.secondCableDepth, 1);
                // сценарий source->surface->cable2->surface->reciever (ref^2)
                addReflectionLineFromCable(cableRefSignal, gunSignal, specSolverOptions.reflection.secondCableDepth + gun.z, 2);
            }
            if (specSolverOptions.reflection.firstCableDepth != 0 && specSolverOptions.reflection.secondCableDepth != 0)
            {
                // сценарий source->cable_min->surface->cable_max->surface->reciever (ref^2)
                addReflectionLineFromCable(cableRefSignal, gunSignal, specSolverOptions.reflection.firstCableDepth + specSolverOptions.reflection.secondCableDepth, 2);
                // сценарий source->surface->cable_min->surface->cable_max->surface->reciever (ref^3)
                addReflectionLineFromCable(cableRefSignal, gunSignal, specSolverOptions.reflection.firstCableDepth + specSolverOptions.reflection.secondCableDepth + gun.z, 3);
            }
            std::transform(gunSignal.begin(), gunSignal.end(), cableRefSignal.begin(), gunSignal.begin(), std::plus<double>());
        }
        std::transform(gunSignal.begin() + valueShift, gunSignal.end(), refSignal.begin(), gunSignal.begin() + valueShift, [&] (double a, double b)
                       {
                           return a + b * specSolverOptions.reflection.refCoef;
                       });
        return std::string();
    }
    catch (const std::exception& e)
    {
        return std::string(e.what());
    }
}

void SpectrumSolver::addReflectionLineFromCable (std::vector<double>& gunCableSignal, const std::vector<double>& gunSignal, const double signalPath, const int power)
{
    double cableTimeShift = 2 * signalPath / specSolverParams.physParams.soundVelocity;
    // соответствующее смещение по элементам вектора
    size_t cableValueShift = static_cast<size_t>(cableTimeShift / specSolverOptions.sigParams.sampleInterval);
    std::vector<double> tmpCableRefSignal(gunSignal.size() - cableValueShift);
    std::copy(gunSignal.begin(), gunSignal.end() - cableValueShift, tmpCableRefSignal.begin());
    std::transform(gunCableSignal.begin() + cableValueShift, gunCableSignal.end(), tmpCableRefSignal.begin(), gunCableSignal.begin() + cableValueShift, [&] (double a, double b)
                   {
                       return a + b * std::pow(specSolverOptions.reflection.refCoef, power);
                   });
}

double SpectrumSolver::computeRefTimeShift (const gund_structs::Gun& gun)
{

    if (checkObservationPoint())
    {
        // т. наблюдения находится в 1 м от пушки
        return 2 * gun.z / specSolverParams.physParams.soundVelocity;
    }
    else
    {
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
        double ghostDistanceToObs = std::sqrt(
            std::pow(gun.x - specSolverOptions.obsPoint.x, 2)
            + std::pow(gun.y - specSolverOptions.obsPoint.y, 2)
            + std::pow(specSolverOptions.obsPoint.z + gun.z, 2)
        );
        return (ghostDistanceToObs - distanceToObs) / specSolverParams.physParams.soundVelocity;
    }
}

double SpectrumSolver::computeTimeShiftForSum (const gund_structs::Gun& gun)
{
    // для точки наблюдения на бесконечности не важна зависимость от расстояния,
    // расчет производится относительно самой глубокой пушки
    if (checkObservationPoint())
    {
        if (gun.z == specSolverTmp.maxDepth)
        {
            return 0.;
        }
        else
        {
            return (specSolverTmp.maxDepth - gun.z) / specSolverParams.physParams.soundVelocity;
        }
    }
    else
    {
        // прямое вычисление длины траектории сигнала
        double distanceToObs = computeDistanceFromGunToObsPoint(gun);
        return distanceToObs / specSolverParams.physParams.soundVelocity;
    }
}

void SpectrumSolver::addNewGunSignal (
    double& minTimeShift,
    std::vector<double>& totalSignal,
    const gund_structs::Gun& gun,
    const std::vector<double>& gunSignal
)
{
    // смещение за счет задержки сигнала до гидрофона
    double gunTimeShift = computeTimeShiftForSum(gun);
    // смещение за счет разного времени старта испускания сигнала
    gunTimeShift += gun.delay;

    if (specSolverOptions.sourceControllerVariation > 0.)
    {
        std::mt19937 engine(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        std::normal_distribution<double> distribution(0., specSolverOptions.sourceControllerVariation);
        // смещение за счет вариации времени пуска от контроллера
        // может быть только положительной - т.к. сигнал на пушку не может прийти раньше чем на контроллер
        gunTimeShift += std::abs(distribution(engine));
    }

    if (specSolverTmp.individualSignals.size() < specSolverParams.gunArray.size())
    {
        // сохранение индивидуальных сигналов (уже после интерполяции, влияния отражения и т. наблюдения)
        // смещение сохраняется отдельно от сигнала
        fillIndividualSignals(gunTimeShift, gunSignal);
    }
    // To Do: проверка размера массива
    if (gunTimeShift == minTimeShift)
    {
        std::transform(totalSignal.begin(), totalSignal.end(), gunSignal.begin(), totalSignal.begin(), std::plus<double>());
    }
    else if (gunTimeShift > minTimeShift)
    {
        // соответствующее относительное смещение по элементам вектора
        size_t valueShift = static_cast<size_t>((gunTimeShift - minTimeShift) / specSolverOptions.sigParams.sampleInterval);
        for (size_t i = valueShift; i < totalSignal.size(); i++)
        {
            totalSignal[i] += gunSignal[i - valueShift];
        }
    }
    else
    {
        // соответствующее относительное смещение по элементам вектора
        size_t valueShift = static_cast<size_t>((minTimeShift - gunTimeShift) / specSolverOptions.sigParams.sampleInterval);
        minTimeShift = gunTimeShift;
        std::vector<double> tmpTotalSignal(totalSignal.size() - valueShift);
        std::copy(totalSignal.begin(), totalSignal.end() - valueShift, tmpTotalSignal.begin());
        totalSignal = gunSignal;
        for (size_t i = valueShift; i < totalSignal.size(); i++)
        {
            totalSignal[i] += tmpTotalSignal[i - valueShift];
        }
    }
}

double SpectrumSolver::computeDistanceBetweenGuns (const gund_structs::Gun& gunA, const gund_structs::Gun& gunB)
{
    double distance = std::sqrt(std::pow(gunA.x - gunB.x, 2) + std::pow(gunA.y - gunB.y, 2) + std::pow(gunA.z - gunB.z, 2));
    return distance;
}

double SpectrumSolver::computeDistanceFromGunToObsPoint (const gund_structs::Gun& gun)
{
    double distance = std::sqrt(std::pow(gun.x - specSolverOptions.obsPoint.x, 2) + std::pow(gun.y - specSolverOptions.obsPoint.y, 2) + std::pow(gun.z - specSolverOptions.obsPoint.z, 2));
    return distance;
}

void SpectrumSolver::fillIndividualSignals (const double gunTimeShift, const std::vector<double>& gunSignal)
{
    specSolverTmp.individualSignals.push_back(IndividualGunSignals(gunTimeShift, std::make_unique<std::vector<double>>(gunSignal)));
}

void SpectrumSolver::updateSumValues (const gund_structs::Gun& gun, std::vector<double>& gunSignal)
{
    // суммарный объем системы
    specSolverTmp.sumValues.sumVolume += gun.volume;

    std::vector<double> tmp_signal;
    if (!checkFilter())
    {
        tmp_signal = gunSignal;
    }
    else
    {
        const std::string& fltFile = specSolverOptions.filter.bandpass.filename;
        if (!fltFile.empty() && specSolverOptions.filter.bandpass.filterData.empty())
        {
            // зачитывание характеристики фильтра
            gund_format_parser::GundalfOutputParser fltParser(dataPath + fltFile);
            specSolverOptions.filter.bandpass.filterData = fltParser.getData();
            specSolverOptions.filter.bandpass.filterSigParams = fltParser.getSigParams();
        }
        // если есть фильтр, то требуется отфильтровывать сигнал всех единичных источников
        fourier_solver::ConcreteFourierSolver solver(specSolverOptions.sigParams);
        solver.solve(gunSignal);

        fourier_solver::FilterFourierSolver filteredSolver(solver, specSolverOptions.filter.bandpass.filterSigParams, specSolverOptions.filter.bandpass.filterData);
        filteredSolver.applyFilter();

        tmp_signal = filteredSolver.getSignal();
    }
    // центр давления основан на значениях peak-to-peak единичных источников
    auto [min, max] = std::minmax_element(tmp_signal.begin(), tmp_signal.end());
    double peakToPeak = *max - *min;
    specSolverTmp.sumValues.individualPeakToPeak.push_back(peakToPeak);

    specSolverTmp.sumValues.pressureCoordSum.x += peakToPeak * gun.x;
    specSolverTmp.sumValues.pressureCoordSum.y += peakToPeak * gun.y;
    specSolverTmp.sumValues.pressureCoordSum.z += peakToPeak * gun.z;
    specSolverTmp.sumValues.peakToPeakSum += peakToPeak;

    if (hashRelErrors.empty())
    {
        // заполнение хэш-таблицы с относительными ошибками
        fillRelErrors();
    }
    gund_structs::ErrorData relErrorData;

    auto it = hashRelErrors.find(gun.type);
    if (it != hashRelErrors.end())
    {
        relErrorData = gund_structs::ErrorData{it->second.peakToPeakErr, it->second.bubblePeakErr, it->second.bubblePeriodErr};
    }
    else
    {
        // значение по умолчанию (относительные ошибки еще не добавлены)
        relErrorData = hashRelErrors.at(gund_structs::GunType::C1500);
    }
    // рассматриваем абсолютные и относительные значения ошибок
    double ppAbsError = peakToPeak * relErrorData.peakToPeakErr;
    double zpAbsError = ppAbsError * 0.5;

    // TODO: убрать думблирование
    int bubble_search_time_shift = static_cast<int>(
        specSolverOptions.bubble_search_start_time / specSolverOptions.sigParams.sampleInterval
    );
    if (static_cast<int>(std::distance(tmp_signal.begin(), max)) > bubble_search_time_shift)
    {
        bubble_search_time_shift = static_cast<int>(
            std::distance(tmp_signal.begin(), max) + 0.04 / specSolverOptions.sigParams.sampleInterval
        );
    }

    auto max_bubble = std::max_element(tmp_signal.begin() + bubble_search_time_shift, tmp_signal.end());
    auto min_bubble = std::min_element(max_bubble, tmp_signal.end());
    double bubPeriod = std::distance(tmp_signal.begin(), max_bubble) * specSolverOptions.sigParams.sampleInterval;
    double bpAbsError = bubPeriod * relErrorData.bubblePeriodErr;

    double prbAbsError;
    double primaryToBubble = (*max - *min) / (*max_bubble - *min_bubble);
    prbAbsError = std::sqrt(
                      std::pow(relErrorData.peakToPeakErr, 2) + std::pow(relErrorData.bubblePeakErr * 0.5, 2)
                  )
                * primaryToBubble;

    specSolverTmp.sumValues.sumSqPeakToPeakError += ppAbsError * ppAbsError;
    specSolverTmp.sumValues.sumSqZeroToPeakError += zpAbsError * zpAbsError;
    specSolverTmp.sumValues.sumSqPrToBubError += prbAbsError * prbAbsError;
    specSolverTmp.sumValues.sumSqBubPeriodError += bpAbsError * bpAbsError;
}

void SpectrumSolver::computeToPeakValues (spectrum_solver_structs::SpectrumResult& result)
{
    std::vector<double> tmp_vec;
    if (checkFilter())
    {
        if (result.filteredSignal.empty())
        {
            throw std::runtime_error("Filtered signal data is empty");
        }
        tmp_vec = result.filteredSignal;
    }
    else
    {
        if (result.signal.empty())
        {
            throw std::runtime_error("Signal data is empty");
        }
        tmp_vec = result.signal;
    }
    auto [min, max] = std::minmax_element(tmp_vec.begin(), tmp_vec.end());
    int bubble_search_time_shift = static_cast<int>(specSolverOptions.bubble_search_start_time / specSolverOptions.sigParams.sampleInterval);
    if (static_cast<int>(std::distance(tmp_vec.begin(), max)) > bubble_search_time_shift)
    {
        bubble_search_time_shift = static_cast<int>(std::distance(tmp_vec.begin(), max) + 0.04 / specSolverOptions.sigParams.sampleInterval);
    }
    auto max_bubble = std::max_element(tmp_vec.begin() + bubble_search_time_shift, tmp_vec.end());
    auto min_bubble = std::min_element(max_bubble, tmp_vec.end());

    double rmsSq = std::accumulate(tmp_vec.begin(), tmp_vec.end(), 0., [] (double sum, double i)
                                   {
                                       return sum + i * i;
                                   })
                 / result.sigParams.sampleNum;

    // проверяем корректность окна пропускания
    if (specSolverOptions.passband.highFrequency >= result.specParams.sampleMax || specSolverOptions.passband.highFrequency <= specSolverOptions.passband.lowFrequency || specSolverOptions.passband.lowFrequency < 0.)
    {
        errors.messages.push_back("Invalid passband limits. Continue simulation with default spectral passband [10 Hz, 50 Hz].");
        specSolverOptions.passband.lowFrequency = 10.;
        specSolverOptions.passband.highFrequency = 50.;
    }
    size_t start_box_index = static_cast<size_t>(specSolverOptions.passband.lowFrequency / (result.specParams.sampleInterval * result.sigParams.sampleNum));
    size_t end_box_index = static_cast<size_t>(specSolverOptions.passband.highFrequency / (result.specParams.sampleInterval * result.sigParams.sampleNum));
    auto [min_box_spec, max_box_spec]
        = std::minmax_element(result.ampSpec.begin() + start_box_index, result.ampSpec.begin() + end_box_index + 1);
    double averBoxSpec = std::accumulate(result.ampSpec.begin() + start_box_index, result.ampSpec.begin() + end_box_index + 1, 0.)
                       / (end_box_index + 1 - start_box_index);

    // значение между максимальным и минимальным значением экстремумов (в Бар-м)
    result.rParams.peakToPeakBar = *max - *min;
    // для нахождения в dB = 20 * log10(val) + 220
    result.rParams.peakToPeakdB = 20 * std::log10(result.rParams.peakToPeakBar) + 220;
    // значение между максимальным экстремумом и нулевым значением
    result.rParams.zeroToPeakBar = *max;
    result.rParams.zeroToPeakdB = 20 * std::log10(result.rParams.zeroToPeakBar) + 220;
    // значение для PB зависит от наличия отражения и фильтрации
    // первый минимум суммарного сигнала - максимум отраженной волны с противоположной фазой, а не начало влияния пузыря
    result.rParams.primaryToBubble = (*max - *min) / (*max_bubble - *min_bubble);
    // время от начала сигнала до второго пика (из-за пузыря)
    result.rParams.bubblePeriodToFirstPeak = std::distance(tmp_vec.begin(), max_bubble) * result.sigParams.sampleInterval;
    // значение rms: rms^2 = 1/T \int^T_0 p(t)^2 dt
    result.rParams.rmsPressureBar = std::sqrt(rmsSq);
    result.rParams.rmsPressuredB = 20 * std::log10(result.rParams.rmsPressureBar) + 220;
    // значения для спектра в дБ
    result.rParams.maximumSpectralValue = *max_box_spec;
    result.rParams.averageSpectralValue = averBoxSpec;
    result.rParams.maximumSpectralRipple = *max_box_spec - *min_box_spec;
}

void SpectrumSolver::computeContribAndErrorValues (spectrum_solver_structs::SpectrumResult& result)
{
    // из-за наличия delay у сигнала, а также из-за погрешности алгоритмических вычислений, суммарный peak-to-peak не равен сумме индивидуальных
    // требуется найти пропорциональный индивидуальному вклад ПИ в получившийся суммарный peak-to-peak
    double realSumPeakToPeak = std::accumulate(specSolverTmp.sumValues.individualPeakToPeak.begin(), specSolverTmp.sumValues.individualPeakToPeak.end(), 0., [] (double sum, double i)
                                               {
                                                   return sum + i;
                                               });
    if (realSumPeakToPeak == 0. || specSolverParams.gunArray.size() != specSolverTmp.sumValues.individualPeakToPeak.size())
    {
        errors.messages.push_back("Individual peak-to-peak calculation is failed.");
    }
    double scale = 1. / realSumPeakToPeak;
    result.singleContrib.resize(specSolverParams.gunArray.size());
    for (size_t i = 0; i < specSolverParams.gunArray.size(); i++)
    {
        result.singleContrib[i].gunName = specSolverParams.gunArray[i].name;
        // при поиске вклада в peak-to-peak итогового сигнала, result.rParams.peakToPeakBar сокращается
        result.singleContrib[i].singlePeakToPeakPct = specSolverTmp.sumValues.individualPeakToPeak[i] * scale * 100.;
    }

    // погрешности (не зависят от наличия фильтра)
    result.rParams.peakToPeakError = std::sqrt(specSolverTmp.sumValues.sumSqPeakToPeakError);
    result.rParams.zeroToPeakError = std::sqrt(specSolverTmp.sumValues.sumSqZeroToPeakError);
    result.rParams.primaryToBubbleError = std::sqrt(specSolverTmp.sumValues.sumSqPrToBubError);
    result.rParams.bubblePeriodToFirstPeakError = std::sqrt(specSolverTmp.sumValues.sumSqBubPeriodError);
}

void SpectrumSolver::computeCenterValues (spectrum_solver_structs::SpectrumResult& result)
{
    result.rParams.sumVolume = specSolverTmp.sumValues.sumVolume;
    // ищем геометрический центр как центр прямоугольника из максимальных и минимальных значений координат активных пушек
    auto& tmpGunArray = specSolverParams.gunArray;
    auto [min_x, max_x] = std::minmax_element(tmpGunArray.begin(), tmpGunArray.end(), [] (const gund_structs::Gun& gunA, const gund_structs::Gun& gunB)
                                              {
                                                  return gunA.x < gunB.x;
                                              });
    auto [min_y, max_y] = std::minmax_element(tmpGunArray.begin(), tmpGunArray.end(), [] (const gund_structs::Gun& gunA, const gund_structs::Gun& gunB)
                                              {
                                                  return gunA.y < gunB.y;
                                              });
    auto [min_z, max_z] = std::minmax_element(tmpGunArray.begin(), tmpGunArray.end(), [] (const gund_structs::Gun& gunA, const gund_structs::Gun& gunB)
                                              {
                                                  return gunA.z < gunB.z;
                                              });
    result.rParams.geometricCenter = gund_structs::Coordinate{
        (max_x->x - min_x->x) * 0.5 + min_x->x,
        (max_y->y - min_y->y) * 0.5 + min_y->y,
        (max_z->z - min_z->z) * 0.5 + min_z->z
    };
    // ищем центр давления как значение накопленной суммы peak-to-peak давления к значению peak-to-peak суммарного сигнала
    // всегда рассматриваем суммарный нефильтрованный сигнал
    std::vector<double> tmp_vec;
    if (checkFilter())
    {
        tmp_vec = result.filteredSignal;
    }
    else
    {
        tmp_vec = result.signal;
    }
    auto [min, max] = std::minmax_element(tmp_vec.begin(), tmp_vec.end());
    // double scale = 1. / ((*max - *min) * specSolverParams.gunArray.size());
    double scale = 1. / specSolverTmp.sumValues.peakToPeakSum;
    // предварительно: пока модель основана на суммировании единичных независимых источников возникает различие для одинаковых значений (сумма единичных p-p_i больше p-p реального)
    double press_x = (min_x->x == max_x->x) ? min_x->x : specSolverTmp.sumValues.pressureCoordSum.x * scale;
    double press_y = (min_y->y == max_y->y) ? min_y->y : specSolverTmp.sumValues.pressureCoordSum.y * scale;
    double press_z = (min_z->z == max_z->z) ? min_z->z : specSolverTmp.sumValues.pressureCoordSum.z * scale;
    result.rParams.pressureCenter = gund_structs::Coordinate{press_x, press_y, press_z};
}

void SpectrumSolver::computeEnergyParams (spectrum_solver_structs::SpectrumResult& result)
{
    std::vector<energy_comp::DelayedGun> dGun;
    std::transform(specSolverParams.gunArray.begin(), specSolverParams.gunArray.end(), std::back_inserter(dGun), [] (const gund_structs::Gun& gun)
                   {
                       energy_comp::DelayedGun temp{gun, gun.delay};
                       return temp;
                   });

    energy_comp::EnergyParams energyParams;
    energyParams.gunArray = dGun;
    energyParams.physParams = {specSolverParams.physParams};

    energy_comp::EnergyOptions energyOptions;
    energyOptions.reflection = specSolverOptions.reflection;
    energyOptions.filter = specSolverOptions.filter;
    energyOptions.sigParams = specSolverOptions.sigParams;

    energy_comp::Energy energySolver(energyParams, energyOptions);
    energySolver.setDataPath(dataPath);
    if (specSolverOptions.diffModel)
        energySolver.prepareSignalsAndLengths(*gun_model_result);
    else
        energySolver.prepareSignalsAndLengths();
    energySolver.solve();

    energy_comp::EnergyResult energyResult = energySolver.getResult();

    result.rParams.totalAcousticEnergy = energyResult.totalAcousticEnergy;
    result.rParams.totalAcousticEfficiency = energyResult.acousticEnergyEffectiveness;
    result.rParams.totalPotentialEnergy = energyResult.totalPotentialEnergy;
    result.rParams.energyCenter = gund_structs::Coordinate{energyResult.energyCenter.x, energyResult.energyCenter.y, energyResult.energyCenter.z};
    // заполнение индивидуальной энергии ПИ
    result.singleContrib.resize(specSolverParams.gunArray.size());
    if (result.singleContrib.size() != energyResult.energy.size())
    {
        throw std::runtime_error("Calculation of individual energy is failed.");
    }
    for (size_t i = 0; i < result.singleContrib.size(); i++)
    {
        if (result.singleContrib[i].gunName.empty())
        {
            result.singleContrib[i].gunName = specSolverParams.gunArray[i].name;
        }
        result.singleContrib[i].singleEnergy = energyResult.energy[i];
    }
    return;
}

std::string SpectrumSolver::computeDiagrams (const std::string& outDir, const std::unique_ptr<gun_model::GunModelResult>& gmr)
{
    try
    {
        std::vector<direct_diag::DelayedGun> dGun;
        std::transform(specSolverParams.gunArray.begin(), specSolverParams.gunArray.end(), std::back_inserter(dGun), [] (const gund_structs::Gun& gun)
                       {
                           direct_diag::DelayedGun temp{gun, gun.delay};
                           return temp;
                       });

        direct_diag::DiagramParams diagParams;
        diagParams.gunArray = dGun;
        diagParams.physParams = specSolverParams.physParams;

        direct_diag::DiagramOptions diagOptions;
        diagOptions.reflection = specSolverOptions.reflection;
        diagOptions.filter = specSolverOptions.filter;
        diagOptions.sigParams = specSolverOptions.sigParams;

        direct_diag::Diagram diagSolver(diagParams, diagOptions);
        diagSolver.setDataPath(dataPath);
        if (specSolverOptions.diffModel)
            diagSolver.prepareSpecs(*gmr);
        else
            diagSolver.prepareSpecs();
        // меняем директорию для сохранения данных
        diagSolver.setDataPath(outDir);
        // вывод двумерных диаграмм направленности
        direct_diag::OutputDiagramOptions ODO;
        ODO.dipIncr = specSolverOptions.dirParams.dipIncr;
        ODO.lowerDB = specSolverOptions.dirParams.lowerDB;
        ODO.higherDB = specSolverOptions.dirParams.higherDB;
        ODO.maxDipAngle = specSolverOptions.dirParams.maxDipAngle;
        ODO.azimut = specSolverOptions.dirParams.inlineAzimutAngle;
        ODO.outputFileName = "inline_directivity.json";
        diagSolver.outputDiagram(ODO);
        ODO.azimut = specSolverOptions.dirParams.crosslineAzimutAngle;
        ODO.outputFileName = "crossline_directivity.json";
        diagSolver.outputDiagram(ODO);
        // вывод одномерных диаграмм направленности
        direct_diag::OutputSignatureDiagramOptions OSDO;
        OSDO.dipIncrForSignRepr = specSolverOptions.dirParams.dipIncrForSignRepr;
        OSDO.azimut = specSolverOptions.dirParams.inlineAzimutAngle;
        OSDO.outputFileName = "inline_directivity.json";
        diagSolver.outputSignatureDiagram(OSDO);
        OSDO.azimut = specSolverOptions.dirParams.crosslineAzimutAngle;
        OSDO.outputFileName = "crossline_directivity.json";
        diagSolver.outputSignatureDiagram(OSDO);
        return std::string();
    }
    catch (const std::exception& e)
    {
        return std::string(e.what());
    }
}

std::string SpectrumSolver::computeAzimutalDiagrams (const std::string& outDir, const std::unique_ptr<gun_model::GunModelResult>& gmr)
{
    try
    {
        if (!specSolverOptions.dirParams.azimutal)
        {
            return std::string();
        }
        std::vector<direct_diag::DelayedGun> dGun;
        std::transform(specSolverParams.gunArray.begin(), specSolverParams.gunArray.end(), std::back_inserter(dGun), [] (const gund_structs::Gun& gun)
                       {
                           direct_diag::DelayedGun temp{gun, gun.delay};
                           return temp;
                       });

        direct_diag::DiagramParams diagParams;
        diagParams.gunArray = dGun;
        diagParams.physParams = specSolverParams.physParams;

        direct_diag::DiagramOptions diagOptions;
        diagOptions.reflection = specSolverOptions.reflection;
        diagOptions.filter = specSolverOptions.filter;
        diagOptions.sigParams = specSolverOptions.sigParams;

        direct_diag::Diagram diagSolver(diagParams, diagOptions);
        diagSolver.setDataPath(dataPath);
        if (specSolverOptions.diffModel)
            diagSolver.prepareSpecs(*gmr);
        else
            diagSolver.prepareSpecs();
        // меняем директорию для сохранения данных
        diagSolver.setDataPath(outDir);
        // вывод азимутальных диаграмм направленности
        direct_diag::OutputAzimutalDiagramOptions ODO;
        ODO.lowerDB = specSolverOptions.dirParams.azimutalLowerDB;
        ODO.higherDB = specSolverOptions.dirParams.azimutalHigherDB;
        ODO.frequency = specSolverOptions.dirParams.azimutalFreq1;
        ODO.outputFileName = "azimutal_directivity_1.json";
        diagSolver.outputAzimutalDiagram(ODO);
        ODO.frequency = specSolverOptions.dirParams.azimutalFreq2;
        ODO.outputFileName = "azimutal_directivity_2.json";
        diagSolver.outputAzimutalDiagram(ODO);
        ODO.frequency = specSolverOptions.dirParams.azimutalFreq3;
        ODO.outputFileName = "azimutal_directivity_3.json";
        diagSolver.outputAzimutalDiagram(ODO);
        ODO.frequency = specSolverOptions.dirParams.azimutalFreq4;
        ODO.outputFileName = "azimutal_directivity_4.json";
        diagSolver.outputAzimutalDiagram(ODO);

        return std::string();
    }
    catch (const std::exception& e)
    {
        return std::string(e.what());
    }
}

bool SpectrumSolver::checkGunCoordinate (size_t vecPosition)
{
    for (size_t i = 0; i < vecPosition; i++)
    {
        if (std::abs(specSolverParams.gunArray[i].x - specSolverParams.gunArray[vecPosition].x) < 0.1 && std::abs(specSolverParams.gunArray[i].y - specSolverParams.gunArray[vecPosition].y) < 0.1 && std::abs(specSolverParams.gunArray[i].z - specSolverParams.gunArray[vecPosition].z) < 0.1)
        {
            return false;
        }
    }
    return true;
}

bool SpectrumSolver::checkObservationCoordinate ()
{
    if (!checkObservationPoint())
    {
        for (size_t i = 0; i < specSolverParams.gunArray.size(); i++)
        {
            if (std::abs(specSolverParams.gunArray[i].x - specSolverOptions.obsPoint.x) < 0.1 && std::abs(specSolverParams.gunArray[i].y - specSolverOptions.obsPoint.y) < 0.1 && std::abs(specSolverParams.gunArray[i].z - specSolverOptions.obsPoint.z) < 0.1)
            {
                return false;
            }
        }
    }
    return true;
}

void SpectrumSolver::fillMaxDepth ()
{
    if (checkObservationPoint())
    {
        auto max = std::max_element(specSolverParams.gunArray.begin(), specSolverParams.gunArray.end(), [] (const gund_structs::Gun& gunA, const gund_structs::Gun& gunB)
                                    {
                                        return gunA.z < gunB.z;
                                    });
        specSolverTmp.maxDepth = max->z;
    }
}

void SpectrumSolver::fillNeighborGuns ()
{
    size_t gunNumber = specSolverParams.gunArray.size();
    if (specSolverTmp.gunsNeighbors.empty())
    {
        specSolverTmp.gunsNeighbors.resize(gunNumber);
    }
    // минимальное расстояние, на котором считаем взаимодействие между источниками нулевым
    double minInfluenceDist = 10.; // magic_const
    // минимальное расстояние, на котором могут находиться источники
    double minDistance = std::sqrt(3 * 0.1 * 0.1); // magic_const
    for (size_t i = 0; i < gunNumber - 1; i++)
    {
        for (size_t j = i + 1; j < gunNumber; j++)
        {
            double dist = computeDistanceBetweenGuns(specSolverParams.gunArray[i], specSolverParams.gunArray[j]);
            if (dist < minDistance)
            {
                specSolverTmp.dublicateGuns.push_back(j);
            }
            // если пушки находятся слишком близко, то они не учитываются (дополнительная валидация)
            else if (dist < minInfluenceDist)
            {
                specSolverTmp.gunsNeighbors[i].push_back(j);
                specSolverTmp.gunsNeighbors[j].push_back(i);
            }
        }
    }
}

std::string SpectrumSolver::removeGunSignal (std::vector<double>& totalSignal, const size_t gunIndex)
{
    // рассматриваем сохраненные ранее значения
    const std::vector<double>& gunSignal = *specSolverTmp.individualSignals[gunIndex].gunSignal;
    double& gunTimeShift = specSolverTmp.individualSignals[gunIndex].timeShift;
    // рассматриваем сдвиг по времени для всего сигнала (в заданной т. наблюдения начало сигнала сдвигается относительно самого "быстрого" сигнала)
    double& resultTimeShift = specSolverTmp.minTimeShift;
    // ПИ игнорируется
    if (gunSignal.empty())
    {
        return std::string();
    }
    if (gunSignal.size() != totalSignal.size())
    {
        return std::string("Incorrectly saved individual signal.");
    }
    // сдвигаем сохраненный сигнал, заполняя остальное нулями
    if (gunTimeShift == resultTimeShift)
    {
        std::transform(totalSignal.begin(), totalSignal.end(), gunSignal.begin(), totalSignal.begin(), std::minus<double>());
    }
    else if (gunTimeShift > resultTimeShift)
    {
        // соответствующее относительное смещение по элементам вектора
        size_t valueShift = static_cast<size_t>((gunTimeShift - resultTimeShift) / specSolverOptions.sigParams.sampleInterval);
        for (size_t i = valueShift; i < totalSignal.size(); i++)
        {
            totalSignal[i] -= gunSignal[i - valueShift];
        }
    }
    else
    {
        // такого быть не может, так как resultTimeShift - минимальный сдвиг
        return std::string("Incorrectly shifted individual signal.");
    }
    // конечный сигнал НЕ смещается, даже если удален наиболее важный для первого пика сигнал
    return std::string();
}

std::variant<SpectrumResult, std::string> SpectrumSolver::computeDropSignal (const std::vector<size_t>& dropGunIndexes, const std::vector<double>& fullSignal)
{
    // сохраняем копию сигнала для последующего изменения
    std::vector<double> signalAfterDrop = fullSignal;
    // создаем объединение соседей
    std::vector<size_t> allNeighbors;
    // убираем сигнал ПИ, участвующей в drop-out
    for (auto gunIndex: dropGunIndexes)
    {
        std::string s = removeGunSignal(signalAfterDrop, gunIndex);
        if (!s.empty())
        {
            return s;
        }

        // добавляем вектор соседей рассматриваемой пушки в объединение всех соседей данного drop-out
        auto gunNeighborsVector = specSolverTmp.gunsNeighbors[gunIndex];
        allNeighbors.insert(allNeighbors.end(), gunNeighborsVector.begin(), gunNeighborsVector.end());
    }
    // чтобы избежать повторного вычитания сигнала соседей, сортируем их и удаляем общих
    std::sort(allNeighbors.begin(), allNeighbors.end());
    auto repeated = std::unique(allNeighbors.begin(), allNeighbors.end());
    allNeighbors.erase(repeated, allNeighbors.end());
    // удаляем источники, участвующие в отсеве из соседей
    auto dropped = std::remove_if(allNeighbors.begin(), allNeighbors.end(), [&] (size_t elem)
                                  {
                                      auto itGunDrop(std::find(dropGunIndexes.begin(), dropGunIndexes.end(), elem));
                                      return itGunDrop != dropGunIndexes.end();
                                  });
    allNeighbors.erase(dropped, allNeighbors.end());

    std::unique_ptr<gun_model::GunModelResult> neighborGunModel;
    if (specSolverOptions.diffModel)
    {
        // создаем вектор только соседних источников
        SpectrumSolverParams tmpParams = specSolverParams;
        tmpParams.gunArray = std::vector<gund_structs::Gun>(allNeighbors.size());
        for (size_t i = 0; i < allNeighbors.size(); i++)
        {
            tmpParams.gunArray[i] = specSolverParams.gunArray[allNeighbors[i]];
        }
        gun_model::GunArraySolver<> solver(gunMap, tmpParams, specSolverOptions, gun_model::GunArraySolver<>::InteractionModel::DEFAULT_INTERACTION);
        auto res = solver.solve();
        if (res != "Solving done well")
        {
            gun_model::GunArraySolver<> solver_close(gunMap, tmpParams, specSolverOptions, gun_model::GunArraySolver<>::InteractionModel::CLOSE_INTERACTION);
            auto res_close = solver_close.solve();
            if (res_close != "Solving done well")
                return "Error from gun_model: " + res_close;
            else
                neighborGunModel = std::make_unique<gun_model::GunModelResult>(solver_close.getResult());
        }
        else
        {
            neighborGunModel = std::make_unique<gun_model::GunModelResult>(solver.getResult());
        }
    }

    // убираем старые сигналы (со старым влиянием между ПИ) и добавляем новые (без влияния выключенного источника)
    // проходим по всем соседям drop-out
    for (size_t i = 0; i < allNeighbors.size(); i++)
    {
        auto neighborIndex = allNeighbors[i];
        // удаляем старый сигнал соседа
        std::string s = removeGunSignal(signalAfterDrop, neighborIndex);
        if (!s.empty())
        {
            return s;
        }

        const gund_structs::Gun& neighborGun = specSolverParams.gunArray[neighborIndex];
        // добавляем обновленный сигнал соседа к результирующему сигналу
        auto modeledSignal = modelIndividualGunSignal(neighborGun);
        if (specSolverOptions.diffModel)
            // порядок источников в allNeighbors совпадает с порядком в neighborGunModel
            modeledSignal = modelIndividualGunSignal_diffModel(neighborGun, i, *neighborGunModel);
        else
            modeledSignal = modelIndividualGunSignal(neighborGun);

        if (std::holds_alternative<std::string>(modeledSignal))
        {
            return std::get<std::string>(modeledSignal);
        }
        auto& neighborSignal = *std::get<std::unique_ptr<std::vector<double>>>(modeledSignal);
        addNewGunSignal(specSolverTmp.minTimeShift, signalAfterDrop, neighborGun, neighborSignal); // To Do: check save signature
    }
    // To Do: ниже копия того, что делается в методе solve
    // cоздать отдельную функцию под это
    SpectrumResult dropResult;
    dropResult.signal = signalAfterDrop;
    dropResult.sigParams = specSolverOptions.output_sigParams;

    // создаем Фурье солвер для расчетов спектра
    fourier_solver::ConcreteFourierSolver solver(specSolverOptions.sigParams);
    // расчет спектра
    solver.solve(signalAfterDrop);
    dropResult.specParams = solver.getSpecParams();

    addFilter(solver, dropResult);
    // вычисление вспомогательных параметров модели (peak-to-peak, rms и тд)
    computeToPeakValues(dropResult);

    return dropResult;
}

DropOut SpectrumSolver::compareDropAndOrigin (const std::vector<size_t>& dropGunIndexes, const SpectrumResult& dropResult, const SpectrumResult& allGunsResult)
{
    DropOut dropOut;

    for (auto gunIndex: dropGunIndexes)
    {
        if (!dropOut.firstGunDrop)
        {
            dropOut.firstGunDrop = &specSolverParams.gunArray[gunIndex];
        }
        else if (!dropOut.secondGunDrop)
        {
            dropOut.secondGunDrop = &specSolverParams.gunArray[gunIndex];
        }
        else if (!dropOut.thirdGunDrop)
        {
            dropOut.thirdGunDrop = &specSolverParams.gunArray[gunIndex];
        }
        else
        {
            throw std::runtime_error("Incorrect data entry.");
        }
    }
    dropOut.peakToPeak = dropResult.rParams.peakToPeakBar;
    dropOut.peakToPeakPercentDrop = 100 * (1. - dropResult.rParams.peakToPeakBar / allGunsResult.rParams.peakToPeakBar);
    dropOut.normXCorrDrop = gund_utility::calculateNormalXCorrelation(dropResult.signal, allGunsResult.signal);

    size_t startBoxIndex = static_cast<size_t>(specSolverOptions.passband.lowFrequency / dropResult.specParams.sampleInterval);
    size_t endBoxIndex = static_cast<size_t>(specSolverOptions.passband.highFrequency / dropResult.specParams.sampleInterval);
    auto maxDropIter = std::max_element(dropResult.ampSpec.begin() + startBoxIndex, dropResult.ampSpec.begin() + endBoxIndex + 1, [&] (const double& a, const double& b)
                                        {
                                            return std::abs(a - allGunsResult.ampSpec[&a - dropResult.ampSpec.data()]) < std::abs(b - allGunsResult.ampSpec[&b - dropResult.ampSpec.data()]);
                                        });

    dropOut.averageAmpDrop = allGunsResult.rParams.averageSpectralValue - dropResult.rParams.averageSpectralValue;
    auto maxDropIndex = std::distance(dropResult.ampSpec.begin(), maxDropIter);
    dropOut.maxAmpDrop = std::abs(*maxDropIter - allGunsResult.ampSpec[maxDropIndex]);
    dropOut.freqMaxAmpDrop = maxDropIndex * dropResult.specParams.sampleInterval;

    if (specSolverOptions.dropParams.minPrimaryToBubble != 0)
    {
        dropOut.primaryToBubble = dropResult.rParams.primaryToBubble;
        dropOut.primaryToBubbleSucceed = dropOut.primaryToBubble > specSolverOptions.dropParams.minPrimaryToBubble;
    }
    else
    {
        dropOut.primaryToBubbleDrop = 100 * (1. - dropResult.rParams.primaryToBubble / allGunsResult.rParams.primaryToBubble);
        dropOut.primaryToBubbleSucceed = dropOut.primaryToBubbleDrop < specSolverOptions.dropParams.maxPrimaryToBubbleDrop;
    }

    dropOut.peakToPeakSucceed = dropOut.peakToPeakPercentDrop < specSolverOptions.dropParams.maxPeakToPeakDrop;
    dropOut.normXCorrSucceed = dropOut.normXCorrDrop > specSolverOptions.dropParams.minNormXCorr;
    dropOut.averageAmpSucceed = dropOut.averageAmpDrop < specSolverOptions.dropParams.maxAverageAmpDrop;
    dropOut.maxAmpSucceed = dropOut.maxAmpDrop < specSolverOptions.dropParams.maxAmpDifference;

    dropOut.succeedDropOut = (dropOut.peakToPeakSucceed && dropOut.primaryToBubbleSucceed && dropOut.normXCorrSucceed && dropOut.averageAmpSucceed && dropOut.maxAmpSucceed);

    return dropOut;
}

std::variant<DropOutAndSpectrumResults, std::string> SpectrumSolver::dropOut ()
{
    auto res = solve();
    spectrum_solver_structs::SpectrumResult allGunsResult;
    // обработка результата
    if (std::holds_alternative<spectrum_solver_structs::SpectrumResult>(res))
    {
        allGunsResult = std::get<spectrum_solver_structs::SpectrumResult>(res);
        if (allGunsResult.rParams.peakToPeakBar == 0 || allGunsResult.rParams.primaryToBubble == 0)
        {
            return std::string("Incorrect calculation result of engeneering parameters. Stop simulation to prevent division by zero.");
        }
    }
    else
    {
        return std::get<std::string>(res);
    }

    DropOutResult result;
    // для последующего параллельного вычисления отказываемся от рекурсии (или нет ??)
    // цикл для одного отключенного ПИ
    if (specSolverOptions.dropParams.mode >= 1)
    {
        for (size_t i = 0; i < specSolverParams.gunArray.size(); i++)
        {
            std::vector<size_t> dropGunsIndexes{i};
            SpectrumResult dropGunResult;
            auto dropRes = computeDropSignal(dropGunsIndexes, allGunsResult.signal);
            if (std::holds_alternative<spectrum_solver_structs::SpectrumResult>(dropRes))
            {
                dropGunResult = std::get<spectrum_solver_structs::SpectrumResult>(dropRes);
            }
            else
            {
                return std::get<std::string>(res);
            }
            auto dropParams = compareDropAndOrigin(dropGunsIndexes, dropGunResult, allGunsResult);
            result.singleDropOut.push_back(dropParams);
        }
        if (result.singleDropOut.empty())
        {
            return std::string("Incorrect single drop-out.");
        }
        result.singleSuccessfulPercent = std::accumulate(result.singleDropOut.begin(), result.singleDropOut.end(), 0., [] (double sum, const DropOut& drop)
                                                         {
                                                             return sum + drop.succeedDropOut;
                                                         })
                                       / result.singleDropOut.size()
                                       * 100;
    }
    // цикл для двух отключенных ПИ
    if (specSolverOptions.dropParams.mode >= 2)
    {
        for (size_t i = 0; i < specSolverParams.gunArray.size() - 1; i++)
        {
            for (size_t j = i + 1; j < specSolverParams.gunArray.size(); j++)
            {
                std::vector<size_t> dropGunsIndexes{i, j};
                SpectrumResult dropGunResult;
                auto dropRes = computeDropSignal(dropGunsIndexes, allGunsResult.signal);
                if (std::holds_alternative<spectrum_solver_structs::SpectrumResult>(dropRes))
                {
                    dropGunResult = std::get<spectrum_solver_structs::SpectrumResult>(dropRes);
                }
                else
                {
                    return std::get<std::string>(res);
                }
                auto dropParams = compareDropAndOrigin(dropGunsIndexes, dropGunResult, allGunsResult);
                result.doubleDropOut.push_back(dropParams);
            }
        }
        if (result.doubleDropOut.empty())
        {
            return std::string("Incorrect double drop-out.");
        }
        result.doubleSuccessfulPercent = std::accumulate(result.doubleDropOut.begin(), result.doubleDropOut.end(), 0., [] (double sum, const DropOut& drop)
                                                         {
                                                             return sum + drop.succeedDropOut;
                                                         })
                                       / result.doubleDropOut.size()
                                       * 100;
    }
    // цикл для трех отключенных ПИ
    if (specSolverOptions.dropParams.mode >= 3)
    {
        for (size_t i = 0; i < specSolverParams.gunArray.size() - 2; i++)
        {
            for (size_t j = i + 1; j < specSolverParams.gunArray.size() - 1; j++)
            {
                for (size_t k = j + 1; k < specSolverParams.gunArray.size(); k++)
                {
                    std::vector<size_t> dropGunsIndexes{i, j, k};
                    SpectrumResult dropGunResult;
                    auto dropRes = computeDropSignal(dropGunsIndexes, allGunsResult.signal);
                    if (std::holds_alternative<spectrum_solver_structs::SpectrumResult>(dropRes))
                    {
                        dropGunResult = std::get<spectrum_solver_structs::SpectrumResult>(dropRes);
                    }
                    else
                    {
                        return std::get<std::string>(res);
                    }
                    auto dropParams = compareDropAndOrigin(dropGunsIndexes, dropGunResult, allGunsResult);
                    result.tripleDropOut.push_back(dropParams);
                }
            }
        }
        if (result.tripleDropOut.empty())
        {
            return std::string("Incorrect triple drop-out.");
        }
        result.tripleSuccessfulPercent = std::accumulate(result.tripleDropOut.begin(), result.tripleDropOut.end(), 0., [] (double sum, const DropOut& drop)
                                                         {
                                                             return sum + drop.succeedDropOut;
                                                         })
                                       / result.tripleDropOut.size()
                                       * 100;
    }
    return std::make_pair(result, allGunsResult);
}

} // namespace spectrum_solver
