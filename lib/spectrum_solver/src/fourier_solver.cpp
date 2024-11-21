#include <fourier_solver.h>

#include <simple_fft/fft.h>

#include <algorithm>

namespace fourier_solver {

    std::vector<double> computeAmp(std::vector<complex_type> spectrum_) {
        if (spectrum_.empty()) {
            throw std::runtime_error("Error: Spectrum data is empty");
        }

        std::vector<double> ampSpec(spectrum_.size());
        for (size_t i = 0; i < spectrum_.size(); i++) {
            ampSpec[i] = 20 * std::log10(std::abs(spectrum_[i]));
        }
        return ampSpec;
    }

    const std::vector<double> ConcreteFourierSolver::getAmp() {
        return computeAmp(spectrum);
    }

    void ConcreteFourierSolver::solve(std::vector<double>& tot_sig) {
        signal.resize(tot_sig.size());
        signal = tot_sig;

        std::vector<double> updSignal(tot_sig.size());
        std::copy(tot_sig.begin(), tot_sig.end(), updSignal.begin());

        // To Do: количество нулей, как и получившийся шаг по времени (не по частоте!) для спектра, дополняют число до степени двойки
        // в случае нулей - количество элементов - степень двойки,
        // случай временного шага - шаг * 10е6 (10е7) кратен степени двойки, либо является таковым - 0.000512 или 0.0007168 (делится на 1024)
        // поэтому, предварительно
        const double updSampleInterval = 0.000512; // magic const

        std::vector<double> zeros(24, 0); // magic const
        // для эффективного использования БПФ требуется количество элементов, кратное степени двойки
        // добавляем нулевых 24 элемента в конец массива
        // UPD: добавить не в конец, а в начало
        updSignal.insert(updSignal.end(), zeros.begin(), zeros.end());

        // расчет спектра для сигнала с добавочными значениями
        spectrum = fastFourierTransform(updSignal);

        specParams.sampleNum = spectrum.size();
        specParams.sampleInterval = 1 / (updSampleInterval * 1000);
        specParams.sampleMax = specParams.sampleInterval * specParams.sampleNum;
        // расчет основных параметров модели
        //computeAmp();
        //computePhase(); // To Do
    }

    std::vector<complex_type> ConcreteFourierSolver::fastFourierTransform(std::vector<double>& tot_sig, double convCoef) {
        // переводим в мкПа
        double scale = sigParams.sampleInterval * convCoef;
        std::vector<double> transformSignal(tot_sig.size());
        std::transform(tot_sig.begin(), tot_sig.end(), transformSignal.begin(),
            [&scale](double element) { return element *= scale; });

        std::vector<complex_type> signalSpectrum(tot_sig.size());
        auto b = simple_fft::FFT(transformSignal, signalSpectrum, tot_sig.size(), this->errStr);
        if (!b) {
            throw std::runtime_error("Error: Signal FFT failed");
        }

        size_t naikwist = tot_sig.size() / 2 + 1;
        // на основании частоты Найквиста уменьшаем массив в два раза (из-за свойств БПФ)
        signalSpectrum.erase(signalSpectrum.begin() + naikwist, signalSpectrum.end());

        return signalSpectrum;
    }

    std::vector<double> ConcreteFourierSolver::inverseFastFourierTransform(std::vector<complex_type>& spec, double convCoef) {
        // To Do: перевод обратно в единицы измерения
        double scale = 1 / (sigParams.sampleInterval * convCoef); // To Do: check
        std::vector<complex_type> transformSpectrum(spec.size());
        std::transform(spec.begin(), spec.end(), transformSpectrum.begin(),
            [&scale](complex_type element) { return element *= scale; });

        // To Do: проверить количество нулей
        std::vector<complex_type> complexSignal;
        auto b = simple_fft::IFFT(transformSpectrum, complexSignal, spec.size(), this->errStr);
        if (!b) {
            throw std::runtime_error("Error: Spectrum IFFT failed");
        }
        // преобразование сигнала опять к реальным значениям 
        std::vector<double> transformSignal(complexSignal.size());
        std::transform(complexSignal.begin(), complexSignal.end(), transformSignal.begin(),
            [](complex_type element) { return element.real(); });
        return transformSignal;
    }

    const std::vector<double> FilterFourierSolver::getAmp() {
        return computeAmp(filteredSpectrum);
    }

    FilterFourierSolver::FilterFourierSolver(FourierSolver& fftSolver_, std::vector<double>& filter) : fftSolver(fftSolver_) {
        // To Do: проверка размеров массива и интерполяция характеристики

        // To Do: 
        double convCoef = 1;
        filterFreqResponse = fftSolver.fastFourierTransform(filter, convCoef);
    }

    std::vector<complex_type> FilterFourierSolver::fastFourierTransform(std::vector<double>& tot_sig, double convCoef) {
        return fftSolver.fastFourierTransform(tot_sig, convCoef);
    }

    std::vector<double> FilterFourierSolver::inverseFastFourierTransform(std::vector<complex_type>& filt_spec, double convCoef) {
        return fftSolver.inverseFastFourierTransform(filt_spec, convCoef);
    }

    void FilterFourierSolver::solve(std::vector<double>& tot_sig) {
        // To Do: не заполнены sigParams

        // расчет спектра для пушек
        fftSolver.solve(tot_sig);

        //spectrum = fftSolver.getSpectrum(); // (?)
        // расчет 

    }

    void FilterFourierSolver::applyFilter() {
        auto& spectrum = fftSolver.getSpectrum();
        if (spectrum.empty()) {
            throw std::runtime_error("Error: spectrum hasn't calculated yet");
        }
        if (spectrum.size() != filterFreqResponse.size()) {
            // To Do: интерполяция значений filterFreqResponse
        }
        filteredSpectrum.resize(spectrum.size());
        // перемножение спектров в частотной области
        for (size_t i = 0; i < spectrum.size(); i++) {
            filteredSpectrum[i] = spectrum[i] * filterFreqResponse[i];
        }
        // восстановление сигнала после обратного Фурье преобразования (размер массива в 2 меньше tot_sig еще из-за прямого преобразования Фурье)
        filteredSignal = inverseFastFourierTransform(filteredSpectrum); // To Do: проверить единицы измерения
    }

} // namespace fourier_solver

