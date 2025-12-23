#include <fourier_solver.h>

#include <simple_fft/fft.h>

#include <algorithm>

namespace fourier_solver
{

std::vector<double> computeAmp (std::vector<complex_type> spectrum_)
{
    if (spectrum_.empty())
    {
        throw std::runtime_error("Error: Spectrum data is empty");
    }

    std::vector<double> ampSpec(spectrum_.size());
    for (size_t i = 0; i < spectrum_.size(); i++)
    {
        ampSpec[i] = 20 * std::log10(std::abs(spectrum_[i]));
    }
    return ampSpec;
}

std::vector<double> computePhase (std::vector<complex_type> spectrum_)
{
    if (spectrum_.empty())
    {
        throw std::runtime_error("Error: Spectrum data is empty");
    }

    std::vector<double> phaseSpec(spectrum_.size()), unwrapper(spectrum_.size());
    for (size_t i = 0; i < spectrum_.size(); i++)
    {
        phaseSpec[i] = std::atan2(spectrum_[i].imag(), spectrum_[i].real());
        if (i > 0)
        {
            unwrapper[i] = unwrapper[i - 1];
            if (std::abs(spectrum_[i]) > 1e-3)
            {
                if (phaseSpec[i] > M_PI / 2. && phaseSpec[i - 1] < -M_PI / 2.)
                    unwrapper[i] -= 2 * M_PI;
                if (phaseSpec[i] < -M_PI / 2. && phaseSpec[i - 1] > M_PI / 2.)
                    unwrapper[i] += 2 * M_PI;
            }
        }
        else
        {
            unwrapper[i] = 0;
        }
    }
    for (size_t i = 0; i < spectrum_.size(); i++)
    {
        phaseSpec[i] += unwrapper[i];
        phaseSpec[i] *= -1;
    }
    return phaseSpec;
}

const std::vector<double> ConcreteFourierSolver::getAmp ()
{
    return computeAmp(spectrum);
}

const std::vector<double> ConcreteFourierSolver::getPhase ()
{
    return computePhase(spectrum);
}

void ConcreteFourierSolver::solve (std::vector<double>& tot_sig)
{
    signal = tot_sig;

    auto newSigSize = gund_utility::findNextPowerOfTwo(tot_sig.size());
    // для эффективного использования БПФ требуется количество элементов, кратное степени двойки
    // добавляем необходимое количество нулевых элементов в конец массива
    // UPD: добавить не в конец, а в начало
    std::vector<double> updSignal(newSigSize, 0);
    std::copy(tot_sig.begin(), tot_sig.end(), updSignal.begin());

    // расчет спектра для сигнала с добавочными значениями
    spectrum = fastFourierTransform(updSignal);

    specParams.sampleNum = static_cast<double>(newSigSize);
    specParams.sampleInterval = 1 / (sigParams.sampleInterval * newSigSize);
    specParams.sampleMax = specParams.sampleInterval * specParams.sampleNum;
}

std::vector<complex_type> ConcreteFourierSolver::fastFourierTransform (std::vector<double>& tot_sig, double convCoef)
{
    // переводим в мкПа
    double scale = sigParams.sampleInterval * convCoef;
    std::vector<double> transformSignal(tot_sig.size());
    std::transform(tot_sig.begin(), tot_sig.end(), transformSignal.begin(), [&scale] (double element)
                   {
                       return element *= scale;
                   });

    std::vector<complex_type> signalSpectrum(tot_sig.size());
    auto b = simple_fft::FFT(transformSignal, signalSpectrum, tot_sig.size(), this->errStr);
    if (!b)
    {
        throw std::runtime_error("Error: Signal FFT failed");
    }

    size_t naikwist = tot_sig.size() / 2 + 1;
    // на основании частоты Найквиста уменьшаем массив в два раза (из-за свойств БПФ)
    signalSpectrum.erase(signalSpectrum.begin() + naikwist, signalSpectrum.end());

    return signalSpectrum;
}

std::vector<double> ConcreteFourierSolver::inverseFastFourierTransform (std::vector<complex_type>& spec, double convCoef)
{
    // переводим обратно в Па
    double scale = 1 / (sigParams.sampleInterval * convCoef);
    std::vector<complex_type> transformSpectrum(spec.size());
    std::transform(spec.begin(), spec.end(), transformSpectrum.begin(), [&scale] (complex_type element)
                   {
                       return element *= scale;
                   });

    // To Do: проверить количество нулей
    std::vector<complex_type> complexSignal;
    complexSignal.resize(spec.size());
    auto b = simple_fft::IFFT(transformSpectrum, complexSignal, spec.size(), this->errStr);
    if (!b)
    {
        throw std::runtime_error("Error: Spectrum IFFT failed");
    }
    // преобразование сигнала опять к реальным значениям
    std::vector<double> transformSignal(complexSignal.size());
    std::transform(complexSignal.begin(), complexSignal.end(), transformSignal.begin(), [] (complex_type element)
                   {
                       return element.real();
                   });
    return transformSignal;
}

const std::vector<double> FilterFourierSolver::getAmp ()
{
    return computeAmp(filteredSpectrum);
}

const std::vector<double> FilterFourierSolver::getPhase ()
{
    return computePhase(filteredSpectrum);
}

FilterFourierSolver::FilterFourierSolver (
    FourierSolver& fftSolver_,
    gund_structs::SignatureParameters& sigParams_,
    std::vector<double>& filter
)
    : fftSolver(fftSolver_), filterSigParams(sigParams_)
{
    // проверка размеров массива и интерполяция характеристики остается на ответственности пользователя
    filterResponse = filter;
    // считаем частотную характеристику фильтра
    double convCoef = 1;
    auto newFiltSize = gund_utility::findNextPowerOfTwo(filter.size());
    // обновление размера фильтра для БПФ
    std::vector<double> updFilter(newFiltSize, 0.);
    std::copy(filter.begin(), filter.end(), updFilter.begin());
    filterFreqResponse = fftSolver.fastFourierTransform(updFilter, convCoef);
}

std::vector<complex_type> FilterFourierSolver::fastFourierTransform (std::vector<double>& tot_sig, double convCoef)
{
    return fftSolver.fastFourierTransform(tot_sig, convCoef);
}

std::vector<double> FilterFourierSolver::inverseFastFourierTransform (std::vector<complex_type>& filt_spec, double convCoef)
{
    return fftSolver.inverseFastFourierTransform(filt_spec, convCoef);
}

void FilterFourierSolver::solve (std::vector<double>& tot_sig)
{
    auto& spec = fftSolver.getSpectrum();
    if (spec.empty())
    {
        // расчет спектра для пушек
        fftSolver.solve(tot_sig);
    }
    applyFilter();
}

void FilterFourierSolver::applyFilter ()
{
    // способ 1: перемножение сигнала в частотной области
    /*auto& spectrum = fftSolver.getSpectrum();
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
    */
    // способ 2: свертка сигнала во временной области
    auto& signal = fftSolver.getSignal();
    filteredSignal = convolveSignal(signal);
    // нормировка массива
    // нормировка массива проводится после свертки для хранения исходной характеристики фильтра (для уменьшения мест возникновения ошибок при некорректной нормировке)
    // аналогичный результат будет достигнут при нормровке характеристики до свертки
    auto max = std::max_element(filterResponse.begin(), filterResponse.end());
    if (*max > 10.)
    {
        double scale = filterSigParams.sampleInterval;
        std::transform(filteredSignal.begin(), filteredSignal.end(), filteredSignal.begin(), [&scale] (double element)
                       {
                           return element *= scale;
                       });
    }
    // для эффективного использования БПФ требуется количество элементов, кратное степени двойки
    auto newSigSize = gund_utility::findNextPowerOfTwo(filteredSignal.size());
    std::vector<double> updSignal(newSigSize, 0.);
    std::copy(filteredSignal.begin(), filteredSignal.end(), updSignal.begin());

    // расчет спектра для сигнала с добавочными значениями
    filteredSpectrum = fastFourierTransform(updSignal);
}

std::vector<double> FilterFourierSolver::convolveSignal (const std::vector<double>& totalSignal)
{
    // в качестве свертки используем модель полной свертки сигнала
    // для соответствия размеру массива totalSignal ограничим итоговый массив его размером
    size_t signalSize = totalSignal.size();
    size_t filterSize = filterResponse.size();
    std::vector<double> convResult(signalSize + filterSize - 1, 0.);
    // процесс свертки
    for (size_t i = 0; i < convResult.size(); i++)
    {
        for (size_t j = 0; j < filterSize; j++)
        {
            if ((i - j) >= 0 && (i - j) < signalSize)
            {
                convResult[i] += totalSignal[i - j] * filterResponse[j];
            }
        }
    }
    return std::vector<double>{convResult.begin(), convResult.begin() + signalSize};
}

} // namespace fourier_solver
