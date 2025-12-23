#pragma once

#include <simple_fft/fft.h>

#include <vector>

#include <gund_structs.h>
#include <gund_utility.h>

namespace fourier_solver
{

// расчет амплитуды спектра на основании преобразования Фурье (дБ в логарифмическом масштабе) - внешнезаданный спектр
std::vector<double> computeAmp(std::vector<complex_type> spectrum);
// расчет фазы спектра на основании преобразования Фурье (радианы) - внешнезаданный спектр
std::vector<double> computePhase(std::vector<complex_type> spectrum);

class FourierSolver
{
public:

    virtual void solve(std::vector<double>& totalSignal) = 0;

    virtual const std::vector<double>& getSignal() = 0;
    virtual const std::vector<complex_type>& getSpectrum() = 0;

    virtual const std::vector<double> getAmp() = 0;
    virtual const std::vector<double> getPhase() = 0;

    virtual std::vector<complex_type> fastFourierTransform(std::vector<double>& totalSignal, double convCoef = 10e10) = 0;
    virtual std::vector<double> inverseFastFourierTransform(std::vector<complex_type>& spectrumSignal, double convCoef = 10e10) = 0;

protected:
};

class ConcreteFourierSolver : public FourierSolver
{
public:

    //    ConcreteFourierSolver() = default;
    ConcreteFourierSolver (gund_structs::SignatureParameters& sigParams_)
        : sigParams(sigParams_)
    {
    }

    void solve(std::vector<double>& totalSignal) override;
    //    // функция для тестирования
    //    // To Do: fft(signal) -> spectrum -> ifft(spectrum) -> signal' -> сравнить
    //    // помнить, что размер массивов уменьшается в два раза при преобразовании Фурье
    //    // bool testFFT();

    const std::vector<double>& getSignal () override
    {
        return signal;
    }
    const std::vector<complex_type>& getSpectrum () override
    {
        return spectrum;
    }
    const gund_structs::SignatureParameters& getSpecParams ()
    {
        return specParams;
    }
    // рассчитывает амплитуду спектра исходного сигнала
    const std::vector<double> getAmp() override;
    // рассчитывает фазу спектра исходного сигнала
    const std::vector<double> getPhase() override;

    std::vector<complex_type> fastFourierTransform(std::vector<double>& totalSignal, double convCoef = 10e10) override;
    std::vector<double> inverseFastFourierTransform(std::vector<complex_type>& spectrumSignal, double convCoef = 10e10) override;

protected:

    // сигнал ПИ
    std::vector<double> signal;
    // спектр сигнала после преобразований Фурье
    std::vector<complex_type> spectrum;

    // параметры модели
    gund_structs::SignatureParameters sigParams;
    gund_structs::SignatureParameters specParams;

private:

    const char* errStr = NULL;
};

//// To Do: функция для интерполяции значений

// функции находятся в public для последовательного использования нескольких фильтров

// decorator
class FilterFourierSolver : public FourierSolver
{
public:

    //    FilterFourierParams() = default;
    FilterFourierSolver(FourierSolver& fftSolver_, gund_structs::SignatureParameters& sigParams_, std::vector<double>& filter);

    // умножение частотной характеристики фильтра на спектр сигнала от ПИ (с нуля, если не заполнены fftSolver)
    void solve(std::vector<double>& totalSignal) override;
    // добавление фильтра к сигналу из fftSolver
    // можно использовать для самостоятельного вызова (при заполненных fftSolver)
    // To Do: сохранить способ фильтрации в частотной области как back, во временной - как forward
    void applyFilter();

    // возвращает отфильтрованный сигнал
    const std::vector<double>& getSignal () override
    {
        return filteredSignal;
    }
    // возвращает спектр отфильтрованного сигнала
    const std::vector<complex_type>& getSpectrum () override
    {
        return filteredSpectrum;
    }
    // рассчитывает амплитуду спектра фильтрованного сигнала
    const std::vector<double> getAmp() override;
    // рассчитывает фазу спектра исходного сигнала
    const std::vector<double> getPhase() override;

    // свертка функции во временной области (для сохранения размера массива)
    std::vector<double> convolveSignal(const std::vector<double>& totalSignal);
    // преобразование Фурье "с нуля" - для заполнениЯ параметров в fftSolver
    std::vector<complex_type> fastFourierTransform(std::vector<double>& totalSignal, double convCoef = 10e10) override;
    // обратное преобразование Фурье
    std::vector<double> inverseFastFourierTransform(std::vector<complex_type>& filteredSpectrumSignal, double convCoef = 10e10) override;

protected:

    FourierSolver& fftSolver;
    // временная характеристика фильтра
    std::vector<double> filterResponse;
    // частотная характеристика фильтра
    std::vector<complex_type> filterFreqResponse;
    // отфильтрованный сигнал
    std::vector<double> filteredSignal;
    // отфильтрованный спектр
    std::vector<complex_type> filteredSpectrum;
    // параметры модели
    gund_structs::SignatureParameters filterSigParams;
};

} // namespace fourier_solver
