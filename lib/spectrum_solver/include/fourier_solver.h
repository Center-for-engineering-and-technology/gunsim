#pragma once

#include <simple_fft/fft_settings.h>

#include <vector>

#include <gund_structs.h>

//const std::vector<complex_type> a;
//std::vector<complex_type> b;
//const size_t size, const char*& error_description;
//auto b = simple_fft::FFT(a, b, size, error_description);

namespace fourier_solver {

    // расчет амплитуды спектра на основании преобразования Фурье (дБ в логарифмическом масштабе) - внешнезаданный спектр
    std::vector<double> computeAmp(std::vector<complex_type> spectrum);

    class FourierSolver {
    public:
        virtual void solve(std::vector<double>& totalSignal) = 0;

        virtual const std::vector<double>& getSignal() = 0;
        virtual const std::vector<complex_type>& getSpectrum() = 0;

        virtual const std::vector<double> getAmp() = 0;

        virtual std::vector<complex_type> fastFourierTransform(std::vector<double>& totalSignal, double convCoef = 10e10) = 0;
        virtual std::vector<double> inverseFastFourierTransform(std::vector<complex_type>& spectrumSignal, double convCoef = 10e10) = 0;
    protected:
        

    };

    class ConcreteFourierSolver : public FourierSolver {
    public:
    //    ConcreteFourierSolver() = default;
        ConcreteFourierSolver(gund_structs::SignatureParameters& sigParams_)
            : sigParams(sigParams_) {}

        void solve(std::vector<double>& totalSignal) override; // To Do: не заполнены sigParams
    //    // функция для тестирования
    //    // To Do: fft(signal) -> spectrum -> ifft(spectrum) -> signal' -> сравнить
    //    // помнить, что размер массивов уменьшается в два раза при преобразовании Фурье
    //    // bool testFFT();

        const std::vector<double>& getSignal() override { return signal; }
        const std::vector<complex_type>& getSpectrum() override { return spectrum; }
        const gund_structs::SignatureParameters& getSpecParams() { return specParams; }
        // рассчитывает амплитуду спектра исходного сигнала
        const std::vector<double> getAmp() override;
        // расчет фазы спектра на основании преобразования Фурье (мб объединить с computeAmp)
        void computePhase();

        std::vector<complex_type> fastFourierTransform(std::vector<double>& totalSignal, double convCoef = 10e10) override; // To Do: Проверить единицы измерения
        std::vector<double> inverseFastFourierTransform(std::vector<complex_type>& spectrumSignal, double convCoef = 10e10) override; // To Do: Проверить перевод единиц измерения
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

    // decorator
    class FilterFourierSolver : public FourierSolver {
    public:
    //    FilterFourierParams() = default;
        FilterFourierSolver(FourierSolver& fftSolver_, std::vector<double>& filter);

        // To Do: умножение частотной характеристики фильтра на спектр сигнала от ПИ (с нуля)
        void solve(std::vector<double>& totalSignal) override;
        // добавление фильтра к сигналу из fftSolver
        void applyFilter();

        // возвращает отфильтрованный сигнал
        const std::vector<double>& getSignal() override { return filteredSignal; }
        // возвращает спектр отфильтрованного сигнала
        const std::vector<complex_type>& getSpectrum() override { return filteredSpectrum; }
        // рассчитывает амплитуду спектра фильтрованного сигнала
        const std::vector<double> getAmp() override;

        // преобразование Фурье "с нуля" - для заполнениЯ параметров в fftSolver
        std::vector<complex_type> fastFourierTransform(std::vector<double>& totalSignal, double convCoef = 10e10) override; // To Do: check convCoef
        // обратное преобразование Фурье
        std::vector<double> inverseFastFourierTransform(std::vector<complex_type>& filteredSpectrumSignal, double convCoef = 10e10) override; // To Do: Проверить перевод единиц измерения
    protected:
        FourierSolver& fftSolver;
        // частотная характеристика фильтра
        std::vector<complex_type> filterFreqResponse;
        // отфильтрованный сигнал
        std::vector<double> filteredSignal;
        // отфильтрованный спектр
        std::vector<complex_type> filteredSpectrum;
    };


} // namespace fourier_solver
