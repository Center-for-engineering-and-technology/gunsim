#pragma once

#include <fstream>

#include <nlohmann/json.hpp>

#include <signature_structs.h>
#include <gund_format_parser.h>
#include <fourier_solver.h>

namespace comparative_analysis
{

class OuterSignatureAnalysis
{
private:

    std::string dataPath;          // папка с файлом сигнатуры
    double dt;                     // шаг по времени
    size_t N;                      // длина сигнатуры
    std::vector<double> times;     // сетка по времени
    std::vector<double> signature; // внешняя сигнатура

    /* расчетные величины */
    std::vector<double> frequencies;
    std::vector<double> amplitude_spectrum;
    std::vector<double> phase_spectrum;

    // TODO: есть дублирование полей с метриками и их расчетом,
    // убрать отсюда и заменить на SignatureMetricsHolder
    double peak_to_peak;
    double db_peak_to_peak;
    double zero_to_peak;
    double db_zero_to_peak;
    double mean_square_pressure;
    double primary_to_bubble;
    double bubble_period;
    double max_spectral_value;
    double average_spectral_value;
    double max_spectral_ripple;

    double peakToPeakCompute() const;
    double zeroToPeakCompte() const;
    double meanSquarePressureCompute() const;
    double primaryToBubbleCompute() const;
    double bubblePeriodCompute() const;
    double maxSpectralValueCompute(const gund_structs::SignatureParameters& specParams) const;
    double averageSpectralValueCompute(const gund_structs::SignatureParameters& specParams) const;
    double maxSpectralRippleCompute(const gund_structs::SignatureParameters& specParams) const;

public:

    OuterSignatureAnalysis() = default;
    void setDataPath (std::string dataPath_)
    {
        dataPath = dataPath_;
    }
    void parseSignature(std::string inputFileName);
    void computations();
    void printResult(std::string outputFileName) const;
};

struct SignatureMetricsHolder
{
    double peak_to_peak{0};
    double db_peak_to_peak{0};
    double zero_to_peak{0};
    double db_zero_to_peak{0};
    double mean_square_pressure{0};
    double primary_to_bubble{0};
    double bubble_period{0};

    double max_spectral_value{0};
    double average_spectral_value{0};
    double max_spectral_ripple{0};
};

template<typename T, typename = gund_structs::enable_if_floating_or_complex_t<T>>
SignatureMetricsHolder ComputeMetrics (const gund_structs::Signature<T>& signature)
{
    if constexpr (gund_structs::is_complex_of_floating_v<T>)
    {
        throw std::logic_error("Metrics computation for signatures which lives in frequency domain not yet implemented");
    }
    return ComputeMetrics(signature.params, signature.data);
}

SignatureMetricsHolder ComputeMetrics(
    const gund_structs::NewSignatureParameters& params,
    const std::vector<double>& signature
);

SignatureMetricsHolder ComputeMetrics(
    const gund_structs::NewSignatureParameters& params,
    const std::vector<std::complex<double>>& signature
);

}; // namespace comparative_analysis
