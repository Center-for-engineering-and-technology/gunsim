#include <comparative_analysis.h>

#include <common.h>

namespace comparative_analysis
{

namespace
{
double Bar2dB (double bar_value)
{
    // Добавляем 220 т.к. приводим из Бар*м к мкПа*м
    // 20 * log10(10^11) = 220
    return 20. * std::log10(bar_value) + 220;
}
} // namespace

// TODO: это, по факту, дублирование кода из spectumprocessor. Перенести все связанное с фурье и спектрами в common!
std::vector<double> computeAmp (const std::vector<complex_type>& spectrum)
{
    if (spectrum.empty())
    {
        return {};
    }

    std::vector<double> ampSpec(spectrum.size());
    for (size_t i = 0; i < spectrum.size(); i++)
    {
        ampSpec[i] = 20 * std::log10(std::abs(spectrum[i]));
    }
    return ampSpec;
}

// TODO: это, по факту, дублирование кода из spectumprocessor. Перенести все связанное с фурье и спектрами в common!
std::vector<double> computePhase (const std::vector<complex_type>& spectrum)
{
    if (spectrum.empty())
    {
        return {};
    }

    std::vector<double> phaseSpec(spectrum.size());
    std::vector<double> unwrapper(spectrum.size());
    for (size_t i = 0; i < spectrum.size(); i++)
    {
        phaseSpec[i] = std::atan2(spectrum[i].imag(), spectrum[i].real());
        if (i > 0)
        {
            unwrapper[i] = unwrapper[i - 1];
            if (std::abs(spectrum[i]) > 1e-3)
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
    for (size_t i = 0; i < spectrum.size(); i++)
    {
        phaseSpec[i] += unwrapper[i];
        phaseSpec[i] *= -1;
    }
    return phaseSpec;
}

void OuterSignatureAnalysis::parseSignature (std::string inputFileName)
{
    gund_format_parser::GundalfOutputParser parser(dataPath + '/' + inputFileName);
    signature = parser.getData();
    std::transform(signature.begin(), signature.end(), signature.begin(), [] (double x)
                   {
                       return x / 1e5;
                   }); // переводим [Pa*m] в [bar*m]
    auto [min, max] = std::minmax_element(signature.begin(), signature.end());
    if (std::distance(signature.begin(), min) < std::distance(signature.begin(), max))
        std::transform(signature.begin(), signature.end(), signature.begin(), [] (double x)
                       {
                           return -x;
                       }); // разворот оси в случае перевернутой сигнатуры
    auto sigParams = parser.getSigParams();
    N = sigParams.sampleNum;
    dt = sigParams.sampleInterval;
    times.resize(signature.size());
    for (size_t i = 0; i < times.size(); ++i)
        times[i] = i * dt;
}

void OuterSignatureAnalysis::computations ()
{
    gund_structs::SignatureParameters params;
    params.sampleNum = N;
    params.sampleInterval = dt;
    params.sampleMax = N * dt;
    fourier_solver::ConcreteFourierSolver solver(params);
    solver.solve(signature);
    auto specParams = solver.getSpecParams();

    amplitude_spectrum = solver.getAmp();
    phase_spectrum = solver.getPhase();
    peak_to_peak = peakToPeakCompute();
    db_peak_to_peak = Bar2dB(peak_to_peak);
    zero_to_peak = zeroToPeakCompte();
    db_zero_to_peak = Bar2dB(zero_to_peak);
    mean_square_pressure = meanSquarePressureCompute();
    primary_to_bubble = primaryToBubbleCompute();
    bubble_period = bubblePeriodCompute();
    max_spectral_value = maxSpectralValueCompute(specParams);
    average_spectral_value = averageSpectralValueCompute(specParams);
    max_spectral_ripple = maxSpectralRippleCompute(specParams);

    frequencies.resize(amplitude_spectrum.size());
    for (size_t i = 0; i < frequencies.size(); ++i)
        frequencies[i] = i * specParams.sampleInterval;
}

double OuterSignatureAnalysis::peakToPeakCompute () const
{
    auto [min, max] = std::minmax_element(signature.begin(), signature.end());
    return *max - *min;
}

double OuterSignatureAnalysis::zeroToPeakCompte () const
{
    auto max = std::max_element(signature.begin(), signature.end());
    return *max;
}

double OuterSignatureAnalysis::meanSquarePressureCompute () const
{
    return std::sqrt(std::accumulate(signature.begin(), signature.end(), 0, [] (double sum, double i)
                                     {
                                         return sum + i * i;
                                     })
                     / static_cast<double>(N));
}

double OuterSignatureAnalysis::primaryToBubbleCompute () const
{
    auto [min, max] = std::minmax_element(signature.begin(), signature.end());
    auto max_bubble = std::max_element(min, signature.end());
    auto min_bubble = std::min_element(max_bubble, signature.end());
    return (*max - *min) / (*max_bubble - *min_bubble);
}

double OuterSignatureAnalysis::bubblePeriodCompute () const
{
    auto min = std::min_element(signature.begin(), signature.end());
    auto max_bubble = std::max_element(min, signature.end());
    return std::distance(signature.begin(), max_bubble) * dt;
}

double OuterSignatureAnalysis::maxSpectralValueCompute (const gund_structs::SignatureParameters& specParams) const
{
    size_t start_box_index = static_cast<size_t>(10. / specParams.sampleInterval);
    size_t end_box_index = static_cast<size_t>(50. / specParams.sampleInterval);
    auto max_box_spec = std::max_element(amplitude_spectrum.begin() + start_box_index, amplitude_spectrum.begin() + end_box_index + 1);
    return *max_box_spec;
}

double OuterSignatureAnalysis::averageSpectralValueCompute (const gund_structs::SignatureParameters& specParams) const
{
    size_t start_box_index = static_cast<size_t>(10. / specParams.sampleInterval);
    size_t end_box_index = static_cast<size_t>(50. / specParams.sampleInterval);
    return std::accumulate(amplitude_spectrum.begin() + start_box_index, amplitude_spectrum.begin() + end_box_index + 1, 0.) / (end_box_index + 1 - start_box_index);
}

double OuterSignatureAnalysis::maxSpectralRippleCompute (const gund_structs::SignatureParameters& specParams) const
{
    size_t start_box_index = static_cast<size_t>(10. / specParams.sampleInterval);
    size_t end_box_index = static_cast<size_t>(50. / specParams.sampleInterval);
    auto [min_box_spec, max_box_spec] = std::minmax_element(amplitude_spectrum.begin() + start_box_index, amplitude_spectrum.begin() + end_box_index + 1);
    return *max_box_spec - *min_box_spec;
}

void OuterSignatureAnalysis::printResult (std::string outputFileName) const
{
    nlohmann::json jsonObject;
    if (amplitude_spectrum.size() != frequencies.size())
        throw std::logic_error("wrong size amplitude");
    if (phase_spectrum.size() != frequencies.size())
        throw std::logic_error("wrong size phase");
    if (signature.size() != times.size())
        throw std::logic_error("wrong size signature");

    jsonObject["signature"]["y-axis (bar*m)"] = nlohmann::json(signature);
    jsonObject["signature"]["x-axis (s)"] = nlohmann::json(times);
    jsonObject["amplitude_spectrum"]["y-axis (dB)"] = nlohmann::json(amplitude_spectrum);
    jsonObject["amplitude_spectrum"]["x-axis (Hz)"] = nlohmann::json(frequencies);
    jsonObject["phase_spectrum"]["y-axis (radian)"] = nlohmann::json(phase_spectrum);
    jsonObject["phase_spectrum"]["x-axis (Hz)"] = nlohmann::json(frequencies);
    jsonObject["peak_to_peak"] = peak_to_peak;
    jsonObject["db_peak_to_peak"] = db_peak_to_peak;
    jsonObject["zero_to_peak"] = zero_to_peak;
    jsonObject["db_zero_to_peak"] = db_zero_to_peak;
    jsonObject["mean_square_pressure"] = mean_square_pressure;
    jsonObject["primary_to_bubble"] = primary_to_bubble;
    jsonObject["bubble_period"] = bubble_period;
    jsonObject["max_spectral_value"] = max_spectral_value;
    jsonObject["average_spectral_value"] = average_spectral_value;
    jsonObject["max_spectral_ripple"] = max_spectral_ripple;

    std::ofstream outputFile(dataPath + '/' + outputFileName);
    outputFile << std::setw(4) << std::setprecision(17) << jsonObject;
    outputFile.close();
}

SignatureMetricsHolder ComputeMetrics (
    const gund_structs::NewSignatureParameters& params,
    const std::vector<double>& signature
)
{
    SignatureMetricsHolder result;

    if (signature.empty())
    {
        return result;
    }

    auto dt = params.sample_interval;

    auto [min, max] = std::minmax_element(signature.begin(), signature.end());

    result.zero_to_peak = *max;
    result.db_peak_to_peak = Bar2dB(result.peak_to_peak);
    result.peak_to_peak = *max - *min;
    result.db_zero_to_peak = Bar2dB(result.zero_to_peak);

    auto max_bubble = std::max_element(min, signature.end());
    auto min_bubble = std::min_element(max_bubble, signature.end());

    result.primary_to_bubble = (*max - *min) / (*max_bubble - *min_bubble);
    result.bubble_period = std::distance(signature.begin(), max_bubble) * dt;

    result.mean_square_pressure = std::sqrt(std::accumulate(signature.begin(), signature.end(), 0, [] (double sum, double i)
                                                            {
                                                                return sum + i * i;
                                                            })
                                            / static_cast<double>(signature.size()));

    auto spectrum = filters_common::RFFTReal(signature);
    auto amp_spec = computeAmp(spectrum);

    // Рассматриваем полосу от 10 до 50 герц
    size_t start_box_index = static_cast<size_t>(10. * dt);
    size_t end_box_index = static_cast<size_t>(50. * dt);

    if (start_box_index < amp_spec.size() && end_box_index < amp_spec.size())
    {
        result.max_spectral_value = *std::max_element(
            amp_spec.begin() + start_box_index,
            amp_spec.begin() + end_box_index + 1
        );

        result.average_spectral_value = std::accumulate(
                                            amp_spec.begin() + start_box_index,
                                            amp_spec.begin() + end_box_index + 1,
                                            0.
                                        )
                                      / (end_box_index + 1 - start_box_index);

        auto [min_box_spec, max_box_spec] = std::minmax_element(
            amp_spec.begin() + start_box_index,
            amp_spec.begin() + end_box_index + 1
        );

        result.max_spectral_ripple = *max_box_spec - *min_box_spec;
    }

    return result;
}
}; // namespace comparative_analysis
