#include <internal_filter.h>

namespace internal_filter
{

void AbstractInternalFilter::DFT::resize (size_t size_)
{
    size = size_;
    exp.resize(size);
    for (size_t k = 0; k < size; ++k)
        exp[k] = std::exp(-2 * M_PI * I / static_cast<double>(size) * static_cast<double>(k));
}

template<typename T>
void AbstractInternalFilter::DFT::apply (const std::vector<T>& signal, std::vector<std::complex<double>>& spectrum)
{

    if (size != signal.size())
        resize(signal.size());

    spectrum.resize(size);
    for (size_t k = 0; k < size; ++k)
    {
        spectrum[k] = 0;
        for (size_t n = 0; n < size; ++n)
            spectrum[k] += signal[n] * exp[(k * n) % size];
    }
}

template<typename T>
void AbstractInternalFilter::DFT::applyInverse (const std::vector<T>& spectrum, std::vector<std::complex<double>>& signal)
{
    if (size != spectrum.size())
        resize(spectrum.size());

    signal.resize(size);
    for (size_t n = 0; n < size; ++n)
    {
        signal[n] = 0;
        for (size_t k = 0; k < size; ++k)
            signal[n] += spectrum[k] / exp[(k * n) % size];
        signal[n] /= static_cast<double>(size);
    }
}

AbstractInternalFilter::Window::Window (WindowType window_type_)
    : window_type(window_type_)
{
    switch (window_type)
    {
        case WindowType::RECTANGULAR:
            w = [] (double n, double m) -> double
            {
                if (std::abs(n) > (m - 1) / 2.)
                    return 0;
                return 1;
            };
            peak_stopband_attenuation = -21;
            transition_bandwith_index = 0.4069;
            break;
        case WindowType::HAMMING:
            w = [] (double n, double m) -> double
            {
                if (std::abs(n) > (m - 1) / 2.)
                    return 0;
                return 0.54 + 0.46 * std::cos(2 * M_PI * n / static_cast<double>(m - 1));
            };
            peak_stopband_attenuation = -44;
            transition_bandwith_index = 0.1470;
            break;
        case WindowType::HANN:
            w = [] (double n, double m) -> double
            {
                if (std::abs(n) > (m - 1) / 2.)
                    return 0;
                return 0.5 + 0.5 * std::cos(2 * M_PI * n / static_cast<double>(m - 1));
            };
            peak_stopband_attenuation = -54;
            transition_bandwith_index = 0.1398;
            break;
        case WindowType::BLACKMAN:
            w = [] (double n, double m) -> double
            {
                if (std::abs(n) > (m - 1) / 2.)
                    return 0;
                return 0.42 + 0.5 * std::cos(2 * M_PI * n / static_cast<double>(m - 1)) + 0.08 * std::cos(4 * M_PI * n / static_cast<double>(m - 1));
            };
            peak_stopband_attenuation = -75;
            transition_bandwith_index = 0.0704;
            break;
    }
}

void AbstractInternalFilter::windowedSincResponce (std::vector<double>& impulse) const
{

    double low_cut_off = M_PI * low_pass / N;
    double high_cut_off = M_PI * high_pass / N;
    double low_transition_bandwidth = 2 * low_cut_off * (std::pow(2, -low_window.peak_stopband_attenuation / low_slope) - 1.) / (std::pow(2, -low_window.peak_stopband_attenuation / low_slope) + 1.);
    double high_transition_bandwidth = 2 * high_cut_off * (std::pow(2, -high_window.peak_stopband_attenuation / high_slope) - 1.) / (std::pow(2, -high_window.peak_stopband_attenuation / high_slope) + 1.);
    int low_m = 2 * static_cast<int>(std::ceil(M_PI / 2. / low_window.transition_bandwith_index / low_transition_bandwidth));
    int high_m = 2 * static_cast<int>(std::ceil(M_PI / 2. / high_window.transition_bandwith_index / high_transition_bandwidth));
    int m = std::max(low_m, high_m);

    std::function<double(double)> sinc = [] (double x) -> double
    {
        if (std::abs(x) < 1e-6)
            return 1 - std::pow(x, 2) / 6. + std::pow(x, 4) / 120. - std::pow(x, 6) / 5040.;
        return std::sin(x) / x;
    };

    impulse.resize(m);
    std::vector<double> low_imp_data(m), high_imp_data(m);
    double low_sum = 0, high_sum = 0;
    for (int i = 0; i < m; ++i)
    {
        low_imp_data[i] = low_window.w(i - (m - 1) / 2., low_m) * sinc(low_cut_off * (i - (m - 1) / 2.)) * low_cut_off;
        high_imp_data[i] = high_window.w(i - (m - 1) / 2., high_m) * sinc(high_cut_off * (i - (m - 1) / 2.)) * high_cut_off;
        low_sum += low_imp_data[i];
        high_sum += high_imp_data[i];
    }

    for (int i = 0; i < m; ++i)
        impulse[i] = high_imp_data[i] / high_sum - low_imp_data[i] / low_sum;
}

void AbstractInternalFilter::printInFile (std::string dataPath) const
{
    std::ofstream output;
    output.open(dataPath + filtername + ".flt");
    output << "# dt = " << dt << "\n";
    output << "# iz = " << iz << "\n";
    output << "# ns = " << imp_data.size() << "\n";
    for (size_t i = 0; i < imp_data.size(); ++i)
        output << imp_data[i] << "\n";
    output.close();
}

void InternalZeroPhaseFilter::computeImpulseResponce ()
{

    windowedSincResponce(imp_data);
    DFT dft;
    std::vector<std::complex<double>> freq_data, tmp_data;
    dft.apply(imp_data, freq_data);
    double delta = std::pow(10, high_window.peak_stopband_attenuation / 20.);
    size_t m = imp_data.size();
    for (size_t i = 0; i < freq_data.size(); ++i)
    {
        double s = (freq_data[i] * std::exp(M_PI * I * static_cast<double>(i * (m - 1)) / static_cast<double>(m))).real();
        freq_data[i] = (s + delta) / (1. + delta) * std::exp(-M_PI * I * static_cast<double>(i * (m - 1)) / static_cast<double>(m));
    }
    dft.applyInverse(freq_data, tmp_data);
    for (size_t i = 0; i < imp_data.size(); ++i)
        imp_data[i] = tmp_data[i].real();
    iz = static_cast<int>(imp_data.size()) / 2;
}

void InternalZeroPhaseFilter::apply (const std::vector<double>& signal, std::vector<double>& filtered_signal) const
{

    size_t m = imp_data.size();
    filtered_signal.resize(signal.size());
    for (size_t i = 0; i < filtered_signal.size(); ++i)
    {
        filtered_signal[i] = 0;
        for (size_t k = 0; k < m; ++k)
        {
            double s = 0;
            if (i + m / 2 - k >= 0 && i + m / 2 - k < signal.size())
                s = signal[i + m / 2 - k];
            filtered_signal[i] += s * imp_data[k];
        }
    }
}

void InternalMinimalPhaseFilter::computeImpulseResponce ()
{

    windowedSincResponce(imp_data);
    DFT dft;
    std::vector<std::complex<double>> freq_data, tmp_data, phase_data;
    dft.apply(imp_data, freq_data);
    double delta = std::pow(10, high_window.peak_stopband_attenuation / 20.);
    for (size_t i = 0; i < freq_data.size(); ++i)
        freq_data[i] = std::log((std::abs(freq_data[i]) + delta) / (delta + 1.));
    dft.apply(freq_data, tmp_data);
    for (size_t i = 0; i < tmp_data.size(); ++i)
        if (i == 0 || i == tmp_data.size() / 2)
            tmp_data[i] = 0;
        else
        {
            if (i < tmp_data.size() / 2)
                tmp_data[i] *= I;
            if (i > tmp_data.size() / 2)
                tmp_data[i] *= -I;
        }
    dft.applyInverse(tmp_data, phase_data);
    for (size_t i = 0; i < freq_data.size(); ++i)
        freq_data[i] = std::exp(freq_data[i] + I * phase_data[i]);
    dft.applyInverse(freq_data, tmp_data);
    imp_data.resize(tmp_data.size() / 2);
    for (size_t i = 0; i < imp_data.size(); ++i)
        imp_data[i] = tmp_data[i].real();
}

void InternalMinimalPhaseFilter::apply (const std::vector<double>& signal, std::vector<double>& filtered_signal) const
{

    size_t m = imp_data.size();
    filtered_signal.resize(signal.size());
    for (size_t i = 0; i < filtered_signal.size(); ++i)
    {
        filtered_signal[i] = 0;
        for (size_t k = 0; k < m; ++k)
        {
            double s = 0;
            if (i - k >= 0 && i - k < signal.size())
                s = signal[i - k];
            filtered_signal[i] += s * imp_data[k];
        }
    }
}
} // namespace internal_filter
