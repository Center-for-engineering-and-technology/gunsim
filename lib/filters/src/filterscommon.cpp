#include "filterscommon.h"

#include <cstring>
#include <simple_fft/fft.h>

#include <common.h>

namespace filters_common
{

std::vector<double> WienerRecursion (size_t impulse_len,
                                     const std::vector<double>& acf,     // длина >= impulse_len
                                     const std::vector<double>& weights, // длина >= impulse_len
                                     WienerRecursionMode mode)
{
    if ((impulse_len > acf.size()) || (impulse_len > weights.size()))
    {
        throw std::runtime_error("impulse_len larget than acf or weights size!");
    }
    double r0 = acf[0];
    if (r0 == 0.0)
    {
        throw std::runtime_error("acf first element is zero!");
    }

    double r1 = (impulse_len > 1) ? acf[1] : 0.0;

    std::vector<double> tmp_poly(impulse_len, 0);
    std::vector<double> fir_out(impulse_len, 0);

    // инициализация
    tmp_poly[0] = 1.0;

    fir_out[0]
        = (mode == WienerRecursionMode::SPIKE) ? (1.0 / r0) : (weights[0] / r0);
    double running = r1 * fir_out[0];
    if (impulse_len == 1)
    {
        return fir_out;
    }

    // порядок 0..impulse_len-2
    for (size_t order = 0; order < impulse_len - 2; ++order)
    {
        const double reflect = -r1 / r0;
        r0 = r0 + r1 * reflect;
        if (r0 == 0.0)
            r0 = std::numeric_limits<double>::epsilon();

        tmp_poly[order + 1] = tmp_poly[0] * reflect;
        auto half = order / 2;
        for (size_t j = 0; j < half; ++j)
        {
            const double a = tmp_poly[j + 1];
            const double b = tmp_poly[order - j];
            tmp_poly[j + 1] = b * reflect + a;
            tmp_poly[order - j] = reflect * a + b;
        }
        if (order % 2 == 1)
        {
            auto mid = (order + 1) / 2;
            tmp_poly[mid] = (reflect + 1.0) * tmp_poly[mid];
        }

        // r1: лаги 0..(order+1)
        r1 = 0.0;
        for (size_t j = 0; j < order + 2; ++j)
        {
            r1 += acf[(order + 2) - j] * tmp_poly[j];
        }

        const double target = (mode == WienerRecursionMode::SPIKE) ? (-running) : (weights[order + 1] - running);
        const double gain = target / r0;

        // fir_out += gain * reversed(tmp_poly[0..order+1])
        for (size_t j = 0; j < order + 1; ++j)
        {
            fir_out[j] += tmp_poly[(order + 1) - j] * gain;
        }
        fir_out[order + 1] = tmp_poly[0] * gain;

        // следующий running — те же лаги
        running = 0.0;
        for (size_t j = 0; j < order + 2; ++j)
        {
            running += acf[(order + 2) - j] * fir_out[j];
        }
    }
    return fir_out;
}

void InvertPowerSeries (const std::vector<double>& h_series, std::vector<double>& g_series)
{
    g_series.resize(h_series.size());
    g_series[0] = 1.0 / h_series[0];
    for (size_t n = 1; n < h_series.size(); ++n)
    {
        double acc = 0.0;
        for (size_t i = 0; i < n; ++i)
            acc += h_series[n - i] * g_series[i];
        g_series[n] = -g_series[0] * acc;
    }
}

size_t TrimByL2Energy (const std::vector<double>& window, size_t window_size, double energy_fraction)
{
    double total = 0.0;
    for (size_t i = 0; i < window_size; ++i)
        total += window[i] * window[i];
    total *= energy_fraction;
    double accum = 0.0;
    for (size_t i = 0; i < window_size; ++i)
    {
        accum += window[i] * window[i];
        if (accum >= total)
            return i;
    }
    return window.size();
}

std::vector<double> MinphaseWindowViaWiener (const std::vector<double>& source_signal, size_t init_len)
{
    if (source_signal.size() == 0 || init_len == 0)
        return {};

    const auto acf_len = init_len + 1;
    std::vector<double> acf(acf_len);
    Autocorrelation(source_signal, acf);
    acf[0] *= 1.001; // лёгкая регуляризация диагонали

    std::vector<double> weights(init_len, 0.0);
    weights[0] = 1.0;

    // std::vector<double> fir(init_len, 0.0);
    // std::vector<double> tmp_poly(init_len, 0.0);
    int status = 0;
    auto fir = WienerRecursion(init_len, acf, weights, WienerRecursionMode::SPIKE);
    if (status != 0)
        throw std::runtime_error("wiener_recursion failed");

    std::vector<double> minphase(init_len, 0.0);
    InvertPowerSeries(fir, minphase);

    const auto cut_len = TrimByL2Energy(minphase, minphase.size(), kSigmaEnergyFraction) + 1;
    minphase.resize(cut_len);

    // масштаб по энергии источника
    double src_e = 0.0;
    double win_e = 0.0;
    for (size_t i = 0; i < source_signal.size(); ++i)
        src_e += source_signal[i] * source_signal[i];
    for (size_t i = 0; i < cut_len; ++i)
        win_e += minphase[i] * minphase[i];
    if (win_e > 0.0)
    {
        const double s = std::sqrt(src_e / win_e);
        for (size_t i = 0; i < cut_len; ++i)
            minphase[i] *= s;
    }
    return minphase;
}

std::vector<double> BuildPefFromSpike (const std::vector<double>& spike_h, size_t gap)
{
    const auto l = spike_h.size();
    std::vector<double> h(gap + l, 0.0);
    h[0] = 1.0;
    std::copy(spike_h.begin(), spike_h.end(), h.begin() + gap);
    for (size_t i = gap; i < h.size(); ++i)
        h[i] = -h[i];
    return h;
}

} // namespace filters_common
