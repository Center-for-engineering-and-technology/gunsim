#include "q_filter.h"

#include <cmath>
#include <stdexcept>

#include <gund_structs.h>
#include <signature_structs.h>

#include <common.h>
#include <filterscommon.h>

namespace filters_common
{
using namespace gund_structs;

std::vector<double> BuildQFilter (
    NewSignatureParameters sig_params,
    const QFilter& config
)
{
    auto q_factor = config.getQValue();
    auto two_way_time_sec = config.getTWTravelTime();
    auto ts_sec = sig_params.sample_interval;

    static constexpr size_t n_fft = 4096;
    const auto df = 1.0 / (ts_sec * static_cast<decltype(ts_sec)>(n_fft));
    const size_t n_unique = RFFTLenFromFFTLen(n_fft);

    std::vector<std::remove_const_t<decltype(df)>> freqs(n_unique);
    for (size_t k = 0; k < n_unique; ++k)
        freqs[k] = static_cast<decltype(df)>(k) * df;

    const auto slope = 0.5 * (two_way_time_sec / q_factor);
    std::vector<std::complex<double>> x(n_unique);
    // A(f) = exp( -2π * (T/Q) * f )
    double max_amp = 0.0;
    for (size_t k = 0; k < n_unique; ++k)
    {
        const double log_amp = -slope * 2.0 * M_PI * freqs[k];
        const double amp = std::exp(log_amp);
        x[k] = std::complex<double>(amp, 0.0);
        if (amp > max_amp)
            max_amp = amp;
    }
    // нормировка к DC
    for (auto& v: x)
        v /= max_amp;

    // IRFFT -> действительная ИХ
    std::vector<double> h(n_unique);
    IRFFTReal(x, h, n_fft);
    // irfft_real(x, h, n_fft);

    // обрезка по накопленной энергии на половине
    const auto half = static_cast<int>(n_fft / 2);
    double e_total = 0.0;
    for (size_t i = 0; i < half; ++i)
        e_total += h[i] * h[i];
    const double e_cut = kSigmaEnergyFraction * e_total;

    double accum = 0.0;
    size_t cutoff = 0;
    for (size_t i = 0; i < half; ++i)
    {
        accum += h[i] * h[i];
        if (accum >= e_cut)
        {
            cutoff = i;
            break;
        }
    }

    auto l = (2 * cutoff) + 1;
    auto max_taps = n_fft;
    if (l > max_taps)
        throw std::runtime_error("max_taps too small");

    // симметричное ядро
    std::vector<double> symmetric;
    symmetric.reserve(l);
    // left_rev
    for (size_t i = cutoff; i > 0; --i)
        symmetric.push_back(h[i]);
    // center
    symmetric.push_back(h[0]);
    // right_fwd
    for (size_t i = 1; i <= cutoff; ++i)
        symmetric.push_back(h[i]);

    // минимум-фаза
    std::vector<double> work(n_fft, 0.0);
    std::copy(symmetric.begin(), symmetric.end(), work.begin());
    auto coeffs = MinphaseWindowViaWiener(work, l);
    // косинусная подрезка краёв (как в Python-порте)
    // ApplyCosineTaper(coeffs, static_cast<int>(l_mp), 40);
    auto taper_len = std::min({cutoff, ( size_t )40, coeffs.size()});
    if (taper_len > 1)
    {
        const double wstep = M_PI / (2.0 * (taper_len - 1));
        for (size_t i = 0; i < taper_len; i++)
        {
            const double w = std::cos(wstep * (i));
            coeffs[coeffs.size() - taper_len + i] *= w; // правый край
            coeffs[taper_len - i - 1] *= w;             // левый край
        }
    }

    // нормировка суммы к 1
    double s
        = 0.0;
    for (size_t i = 0; i < coeffs.size(); ++i)
        s += coeffs[i];
    if (s == 0.0)
        throw std::runtime_error("Q-filter zero sum");
    for (size_t i = 0; i < coeffs.size(); ++i)
        coeffs[i] /= s;

    return coeffs;
}
} // namespace filters_common
