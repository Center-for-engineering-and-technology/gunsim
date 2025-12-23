#include "bandpassfilter.h"

#include <simple_fft/fft.h>

#include <gund_structs.h>
#include <signature_structs.h>

#include <common.h>
#include <filterscommon.h>

namespace filters_common
{

struct FilterResult
{
    std::vector<double> h; // временная характеристика (FIR ядро)
    std::vector<double> H; // частотная маска (амплитуды 0..Nyquist)
    size_t L;              // длина ядра
};

FilterResult build_bandpass_filter (double f_low, double slope_low, double f_high, double slope_high, double ts, size_t n_fft, bool want_time)
{
    if (n_fft == 0 || ts <= 0.0)
        throw std::runtime_error("oc_filtgn: invalid ts or n_fft");
    if (n_fft & (n_fft - 1))
        throw std::runtime_error("oc_filtgn: n_fft must be power of 2");

    constexpr double K = 0.1723308333814104;
    constexpr double LN2_HALF = 0.3465735902799726;
    constexpr double CL_MIN = -10.0;
    constexpr double CL_MAX = 46.0;

    const double df = 1.0 / (ts * static_cast<double>(n_fft));
    const size_t half = n_fft / 2;

    std::vector<double> H(half + 1);
    H[0] = (f_low != 0.0) ? 0.00004539992976248485 : 1.0;

    const double abs_slope_low = std::abs(slope_low);
    const double abs_slope_high = std::abs(slope_high);

    const double log_low = std::log(f_low);
    const double log_high = std::log(f_high);

    for (size_t k = 1; k <= half; ++k)
    {
        const double f = k * df;
        double terms = 0.0;
        const double log_f = std::log(f);
        if (f_low != 0.0)
            terms += std::exp(std::min(std::max((log_low - log_f) * (K * abs_slope_low), CL_MIN), CL_MAX));
        if (f_high != 0.0)
            terms += std::exp(std::min(std::max((log_f - log_high) * (K * abs_slope_high), CL_MIN), CL_MAX));

        double x = -LN2_HALF * terms;
        x = std::min(std::max(x, CL_MIN), CL_MAX);
        H[k] = std::exp(x);
    }

    const double vmax = *std::max_element(H.begin(), H.end());
    if (vmax <= 0.0)
        throw std::runtime_error("oc_filtgn: invalid normalization (vmax <= 0)");
    for (auto& v: H)
        v /= vmax;

    if (!want_time)
        return {{}, H, n_fft};

    // --- Inverse FFT ---
    std::vector<std::complex<double>> H_complex(H.size());
    for (size_t i = 0; i < H.size(); ++i)
        H_complex[i] = std::complex<double>(H[i], 0.0);

    std::vector<double> h_time(n_fft);
    filters_common::IRFFTReal(H_complex, h_time, n_fft);

    // --- Обрезка по энергии ---
    auto L = filters_common::TrimByL2Energy(h_time, h_time.size() / 2, filters_common::kSigmaEnergyFraction);
    if (L < 0)
        L = 0;
    if (L > half - 1)
        L = half - 1;

    // --- Нормировка по максимуму спектра ---
    std::vector<std::complex<double>> spectrum;
    filters_common::RFFTReal(h_time, spectrum, n_fft);

    double Amax = 0.0;
    for (auto& c: spectrum)
        Amax = std::max(Amax, std::abs(c));
    if (Amax > 0.0)
        for (auto& v: h_time)
            v /= Amax;

    // --- Перестановка и зеркалирование ---
    const double h0 = h_time[0];
    if (L > 1)
    {
        for (size_t i = 0; i < L; ++i)
            h_time[i] = h_time[n_fft - L + i];
        h_time[L] = h0;
        for (size_t i = 0; i < L; ++i)
            h_time[L + 1 + i] = h_time[L - 1 - i];
    }

    h_time.resize(2 * L + 1);

    // --- Косинусная подбривка ---
    filters_common::ApplyCosineTaper(h_time, L, 10);

    return {h_time, H, L};
}

std::vector<double> BuildBandPassFilter (
    gund_structs::NewSignatureParameters sig_params,
    size_t samples_count,
    const gund_structs::NewBandPassFilter& config
)
{
    int n = ( int )samples_count;

    const double ts = sig_params.sample_interval;
    const size_t base_nfft = std::max<size_t>(0x800, static_cast<size_t>(std::pow(2, std::ceil(std::log2(n) + 1))));

    // 1. Базовый фильтр
    auto filt = build_bandpass_filter(config.getLowFreqCut(), config.getLowSlope(), config.getHighFreqCut(), config.getHighSlope(), ts, base_nfft, true);

    // 2. Опционально — минимально-фазовая аппроксимация
    std::vector<double> kernel = filt.h;

    if (config.minphase)
    {
        auto hmin = MinphaseWindowViaWiener(kernel, kernel.size());
        kernel = std::move(hmin);
    }

    // 3. Усечение / дополнение до нужной длины
    std::vector<double> result(n, 0.0);
    const size_t copy_len = std::min(( size_t )n, kernel.size());
    std::copy_n(kernel.begin(), copy_len, result.begin());
    return result;
}

} // namespace filters_common
