#include "wienerfilter.h"

#include <gund_structs.h>
#include <signature_structs.h>

#include <common.h>
#include <filterscommon.h>

namespace filters_common
{
std::vector<double> BuildWienerFilter (
    const std::vector<double>& sig_values,
    const gund_structs::WienerFilter& config
)
{
    auto wiener_len = config.getLen();
    auto wiener_gap = config.getGap();
    auto wiener_wl = config.getWhiteLightPercentage();
    auto spike_mode = config.spike_mode;
    constexpr double reg_c = 0.01;

    if (wiener_len <= 0)
        throw std::invalid_argument("wiener_len must be > 0");
    if (wiener_gap < 0)
        throw std::invalid_argument("wiener_gap must be >= 0");
    if (wiener_gap > wiener_len)
        throw std::invalid_argument("wiener_gap must be <= wiener_len");

    const auto L = wiener_len;

    // 1) автокорреляция длиной L
    std::vector<double> R(L, 0.0);
    Autocorrelation(sig_values, R);

    // 2) регуляризация диагонали
    R[0] *= (1.0 + reg_c * wiener_wl);

    // 3) target-вектор weights
    const auto mode = spike_mode ? WienerRecursionMode::SPIKE : WienerRecursionMode::GAPPED;
    std::vector<double> weights(L, 0.0);
    if (spike_mode)
    {
        weights[0] = sig_values[0];
    }
    else
    {
        const auto j_lo = wiener_gap;
        const auto j_hi = L - wiener_gap;
        if (j_hi > j_lo)
        {
            std::copy(R.begin() + j_lo, R.begin() + j_hi, weights.begin() + j_lo);
        }
    }

    // 4) recursion -> spike ядро
    auto spike = WienerRecursion(L, R, weights, mode);

    // 5) итоговый фильтр
    std::vector<double> h;
    if (spike_mode)
        h = spike;
    else
        h = BuildPefFromSpike(spike, wiener_gap);

    return h;
}
} // namespace filters_common
