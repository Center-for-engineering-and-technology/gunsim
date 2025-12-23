#pragma once

#include <vector>

namespace gund_structs
{
struct NewSignatureParameters;
class NewBandPassFilter;
} // namespace gund_structs

namespace filters_common
{
/*
 * @brief Строит ИХ полосового фильтра
 */
std::vector<double> BuildBandPassFilter(
    gund_structs::NewSignatureParameters sig_params,
    size_t samples_count,
    const gund_structs::NewBandPassFilter& config
);
} // namespace filters_common
