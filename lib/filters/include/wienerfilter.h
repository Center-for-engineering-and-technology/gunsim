#pragma once

#include <vector>

namespace gund_structs
{
struct NewSignatureParameters;
class WienerFilter;
} // namespace gund_structs

namespace filters_common
{
/*
 * @brief Строит ИХ фильтра Винера
 */
std::vector<double> BuildWienerFilter(
    const std::vector<double>& sig_values,
    const gund_structs::WienerFilter& config
);
} // namespace filters_common
