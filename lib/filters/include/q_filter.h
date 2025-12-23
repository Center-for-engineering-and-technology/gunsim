#pragma once

#include <vector>

namespace gund_structs
{
struct NewSignatureParameters;
class QFilter;
} // namespace gund_structs

namespace filters_common
{
/*
 * @brief Строит ИХ Q фильтра
 */
std::vector<double> BuildQFilter(
    gund_structs::NewSignatureParameters sig_params,
    const gund_structs::QFilter& config
);
} // namespace filters_common
