#include "filter.h"

#include <bandpassfilter.h>
#include <common.h>
#include <filterscommon.h>
#include <q_filter.h>
#include <wienerfilter.h>

namespace filters_common
{

gund_structs::RealSignature constructFilter (
    const gund_structs::QFilter& config,
    const gund_structs::NewSignatureParameters& params
)
{
    gund_structs::RealSignature result;
    result.params = params;
    result.data = BuildQFilter(params, config);
    return result;
}

gund_structs::RealSignature constructFilter (
    const gund_structs::NewBandPassFilter& config,
    size_t sample_num,
    const gund_structs::NewSignatureParameters& params
)
{
    gund_structs::RealSignature result;
    result.params = params;
    result.data = BuildBandPassFilter(params, sample_num, config);
    return result;
}

gund_structs::RealSignature constructFilter (
    const gund_structs::AntiAliasFilter& config,
    size_t sample_num,
    const gund_structs::NewSignatureParameters& params
)
{
    return constructFilter(static_cast<const gund_structs::NewBandPassFilter&>(config), sample_num, params);
}

gund_structs::RealSignature ApplyFilter (
    const gund_structs::RealSignature& signal,
    const gund_structs::RealSignature& filter_kernel
)
{
    gund_structs::RealSignature filtered_signal;
    filtered_signal.params = signal.params;

    filtered_signal.data = ApplyFirFilter(signal.data, filter_kernel.data);

    return filtered_signal;
}

gund_structs::RealSignature
    constructFilter (
        const gund_structs::WienerFilter& config,
        const gund_structs::RealSignature& sig
    )
{
    gund_structs::RealSignature result;
    result.params = sig.params;
    result.data = BuildWienerFilter(sig.data, config);
    return result;
}

} // namespace filters_common
