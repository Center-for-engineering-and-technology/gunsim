#pragma once

#include <signature_structs.h>

namespace filters_common
{

gund_structs::RealSignature
    constructFilter(
        const gund_structs::QFilter& config,
        const gund_structs::NewSignatureParameters& params
    );

gund_structs::RealSignature
    constructFilter(
        const gund_structs::WienerFilter& config,
        const gund_structs::RealSignature& sig
    );

gund_structs::RealSignature
    constructFilter(
        const gund_structs::NewBandPassFilter& config,
        size_t sample_num,
        const gund_structs::NewSignatureParameters& params
    );

gund_structs::RealSignature
    constructFilter(
        const gund_structs::AntiAliasFilter& config,
        size_t sample_num,
        const gund_structs::NewSignatureParameters& params
    );

gund_structs::RealSignature ApplyFilter(
    const gund_structs::RealSignature& signal,
    const gund_structs::RealSignature& filter_kernel
);
} // namespace filters_common
