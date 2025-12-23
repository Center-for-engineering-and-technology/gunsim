#pragma once

#include <vector>

#include <signature_structs.h>

namespace gun_model
{

gund_structs::RealSignature CalculateNotionalSignaturesSum(
    size_t output_signature_size,
    const std::vector<gund_structs::RealSignature> m_signatures_to_sum,
    const gund_structs::Reflection& reflection,
    const gund_structs::ObservationPoint& observation_point,
    const gund_structs::PhysicalParameters& physical_params
);

struct NotionalSignaturesEstimatorConfig
{
    double reg{1e-3};
    double fmin{0.0};
    double fmax{-1.0};

    bool remove_dc{true};

    std::size_t pad_factor{2};
};

std::vector<gund_structs::RealSignature> EstimateNotional(
    const std::vector<gund_structs::ObservationPoint>& guns_xyz,
    const std::vector<gund_structs::RealSignature>& nfh_time_all,
    const gund_structs::PhysicalParameters& physical_params,
    const gund_structs::Reflection& reflection,
    const NotionalSignaturesEstimatorConfig& config
);

} // namespace gun_model
