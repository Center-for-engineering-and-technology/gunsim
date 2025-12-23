#include "notional_source_processor.h"

#include <cmath>
#include <numeric>

#include <Eigen/Dense>

#include <common.h>

namespace gun_model
{
using Vec3 = Eigen::Vector3d;
using Cx = std::complex<double>;

namespace
{

gund_structs::ObservationPoint getGhostPoint (
    const gund_structs::ObservationPoint& signal_source,
    double sea_surface_level_z
)
{
    if (signal_source.observationType == gund_structs::ObservationPoint::ObservationPointType::INFINITE)
        return signal_source;

    gund_structs::ObservationPoint result(signal_source);
    result.z = 2 * sea_surface_level_z - signal_source.z;
    return result;
}

double L2DistanceSquared (
    const gund_structs::ObservationPoint& a,
    const gund_structs::ObservationPoint& b
)
{
    const auto dx = a.x - b.x;
    const auto dy = a.y - b.y;
    const auto dz = a.z - b.z;
    return (dx * dx) + (dy * dy) + (dz * dz);
}

double L2Distance (
    const gund_structs::ObservationPoint& a,
    const gund_structs::ObservationPoint& b
)
{
    return std::sqrt(L2DistanceSquared(a, b));
}

gund_structs::RealSignature signatureTimeShift (
    const gund_structs::RealSignature& signature,
    double timeDelay
)
{
    auto n_fft = signature.data.size();

    if (n_fft < filters_common::kMinFftLen)
    {
        n_fft = filters_common::kMinFftLen;
    }

    size_t n_rfft = filters_common::RFFTLenFromFFTLen(n_fft);

    std::vector<std::complex<double>> m_spectrum(n_rfft, 0);

    auto dt = signature.params.sample_interval; // Sec
    auto sample_rate = 1. / dt;                 // Hz

    auto df_rfft = sample_rate / n_rfft;

    const auto& sig_data = signature.data;

    std::vector<double> tmp_signal(n_fft, 0);
    std::copy(sig_data.begin(), sig_data.end(), tmp_signal.begin());

    auto rfft_sig
        = filters_common::RFFTReal(tmp_signal);

    assert(rfft_sig.size() == m_spectrum.size());

    auto phase_mult = -M_PI * timeDelay;

    for (size_t j = 0; j < m_spectrum.size(); ++j)
    {
        auto phase = phase_mult * df_rfft * j;

        std::complex<double> sig_spectum_val = rfft_sig[j];

        std::complex<double>
            shift(std::cos(phase), std::sin(phase));
        m_spectrum[j] = sig_spectum_val * shift;
    }
    gund_structs::RealSignature result(signature.params, filters_common::IRFFTReal(m_spectrum));
    result.data.resize(sig_data.size());
    return result;
}

void CalculateDelays (
    const gund_structs::ObservationPoint& obs_point,
    const gund_structs::ObservationPoint& source_point,
    const gund_structs::PhysicalParameters& physical_params,
    double reflection_coeff,
    double& ampl,
    double& delay,
    double& ghost_ampl,
    double& ghost_delay
)
{
    bool use_ghost = (0. != reflection_coeff);

    if (obs_point.observationType == gund_structs::ObservationPoint::ObservationPointType::INFINITE)
    {
        Vec3 obs_normal(0., 0., 1.);

        double plane_xy = source_point.x * obs_normal.x() + source_point.y * obs_normal.y();

        ampl = 1.;
        delay = (plane_xy - source_point.z * obs_normal.z()) / physical_params.soundVelocity;

        if (use_ghost)
        {
            ghost_ampl = 1.;
            ghost_delay = (plane_xy + source_point.z * obs_normal.z()) / physical_params.soundVelocity;
        }
        else
        {
            ghost_ampl = 0.;
            ghost_delay = 0.;
        }
    }
    else
    {
        auto distance = L2Distance(
                            obs_point,
                            source_point
                        );

        auto ghost_distance = L2Distance(
                                  obs_point,
                                  // Считаем отражение от поверхности моря, поээтому берем 0.
                                  // Для отажений от кабеля - указали бы их глубину
                                  getGhostPoint(source_point, 0.)
                              );

        ampl = 1. / distance;
        delay = distance / physical_params.soundVelocity;

        if (use_ghost)
        {
            ghost_ampl = 1. / ghost_distance;
            ghost_delay = ghost_distance / physical_params.soundVelocity;
        }
        else
        {
            ghost_ampl = 0.;
            ghost_delay = 0.;
        }
    }
}

} // namespace

gund_structs::RealSignature CalculateNotionalSignaturesSum (
    size_t output_signature_size,
    const std::vector<gund_structs::RealSignature> m_signatures_to_sum,
    const gund_structs::Reflection& reflection,
    const gund_structs::ObservationPoint& observation_point,
    const gund_structs::PhysicalParameters& physical_params
)
{

    if (m_signatures_to_sum.empty())
    {
        throw std::runtime_error("Notional sources signatures to sum is empty");
    }

    gund_structs::NewSignatureParameters result_params;
    result_params.observation_point = observation_point;
    result_params.sample_interval = m_signatures_to_sum[0].params.sample_interval;

    gund_structs::RealSignature result(result_params);

    auto dt = result_params.sample_interval; // Sec
    auto sample_rate = 1. / dt;              // Hz

    const auto m_reflection_coefficient = reflection.refCoef;

    size_t notional_sources_count = m_signatures_to_sum.size();
    size_t n_fft{0};

    std::vector<double> ampl(notional_sources_count);
    std::vector<double> delays(notional_sources_count);

    std::vector<double> ghost_ampl(notional_sources_count);
    std::vector<double> ghost_delays(notional_sources_count);

    const auto& result_observation_point = result_params.observation_point;

    for (size_t i = 0; i < notional_sources_count; ++i)
    {
        const auto& p_sig = m_signatures_to_sum[i];
        const auto& sig_params = p_sig.params;
        const auto& sig_data = p_sig.data;

        CalculateDelays(
            result_observation_point,
            sig_params.observation_point,
            physical_params,
            reflection.refCoef,
            ampl[i],
            delays[i],
            ghost_ampl[i],
            ghost_delays[i]
        );

        const size_t sig_data_size = sig_data.size();
        const size_t delay_in_indexes = static_cast<size_t>(std::abs((delays[i] / dt)));
        const size_t ghost_delay_in_indexes = static_cast<size_t>(std::abs((ghost_delays[i] / dt)));

        n_fft = std::max({delay_in_indexes + sig_data_size, ghost_delay_in_indexes + sig_data_size, n_fft});
    }

    auto min_delay = std::min(
        {*std::min_element(delays.begin(), delays.end()),
         *std::min_element(ghost_delays.begin(), ghost_delays.end())}
    );

    std::transform(delays.begin(), delays.end(), delays.begin(), [&min_delay] (auto a)
                   {
                       return a - min_delay;
                   });
    std::transform(ghost_delays.begin(), ghost_delays.end(), ghost_delays.begin(), [&min_delay] (auto a)
                   {
                       return a - min_delay;
                   });

    if (n_fft < filters_common::kMinFftLen)
    {
        n_fft = filters_common::kMinFftLen;
    }
    else
    {
        n_fft = filters_common::NextPowerOfTwo(n_fft);
    }

    size_t n_rfft = filters_common::RFFTLenFromFFTLen(n_fft);

    std::vector<std::complex<double>> m_spectrum(n_rfft, 0);

    auto df_rfft = sample_rate / n_rfft;

    for (size_t i = 0; i < m_signatures_to_sum.size(); ++i)
    {
        const auto& p_sig = m_signatures_to_sum[i];
        const auto& sig_data = p_sig.data;

        std::vector<double> tmp_signal(n_fft, 0);
        std::copy(sig_data.begin(), sig_data.end(), tmp_signal.begin());

        auto rfft_sig
            = filters_common::RFFTReal(tmp_signal);

        assert(rfft_sig.size() == m_spectrum.size());

        auto phase_mult = -M_PI * delays[i];

        auto ghost_phase_mult = -M_PI * ghost_delays[i];

        for (size_t j = 0; j < m_spectrum.size(); ++j)
        {
            auto phase = phase_mult * df_rfft * j;

            std::complex<double> sig_spectum_val = rfft_sig[j];

            std::complex<double>
                shift(std::cos(phase), std::sin(phase));
            m_spectrum[j] += ampl[i] * sig_spectum_val * shift;

            if (m_reflection_coefficient != 0.0)
            {
                auto ghost_phase = ghost_phase_mult * df_rfft * j;
                std::complex<double>
                    shift_mirror(std::cos(ghost_phase), std::sin(ghost_phase));

                m_spectrum[j] += ghost_ampl[i] * m_reflection_coefficient * sig_spectum_val * shift_mirror;
            }
        }
    }

    result.data = filters_common::IRFFTReal(m_spectrum);

    if (reflection.refCoef != 0. && observation_point.observationType == gund_structs::ObservationPoint::ObservationPointType::INFINITE)
    {
        if (reflection.firstCableDepth != 0)
        {
            auto delay = (reflection.firstCableDepth + reflection.firstCableDepth) / (physical_params.soundVelocity);
            auto firstCableGhost = signatureTimeShift(result, delay);

            std::transform(
                result.data.begin(),
                result.data.end(),
                firstCableGhost.data.begin(),
                result.data.begin(),
                [&] (auto a, auto b)
                {
                    return a + reflection.refCoef * b;
                }
            );
            if (reflection.secondCableDepth != 0)
            {
                delay = (reflection.secondCableDepth + reflection.secondCableDepth) / (physical_params.soundVelocity);
                auto secondCableGhost = signatureTimeShift(result, delay);

                std::transform(
                    result.data.begin(),
                    result.data.end(),
                    secondCableGhost.data.begin(),
                    result.data.begin(),
                    [&] (auto a, auto b)
                    {
                        return a + reflection.refCoef * b;
                    }
                );
            }
        }
    }

    result.data.resize(output_signature_size);
    return result;
}

std::vector<gund_structs::RealSignature> EstimateNotional (
    const std::vector<gund_structs::ObservationPoint>& guns_xyz,
    const std::vector<gund_structs::RealSignature>& nfh_time_all,
    const gund_structs::PhysicalParameters& physical_params,
    const gund_structs::Reflection& /* reflection */,
    const NotionalSignaturesEstimatorConfig& config
)
{
    std::vector<gund_structs::RealSignature> result;
    if (guns_xyz.empty() || nfh_time_all.empty())
    {
        return result;
    }

    size_t Ng = guns_xyz.size();
    std::vector<double> gun_weights_(Ng, 1.0);

    const size_t H = nfh_time_all.size();

    std::vector<std::vector<int>> guns_per_hydro(H);

    for (size_t g = 0; g < Ng; ++g)
    {
        double best_d2 = std::numeric_limits<double>::max();
        size_t best_h = 0;

        for (size_t h = 0; h < H; ++h)
        {

            const auto& hydro_position = nfh_time_all[h].params.observation_point;
            double d2 = (( Vec3 )guns_xyz[g] - ( Vec3 )hydro_position).squaredNorm();
            if (d2 < best_d2)
            {
                best_d2 = d2;
                best_h = h;
            }
        }
        guns_per_hydro[best_h].push_back(static_cast<int>(g));
    }

    double fs_ = std::numeric_limits<double>::quiet_NaN();

    std::vector<size_t> cluster_hydro_indices_;
    std::vector<gund_structs::ObservationPoint> cluster_hydro_xyz_;
    std::vector<std::vector<int>> clusters_guns_;

    std::vector<gund_structs::RealSignature> cluster_notional_;

    for (size_t h = 0; h < H; ++h)
    {
        if (!guns_per_hydro[h].empty())
        {
            cluster_hydro_indices_.push_back(h); // индекс гидрофона
            cluster_hydro_xyz_.push_back(nfh_time_all[h].params.observation_point);
            clusters_guns_.push_back(guns_per_hydro[h]); // пушки этого кластера
        }

        if (std::isnan(fs_))
        {
            fs_ = 1. / nfh_time_all[h].params.sample_interval;
        }
    }

    const size_t H_all = nfh_time_all.size();
    const size_t Nt = std::transform_reduce(
        nfh_time_all.begin(),
        nfh_time_all.end(),
        0,
        [] (auto a, auto b)
        {
            return a > b ? a : b;
        },
        [] (auto x)
        {
            return x.size();
        }
    );

    const size_t K = clusters_guns_.size();

    if (H_all != nfh_time_all.size())
        throw std::runtime_error("nfh_time_all rows must match hydro_xyz_all size");

    Eigen::MatrixXd nfh_used(H_all, Nt);
    for (size_t h = 0; h < H_all; ++h)
    {

        for (size_t nt = 0; nt < Nt; ++nt)
        {
            nfh_used.row(h)[nt] = nfh_time_all[h].data[nt];
        }

        if (config.remove_dc)
        {
            double m = nfh_used.row(h).mean();
            nfh_used.row(h).array() -= m;
        }
    }

    size_t target_len = config.pad_factor * Nt;
    size_t n_fft = filters_common::NextPowerOfTwo(target_len);
    size_t n_freq = n_fft / 2 + 1;

    Eigen::MatrixXcd NFH_f(H_all, n_freq);
    {
        std::vector<double> tmp_in(n_fft);
        std::vector<std::complex<double>> tmp_out(n_freq);

        for (size_t h = 0; h < H_all; ++h)
        {
            std::fill(tmp_in.begin(), tmp_in.end(), 0.0);
            for (size_t t = 0; t < Nt; ++t)
                tmp_in[t] = nfh_used(h, t);

            tmp_out = filters_common::RFFTReal(tmp_in);

            for (size_t f = 0; f < n_freq; ++f)
                NFH_f(h, f) = tmp_out[f];
        }
    }

    Ng = guns_xyz.size();
    Eigen::MatrixXd dist_hg(H_all, Ng);
    for (size_t h = 0; h < H_all; ++h)
    {
        for (size_t g = 0; g < Ng; ++g)
        {
            Vec3 d = ( Vec3 )nfh_time_all[h].params.observation_point - ( Vec3 )guns_xyz[g];
            double r = d.norm();
            if (r < 1e-3)
                r = 1e-3;
            dist_hg(h, g) = r;
        }
    }

    auto fmax = config.fmax;
    if (fmax <= 0.0)
        fmax = 0.5 * fs_;

    Eigen::MatrixXcd S_cl_f(K, n_freq);
    S_cl_f.setZero();

    for (size_t f = 0; f < n_freq; ++f)
    {
        double freq = static_cast<double>(f) * fs_ / static_cast<double>(n_fft);
        double omega = 2.0 * M_PI * freq;
        double k_wave = omega / physical_params.soundVelocity;

        Eigen::MatrixXcd Gf(H_all, K);
        Gf.setZero();

        for (size_t k_cl = 0; k_cl < K; ++k_cl)
        {
            const auto& guns_k = clusters_guns_[k_cl];
            Eigen::VectorXcd g_sum(H_all);
            g_sum.setZero();

            for (int g_idx: guns_k)
            {
                double w = gun_weights_[g_idx]; // вес этой пушки
                Eigen::VectorXcd g_j(H_all);
                for (size_t h = 0; h < H_all; ++h)
                {
                    double rr = dist_hg(h, g_idx);
                    Cx val = std::exp(Cx(0, -k_wave * rr)) / (4.0 * M_PI * rr);
                    g_j(h) = w * val; // взвешиваем вклад
                }
                g_sum += g_j;
            }

            Gf.col(k_cl) = g_sum;
        }

        Eigen::VectorXcd p = NFH_f.col(f); // (H_all,)

        // нормальные уравнения
        Eigen::MatrixXcd A = Gf.adjoint() * Gf; // (K x K)
        A += config.reg * Eigen::MatrixXcd::Identity(K, K);
        Eigen::VectorXcd b = Gf.adjoint() * p; // (K,)

        Eigen::VectorXcd s = A.ldlt().solve(b);
        S_cl_f.col(f) = s;
    }

    Eigen::MatrixXd S_cl_time(K, Nt);
    {
        std::vector<std::complex<double>> tmp_in(n_freq);
        std::vector<double> tmp_out(n_fft);

        for (size_t k = 0; k < K; ++k)
        {
            for (size_t f = 0; f < n_freq; ++f)
                tmp_in[f] = S_cl_f(k, f);

            tmp_out = filters_common::IRFFTReal(tmp_in);

            for (size_t t = 0; t < Nt; ++t)
                S_cl_time(k, t) = tmp_out[t];
        }
    }

    for (size_t h = 0; h < H_all; ++h)
    {
        cluster_notional_.emplace_back(nfh_time_all[0].params, Nt);

        for (size_t nt = 0; nt < Nt; ++nt)
        {
            cluster_notional_.back().data[nt] = nfh_used.row(h)[nt];
        }
    }
    return cluster_notional_;
}

gund_structs::ObservationPoint Vec3ToObsPoint (const Vec3& a)
{
    gund_structs::ObservationPoint result;

    result.x = a[0];
    result.y = a[1];
    result.z = a[2];
    return result;
}

} // namespace gun_model
