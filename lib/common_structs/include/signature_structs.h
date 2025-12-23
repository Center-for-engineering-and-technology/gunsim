#pragma once

#include <complex>
#include <type_traits>
#include <vector>

#include <gund_structs.h>

namespace gund_structs
{

template<typename T>
constexpr bool is_complex_of_floating_v = false;

template<typename U>
constexpr bool is_complex_of_floating_v<std::complex<U>> = std::is_floating_point_v<U>;

template<typename T>
using enable_if_floating_or_complex_t = std::enable_if_t<
    std::is_floating_point_v<T> || is_complex_of_floating_v<T>>;

struct NewSignatureParameters
{
    NewSignatureParameters() = default;
    NewSignatureParameters (
        double sample_interval,
        ObservationPoint observation_point = {}
    )
        : sample_interval(sample_interval)
        , observation_point(observation_point)
    {
    }

    double sample_interval{0.0005}; // сек
    gund_structs::ObservationPoint observation_point{};
};

inline SignatureParameters
    ConvertNewToOldSignatureParameters (const NewSignatureParameters& new_params, size_t sample_num)
{
    SignatureParameters res;
    res.sampleNum = static_cast<double>(sample_num);
    res.sampleInterval = new_params.sample_interval;
    res.sampleMax = sample_num * res.sampleInterval;

    return res;
}

inline NewSignatureParameters
    ConvertOldToNewSignatureParameters (const SignatureParameters& params, ObservationPoint obs_point = {})
{
    NewSignatureParameters res(
        params.sampleInterval,
        obs_point
    );

    return res;
}

template<typename T, typename = enable_if_floating_or_complex_t<T>>
struct Signature
{
    Signature() = default;
    Signature(const Signature& other) = default;
    Signature& operator=(const Signature& other)
    {
        params = other.params;
        data = other.data;
        return *this;
    };

    Signature(Signature&& other) = default;
    Signature& operator=(Signature&& other)
    {
        params = std::move(other.params);
        data = std::move(other.data);
        return *this;
    };

    Signature(const NewSignatureParameters& params, const std::vector<T>& data)
        : params(params), data(data) {};

    Signature(NewSignatureParameters&& params, std::vector<T>&& data)
        : params(std::move(params)), data(std::move(data)) {};

    Signature(const NewSignatureParameters& params, size_t size)
        : params(params), data(size) {};

    Signature(NewSignatureParameters&& params, size_t size)
        : params(std::move(params)), data(size) {};

    Signature(const NewSignatureParameters& params)
        : params(params) {};

    Signature(NewSignatureParameters&& params)
        : params(std::move(params)) {};

    NewSignatureParameters params{};
    std::vector<T> data{};

    double getTimeSampleInterval ()
    {
        return params.sample_interval;
    }

    double getFreqSampleInterval ()
    {
        if (data.empty() || params.sample_interval == 0.0)
        {
            return std::numeric_limits<double>::infinity();
        }

        1 / (params.sample_interval * data.size());
    }

    size_t size ()
    {
        return data.size();
    }

    double getSampleMax ()
    {
        return params.sample_interval * data.size();
    }
};

using RealSignature = Signature<double>;
using ComplexSignature = Signature<std::complex<double>>;

template<typename T>
class NewToOldSignatureConverter
{
public:

    using signature_type
        = Signature<T>;

    NewToOldSignatureConverter (const signature_type& signature)
        : m_signature(signature), m_params(ConvertNewToOldSignatureParameters(signature.params))
    {
    }

    const SignatureParameters& params () const
    {
        return m_params;
    }

    const std::vector<double>&
        data () const
    {
        return m_signature.data.raw();
    }

private:

    const signature_type& m_signature;
    SignatureParameters m_params;
};

} // namespace gund_structs
