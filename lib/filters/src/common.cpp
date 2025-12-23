#include "common.h"

#include <cstring>
#include <simple_fft/fft.h>

namespace filters_common
{

void RFFTReal (
    const std::vector<double>& x_time,
    std::vector<std::complex<double>>& X_unique, // len = n_fft/2+1
    size_t n_fft
)
{
    using namespace simple_fft;
    using complex_t = std::complex<double>;

    std::vector<complex_t> in(n_fft, complex_t{0.0, 0.0});
    for (size_t i = 0; i < x_time.size() && i < n_fft; ++i)
    {
        in[i] = complex_t{x_time[i], 0.0};
    }

    std::vector<complex_t> out;
    out.resize(n_fft);

    const char* errStr = NULL;

    if (!simple_fft::FFT(in, out, n_fft, errStr))
    {
        throw std::runtime_error("simple_fft::FFT forward failed");
    }

    const size_t n_unique = RFFTLenFromFFTLen(n_fft);
    X_unique.assign(out.begin(), out.begin() + n_unique);
}

void IRFFTReal (
    const std::vector<std::complex<double>>& X_unique,
    std::vector<double>& y_time,
    size_t n_fft
)
{
    using namespace simple_fft;
    using complex_t = std::complex<double>;

    const size_t n_unique = X_unique.size();
    if (n_fft != 2 * (n_unique - 1))
        throw std::invalid_argument("irfft_real: n_fft not consistent with n_unique");

    std::vector<complex_t> spec(n_fft);
    for (size_t k = 0; k < n_unique; ++k)
    {
        spec[k] = X_unique[k];
    }
    for (size_t k = n_unique; k < n_fft; ++k)
    {
        size_t m = n_fft - k;
        spec[k] = std::conj(spec[m]);
    }

    std::vector<complex_t> out(n_fft);

    const char* errStr = NULL;
    if (!simple_fft::IFFT(spec, out, n_fft, errStr))
    {
        throw std::runtime_error("simple_fft::FFT inverse failed");
    }

    y_time.resize(n_fft);
    for (size_t i = 0; i < n_fft; ++i)
    {
        y_time[i] = out[i].real();
    }
}

void StridedConvolution (const std::vector<double>& signal, int signal_offset, int signal_stride, const std::vector<double>& kernel, int kernel_offset, int kernel_stride, std::vector<double>& out, int out_offset, int out_stride)
{
    for (size_t i = 0; i < out.size(); ++i)
    {
        out[out_offset + (i * out_stride)] = 0.0;
    }

    for (size_t k = 0; k < kernel.size(); ++k)
    {
        const auto coeff = kernel[kernel_offset + (k * kernel_stride)];
        const auto s0 = signal_offset + (k * signal_stride);

        for (size_t i = 0; i < out.size() - out_offset; ++i)
        {
            out[out_offset + (i * out_stride)] += coeff * signal[s0 + (i * signal_stride)];
        }
    }
}

std::vector<double> ApplyFirFilter (
    const std::vector<double>& signal,
    const std::vector<double>& kernel,
    const char* mode
)
{
    const auto n = signal.size();
    const auto m = kernel.size();

    const bool force_direct = std::string(mode) == "direct";
    const bool force_fft = std::string(mode) == "fft";
    const bool use_direct = force_direct || (!force_fft && m < kShortKernelThreshold);

    if (use_direct)
    {
        // нули слева, затем сигнал; ядро развернуто (y[n]=sum h[k] x[n-k])
        std::vector<double> tmp(n + m - 1, 0.0);
        std::copy(signal.begin(), signal.end(), tmp.begin() + (m - 1));
        std::vector<double> h_rev(kernel.rbegin(), kernel.rend());

        std::vector<double> out(n);
        StridedConvolution(
            tmp,
            0,
            1,
            kernel,
            static_cast<int>(m) - 1,
            -1,
            out,
            0,
            1
        );
        return out;
    }

    // FFT-ветка
    size_t n_fft = std::max<size_t>(kMinFftLen, NextPowerOfTwo(n + m));
    const size_t n_unique = RFFTLenFromFFTLen(n_fft);

    // заполнение нулями
    std::vector<double> x(n_fft, 0.0), h(n_fft, 0.0);
    std::copy(signal.begin(), signal.end(), x.begin());
    std::copy(kernel.begin(), kernel.end(), h.begin());

    // RFFT обеих
    std::vector<std::complex<double>> X(n_unique), H(n_unique);
    RFFTReal(x, X, n_fft);
    RFFTReal(h, H, n_fft);

    // умножение в спектре
    std::vector<std::complex<double>> Y(n_unique);
    for (size_t k = 0; k < n_unique; ++k)
        Y[k] = X[k] * H[k];

    // IRFFT и обрезка
    std::vector<double> y;
    IRFFTReal(Y, y, n_fft);
    y.resize(n);
    return y;
}

void Autocorrelation (const std::vector<double>& signal, std::vector<double>& out)
{
    std::vector<double> padded(signal.size() + out.size(), 0);
    std::copy(signal.begin(), signal.end(), padded.begin());

    StridedConvolution(padded, 0, 1, signal, 0, 1, out, 0, 1);
}

void ApplyCosineTaper (std::vector<double>& h, size_t taper_len, size_t max_edge)
{
    taper_len = std::min(std::min(taper_len, max_edge), h.size());
    if (taper_len <= 1)
        return;
    const double wstep = M_PI / (2.0 * (taper_len));
    for (int i = 0; i < static_cast<int>(taper_len); i++)
    {
        const double w = std::cos(wstep * (i - 1));
        h[h.size() - taper_len + i] *= w; // правый край
        h[taper_len - i - 1] *= w;        // левый край
    }
}
} // namespace filters_common
