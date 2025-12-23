#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

#include <complex>
#include <vector>

namespace filters_common
{

// Отсечка по L2 норме, ~ 5 сигма
constexpr double kSigmaEnergyFraction = 0.99995;
// Длина ядра фильтра при которой свертку делаем во временном пространстве - если больше то через FFT
constexpr size_t kShortKernelThreshold = 10;
// Минимальная длина для FFT
constexpr size_t kMinFftLen = 2048;

/*
 * @brief Возвращает длину сигнала после RFFT
 */
inline size_t RFFTLenFromFFTLen (size_t fft_len)
{
    return fft_len / 2 + 1;
}

inline size_t FFTLenFromRFFTLen (size_t rfft_len)
{
    return (rfft_len > 0) ? (rfft_len - 1) * 2 : 0;
}

/*
 * @brief Возвращает 2 возведенную в наименьшую степень x так, что pow(2, x) > n
 */
inline size_t NextPowerOfTwo (size_t n)
{
    if (n <= 1)
        return 1;
    size_t p = 1;
    while (p < n)
        p <<= 1;
    return p;
}

/*
 * @brief FFT для вещественного сигнала
 */
void RFFTReal(
    const std::vector<double>& x_time,
    std::vector<std::complex<double>>& X_unique, // len = n_fft/2+1
    size_t n_fft
);

inline std::vector<std::complex<double>> RFFTReal (
    const std::vector<double>& x_time
)
{
    auto n_fft = NextPowerOfTwo(x_time.size());
    std::vector<std::complex<double>> X(RFFTLenFromFFTLen(n_fft));
    RFFTReal(x_time, X, n_fft);
    return X;
}

/*
 * @brief IFFT для вечественного сигнала
 */
void IRFFTReal(
    const std::vector<std::complex<double>>& X_unique,
    std::vector<double>& y_time,
    size_t n_fft
);

inline std::vector<double> IRFFTReal (
    const std::vector<std::complex<double>>& X_unique
)
{
    auto n_fft = FFTLenFromRFFTLen(X_unique.size());
    std::vector<double> x(n_fft);
    IRFFTReal(X_unique, x, n_fft);
    return x;
}

/*
 * @brief Явная свертка за O(n^2)
 */
void StridedConvolution(
    const std::vector<double>& signal,
    int signal_offset,
    int signal_stride,
    const std::vector<double>& kernel,
    int kernel_offset,
    int kernel_stride,
    std::vector<double>& out,
    int out_offset,
    int out_stride
);

/*
 * @brief Сворачивает сигнал с ядром фильтра
 */
std::vector<double> ApplyFirFilter(
    const std::vector<double>& signal,
    const std::vector<double>& kernel,
    const char* mode = "auto"
);

/*
 * @brief Возвращает автокорреляцию сигнала длиной out.size()!
 */
void Autocorrelation(const std::vector<double>& signal, std::vector<double>& out);

/*
 * @brief Косинусная подрезка краев
 */
void ApplyCosineTaper(std::vector<double>& h, size_t taper_len, size_t max_edge);

} // namespace filters_common
