#pragma once

#define _USE_MATH_DEFINES

#include <cmath>
#include <complex>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include <gund_structs.h>

namespace internal_filter
{

constexpr std::complex<double> I = {0, 1};

class AbstractInternalFilter
{ // abstract band pass internal filter structure
protected:

    struct DFT
    { // Discrete Fourier transform tools

        size_t size;
        std::vector<std::complex<double>> exp;

        void resize(size_t size_);

        explicit DFT (size_t size_ = 0)
        {
            resize(size_);
        }

        template<typename T>
        void apply(const std::vector<T>& signal, std::vector<std::complex<double>>& spectrum);

        template<typename T>
        void applyInverse(const std::vector<T>& spectrum, std::vector<std::complex<double>>& signal);
    };

    std::string filtername;
    double dt;         // [s] sample interval
    int iz;            // initial element
    double N;          // [Hz] Nyquist frequency
    double low_pass;   // [Hz]
    double high_pass;  // [Hz]
    double low_slope;  // [dB/oct]
    double high_slope; // [dB/oct]

    struct Window
    { // the amount of most common window types
        enum struct WindowType
        {
            RECTANGULAR,
            HANN,
            HAMMING,
            BLACKMAN
        } window_type;
        std::function<double(double, double)> w;
        double peak_stopband_attenuation = 0;
        double transition_bandwith_index = 0;

        explicit Window(WindowType window_type_);
    } low_window, high_window;

    std::vector<double> imp_data; // filter's impulse response

    void windowedSincResponce(std::vector<double>& impulse) const; // construction the symmetric sinc response by sertain filter's parameters

    virtual void computeImpulseResponce() = 0;

public:

    AbstractInternalFilter (std::string filtername_, double sampleInterval_, gund_structs::BandpassFilter::InternalFilter filterParams, Window::WindowType low_window_type_, Window::WindowType high_window_type_)
        : filtername(filtername_), dt(sampleInterval_), iz(0), N(1 / 2. / dt), low_pass(filterParams.lowCutOffFrec), high_pass(filterParams.highCutOffFrec), low_slope(filterParams.lowCutOffSlope), high_slope(filterParams.highCutOffSlope), low_window(low_window_type_), high_window(high_window_type_)
    {
    }

    virtual void apply(const std::vector<double>& signal, std::vector<double>& filtered_signal) const = 0;

    const std::vector<double>& getImpulse () const
    {
        return imp_data;
    }

    void printInFile(std::string dataPath) const;
};

class InternalZeroPhaseFilter : public AbstractInternalFilter
{
private:

    void computeImpulseResponce() override;

public:

    InternalZeroPhaseFilter (std::string filtername_, double sampleInterval_, gund_structs::BandpassFilter::InternalFilter filterParams, Window::WindowType low_window_type_ = Window::WindowType::HAMMING, Window::WindowType high_window_type_ = Window::WindowType::HAMMING)
        : AbstractInternalFilter(filtername_, sampleInterval_, filterParams, low_window_type_, high_window_type_)
    {
        computeImpulseResponce();
    }

    void apply(const std::vector<double>& signal, std::vector<double>& filtered_signal) const override;
};

class InternalMinimalPhaseFilter : public AbstractInternalFilter
{
private:

    void computeImpulseResponce() override;

public:

    InternalMinimalPhaseFilter (std::string filtername_, double sampleInterval_, gund_structs::BandpassFilter::InternalFilter filterParams, Window::WindowType low_window_type_ = Window::WindowType::HAMMING, Window::WindowType high_window_type_ = Window::WindowType::HAMMING)
        : AbstractInternalFilter(filtername_, sampleInterval_, filterParams, low_window_type_, high_window_type_)
    {
        computeImpulseResponce();
    }

    void apply(const std::vector<double>& signal, std::vector<double>& filtered_signal) const override;
};
} // namespace internal_filter
