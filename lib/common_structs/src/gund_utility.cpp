#define _USE_MATH_DEFINES
#include <gund_utility.h>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#else
#include <limits.h>
#include <unistd.h>
#endif

namespace gund_utility
{

namespace
{
constexpr int MAX_PATH_LOCAL = 260;
} // namespace

static unsigned char* upper_table = 0;

unsigned char locale_toupper (unsigned char c)
{
    if (c < 128 || !upper_table)
    {
        int toupper_c = toupper(c);
        return static_cast<unsigned char>(toupper_c);
    }
    return upper_table[c - 128];
}

void toUpperCase (std::string& val)
{
    std::transform(val.begin(), val.end(), val.begin(), locale_toupper);
}

std::string getExecutableDir ()
{
    std::string exePath;

#ifdef _WIN32
    char buffer[MAX_PATH_LOCAL];
    if (GetModuleFileNameA(NULL, buffer, MAX_PATH_LOCAL))
    {
        exePath = buffer;
    }
    else
    {
        throw std::runtime_error("Failed to get the path of the executable.");
    }
#else
    char buffer[MAX_PATH_LOCAL];
    ssize_t len = readlink("/proc/self/exe", buffer, sizeof(buffer) - 1);
    if (len != -1)
    {
        buffer[len] = '\0';
        exePath = buffer;
    }
    else
    {
        throw std::runtime_error("Failed to get the path of the executable.");
    }
#endif

    auto path = std::filesystem::path(exePath).parent_path();
    return path.generic_string();
}

gund_structs::Coordinate changeCoordinates (const double& r, const double& dip, const double& azimuth)
{
    // углы передаются в градусах, переводим в радианы
    double x = r * sin(dip * M_PI / 180) * cos(azimuth * M_PI / 180);
    double y = r * sin(dip * M_PI / 180) * sin(azimuth * M_PI / 180);
    double z = r * cos(dip * M_PI / 180);

    return gund_structs::Coordinate{x, y, z};
}

size_t findNextPowerOfTwo (size_t n)
{
    if (n == 0 || n == 1)
    {
        return 1;
    }
    // уменьшаем число для обработки случая, когда n - уже степень двойки
    n = n - 1;
    // сбрасываем крайний правый бит, пока не останется ближайшая меньшая степень двойки
    while (n & (n - 1))
    {
        n = n & (n - 1);
    }
    return n << 1;
}

double calculateNormalXCorrelation (const std::vector<double>& vecA, const std::vector<double> vecB)
{
    size_t size = std::min(vecA.size(), vecB.size());
    double meanA = std::accumulate(vecA.begin(), vecA.begin() + size, 0.) / size;
    double meanB = std::accumulate(vecB.begin(), vecB.begin() + size, 0.) / size;
    double xMult = 0, vecASq = 0, vecBSq = 0;
    for (size_t i = 0; i < size; i++)
    {
        xMult += (vecA[i] - meanA) * (vecB[i] - meanB);
        vecASq += (vecA[i] - meanA) * (vecA[i] - meanA);
        vecBSq += (vecB[i] - meanB) * (vecB[i] - meanB);
    }
    double normXCorr = xMult / (std::sqrt(vecASq) * std::sqrt(vecBSq));
    return normXCorr;
}

double calculateSeaWaterSoundSpeed (double CelsiusTemperature)
{
    auto T = CelsiusTemperature;

    const double c0 = 1611.199715;
    const double c1 = -26.196465;
    const double c2 = -2.420430;
    const double c3 = 0.126685;
    double poly = c0 + (c1 + (c2 + c3 * T) * T) * T;

    // Коэффициенты для формулы Юнеско
    const double C00 = 1402.388;
    const double C01 = 5.03711;
    const double C02 = -5.80852e-2;
    const double C03 = 3.3420e-4;
    const double C04 = -1.4780e-6;
    const double C05 = 3.1464e-9;

    const double C10 = 0.153563;
    const double C11 = 6.8982e-4;
    const double C12 = -8.1788e-6;
    const double C13 = 1.3621e-7;
    const double C14 = -6.1185e-10;

    const double C20 = 3.1260e-5;
    const double C21 = -1.7107e-6;
    const double C22 = 2.5974e-8;
    const double C23 = -2.5335e-10;
    const double C24 = 1.0405e-12;

    const double C30 = -9.7729e-9;
    const double C31 = 3.8504e-10;
    const double C32 = -2.3643e-12;

    const double A00 = 1.389;
    const double A01 = -1.262e-2;
    const double A02 = 7.164e-5;
    const double A03 = 2.006e-6;
    const double A04 = -3.21e-8;

    const double A10 = 9.4742e-5;
    const double A11 = -1.2580e-5;
    const double A12 = -6.4885e-8;
    const double A13 = 1.0507e-8;
    const double A14 = -2.0122e-10;

    const double A20 = -3.9064e-7;
    const double A21 = 9.1041e-9;
    const double A22 = -1.6002e-10;
    const double A23 = 7.988e-12;

    const double A30 = 1.100e-10;
    const double A31 = 6.649e-12;
    const double A32 = -3.389e-13;

    const double B00 = -1.922e-2;
    const double B01 = -4.42e-5;

    const double B10 = 7.3637e-5;
    const double B11 = 1.7945e-7;

    const double D00 = 1.727e-3;
    const double D10 = -7.9836e-6;

    double T2 = T * T;
    double T3 = T2 * T;
    double T4 = T3 * T;
    double T5 = T4 * T;

    double P = 115.1;
    double P2 = P * P;
    double P3 = P2 * P;

    double S = 40;
    double S15 = pow(S, 1.5);
    double S2 = S * S;

    double Cw = C00
              + C01
              * T
              + C02
              * T2
              + C03
              * T3
              + C04
              * T4
              + C05
              * T5;
    Cw += (C10 + C11 * T + C12 * T2 + C13 * T3 + C14 * T4) * P;
    Cw += (C20 + C21 * T + C22 * T2 + C23 * T3 + C24 * T4) * P2;
    Cw += (C30 + C31 * T + C32 * T2) * P3;

    double A = A00
             + A01
             * T
             + A02
             * T2
             + A03
             * T3
             + A04
             * T4;
    A += (A10 + A11 * T + A12 * T2 + A13 * T3 + A14 * T4) * P;
    A += (A20 + A21 * T + A22 * T2 + A23 * T3) * P2;
    A += (A30 + A31 * T + A32 * T2) * P3;

    double B = B00 + B01 * T + (B10 + B11 * T) * P;

    double D = D00 + D10 * P;

    double unesco = Cw + A * S + B * S15 + D * S2;

    double weight;

    // Границы области сшивки фита по данным и модели Юнеско
    constexpr double transition_start = 4.6;
    constexpr double transition_end = 18.0;
    if (T <= transition_start)
    {
        weight = 1.0;
    }
    else if (T >= transition_end)
    {
        weight = 0.0;
    }
    else
    {
        double t_norm = (T - transition_start) / (transition_end - transition_start);
        weight = 1.0 - t_norm * t_norm * (3.0 - 2.0 * t_norm);
    }

    return weight * poly + (1.0 - weight) * unesco;
}

double calculateSeaWaterTemperature (double SoundSpeed)
{

    auto data_fitted_poly_derivative = [] (double T, double h = 1e-8)
    {
        // Считаем численно производную
        return (calculateSeaWaterSoundSpeed(T + h) - calculateSeaWaterSoundSpeed(T)) / (h);
    };

    auto inverse_f = [&data_fitted_poly_derivative] (double y_target, double x_start, double tolerance = 1e-10, int max_iterations = 1000)
    {
        double x = x_start;

        for (int i = 0; i < max_iterations; ++i)
        {
            double fx = calculateSeaWaterSoundSpeed(x);
            double fpx = data_fitted_poly_derivative(x);

            double x_new = x - (fx - y_target) / fpx;

            if (std::abs(x_new - x) < tolerance)
            {
                return x_new;
            }

            x = x_new;
        }

        return x;
    };

    // Берем середину диапазона значений, которые были предоставлены МАГЭ
    double T0 = 0.5 * (2.5 + 5);
    return inverse_f(SoundSpeed, T0);
}

} // namespace gund_utility
