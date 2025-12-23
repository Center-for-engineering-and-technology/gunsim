#pragma once

#include <vector>

namespace filters_common
{

// Режимы для Фильтра винера
enum class WienerRecursionMode
{
    SPIKE,
    GAPPED,
};

//==========================
// Винеровская рекурсия
//==========================
std::vector<double> WienerRecursion(size_t impulse_len,
                                    const std::vector<double>& acf,     // длина >= impulse_len
                                    const std::vector<double>& weights, // длина >= impulse_len
                                    WienerRecursionMode mode);

/*
 * @brief Инверсия степенного ряда
 * g(z) = 1 / H(z), по коэффициентам
 */
void InvertPowerSeries(const std::vector<double>& h_series, std::vector<double>& g_series);

/*
 * Обрезка по L2-норме
 */
size_t TrimByL2Energy(const std::vector<double>& window, size_t window_size, double energy_fraction);

/*
 * @brief Построение минимум фазового сигнала из амплитудного спектра
 */
std::vector<double> MinphaseWindowViaWiener(const std::vector<double>& source_signal, size_t init_len);

//==========================
// Сборка PEF из spike-ядра
//==========================
/*
 * @brief Сборка Prediction Error Filter из ядра с учетом gap
 * @return Вектор вида [1.0, 0.0, ... , 0.0, -spike_h[0], ..., -spike_h[spike_h.size()-1]]. Нулей после 1.0 ровно gap штук
 */
std::vector<double> BuildPefFromSpike(const std::vector<double>& spike_h, size_t gap);

} // namespace filters_common
