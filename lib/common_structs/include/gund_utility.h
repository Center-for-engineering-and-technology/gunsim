#pragma once

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <numeric>
#include <stdexcept>
#include <string>

#include <gund_structs.h>

namespace gund_utility
{

void toUpperCase(std::string& val);

// получаем адрес сборки
std::string getExecutableDir();
// перевод сферических координат в декартовы
gund_structs::Coordinate changeCoordinates(const double& r, const double& dip, const double& azimuth);
// находим ближайшую степень двойки большую или равную рассматриваемому числу
size_t findNextPowerOfTwo(size_t n);
// вычисление нормализованной кросс-кореляции между векторами
double calculateNormalXCorrelation(const std::vector<double>& vecA, const std::vector<double> vecB);
// Возвращает скорость звука в воде по данным от заказчика + экстраполяция моделью вилсона
double calculateSeaWaterSoundSpeed(double CelsiusTemperature);
double calculateSeaWaterTemperature(double SoundSpeed);

} // namespace gund_utility
