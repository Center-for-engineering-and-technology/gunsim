#pragma once

#include <algorithm>
#include <cctype>
#include <string>
#include <stdexcept>
#include <filesystem>

namespace gund_utility {

    void toUpperCase(std::string& val);

    // получаем адрес сборки
    std::string getExecutableDir();
    // находим ближайшую степень двойки большую или равную рассматриваемому числу
    size_t findNextPowerOfTwo(size_t n);

} // namespace gund_utility
