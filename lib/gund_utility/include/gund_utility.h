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

} // namespace gund_utility