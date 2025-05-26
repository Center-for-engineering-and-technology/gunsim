#include <gund_utility.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <limits.h>
#endif

namespace gund_utility
{

namespace
{
constexpr int MAX_PATH_LOCAL = 260;
} // namespace

static unsigned char* upper_table = 0;

unsigned char locale_toupper(unsigned char c) {
    if (c < 128 || !upper_table) {
        int toupper_c = toupper(c);
        return static_cast<unsigned char>(toupper_c);
    }
    return upper_table[c - 128];
}

void toUpperCase(std::string& val) {
    std::transform(val.begin(), val.end(), val.begin(), locale_toupper);
}

std::string getExecutableDir() {
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
    int len = readlink("/proc/self/exe", buffer, sizeof(buffer) - 1);
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

size_t findNextPowerOfTwo(size_t n) {
    if (n == 0 || n == 1) { return 1; }
    // уменьшаем число для обработки случая, когда n - уже степень двойки
    n = n - 1;
    // сбрасываем крайний правый бит, пока не останется ближайшая меньшая степень двойки
    while (n & (n - 1)) {
        n = n & (n - 1);
    }
    return n << 1;
}

} // namespace gund_utility