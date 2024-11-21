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
    size_t len = readlink("/proc/self/exe", buffer, sizeof(buffer) - 1);
    if (len!= -1) 
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

} // namespace gund_utility