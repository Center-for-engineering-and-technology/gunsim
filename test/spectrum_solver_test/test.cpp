#define BOOST_TEST_MODULE spectrum_solver_test
#pragma warning( push )
#pragma warning( disable : 4265 )
#include <boost/test/unit_test.hpp>
#pragma warning( pop )

#include <filesystem>

#include <gund_format_parser.h>
#include <gund_json_parser.h>
#include <spectrum_solver.h>

namespace spectrum_solver_test {

using namespace spectrum_solver;
using namespace spectrum_solver_structs;
using namespace gund_format_parser;
using namespace gund_json_parser;

std::filesystem::path TestDir() {
    return std::filesystem::path(__FILE__).parent_path() / "test_data";
}

// файл без нужного разрешения
std::string createFilename(std::string& base) {
    const auto file = TestDir() / (base);
    return file.generic_string();
}

BOOST_AUTO_TEST_CASE(solver_check) {

    // зачитывание данных из json и файла
    std::string test_name = "solver_check";
    std::string filename = createFilename(test_name);

    std::string json = filename + ".json";
    SpectrumSolverParams params;
    SpectrumSolverOptions options;
    parseInitJson(json, params, options);
    // адрес данных относительно исходников
    const auto dataPath = std::filesystem::absolute(std::filesystem::path(__FILE__).parent_path() / "../../data/");

    // расчет
    SpectrumSolver solver(params, options);
    solver.setDataPath(dataPath.generic_string());
    auto result = solver.solve();

    // выгрузка параметров
}

} // namespace spectrum_solver_test
