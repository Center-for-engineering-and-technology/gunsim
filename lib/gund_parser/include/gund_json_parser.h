#pragma once

#include <fstream>
#include <string>
#include <nlohmann/json.hpp>

#include <gund_utility.h>
#include <spectrum_solver_structs.h>

namespace gund_json_parser {
    // преобразование string с типом пушки в enum gunType
    gund_structs::GunType convertGunType(std::string& type);
    // преобразование enum gunType в строку
    std::string reconvertGunType(gund_structs::GunType const& gunType);

    // чтение данных для инициализации солвера
    bool parseInitJson(std::string const& inputFilename, spectrum_solver_structs::SpectrumSolverParams& params, spectrum_solver_structs::SpectrumSolverOptions& options);
    // запись данных в json
    void writeJson(nlohmann::json& jsonObject, spectrum_solver_structs::SpectrumSolverParams const& params, spectrum_solver_structs::SpectrumSolverOptions const& options);

} // namespace gund_json_parser

