#include <gun_parameters_table.h>

namespace gun_parameters_table
{

void GunTable::parseParams (const nlohmann::json& jsonObject)
{
    // fixed params parsing
    if (jsonObject.find("t_open") != jsonObject.end())
    {
        jsonObject["t_open"].get_to(fixed.t_open);
    }
    else
    {
        throw std::runtime_error("wrong input file content: t_open");
    }
    // depended params parsing
    if (jsonObject.find("A_hvp") != jsonObject.end() && checkDimArray(jsonObject["A_hvp"]))
    {
        depended.A_hvp = jsonObject["A_hvp"].get<std::vector<std::vector<std::vector<double>>>>();
    }
    else
    {
        throw std::runtime_error("wrong input file content: A_hvp");
    }
    if (jsonObject.find("M_hc_hvp") != jsonObject.end() && checkDimArray(jsonObject["M_hc_hvp"]))
    {
        depended.M_hc_hvp = jsonObject["M_hc_hvp"].get<std::vector<std::vector<std::vector<double>>>>();
    }
    else
    {
        throw std::runtime_error("wrong input file content: M_hc_hvp");
    }
    if (jsonObject.find("alpha_hvp") != jsonObject.end() && checkDimArray(jsonObject["alpha_hvp"]))
    {
        depended.alpha_hvp = jsonObject["alpha_hvp"].get<std::vector<std::vector<std::vector<double>>>>();
    }
    else
    {
        throw std::runtime_error("wrong input file content: alpha_hvp");
    }
    if (jsonObject.find("alpha_c_hvp") != jsonObject.end() && checkDimArray(jsonObject["alpha_c_hvp"]))
    {
        depended.alpha_c_hvp = jsonObject["alpha_c_hvp"].get<std::vector<std::vector<std::vector<double>>>>();
    }
    else
    {
        throw std::runtime_error("wrong input file content: alpha_c_hvp");
    }
    if (jsonObject.find("by_hvp") != jsonObject.end() && checkDimArray(jsonObject["by_hvp"]))
    {
        depended.by_hvp = jsonObject["by_hvp"].get<std::vector<std::vector<std::vector<double>>>>();
    }
    else
    {
        throw std::runtime_error("wrong input file content: by_hvp");
    }
    if (jsonObject.find("alpha_mu_hvp") != jsonObject.end() && checkDimArray(jsonObject["alpha_mu_hvp"]))
    {
        depended.alpha_mu_hvp = jsonObject["alpha_mu_hvp"].get<std::vector<std::vector<std::vector<double>>>>();
    }
    else
    {
        throw std::runtime_error("wrong input file content: alpha_mu_hvp");
    }
    if (jsonObject.find("alpha_clust_hvp") != jsonObject.end() && checkDimArray(jsonObject["alpha_clust_hvp"]))
    {
        depended.alpha_clust_hvp = jsonObject["alpha_clust_hvp"].get<std::vector<std::vector<std::vector<double>>>>();
    }
    else
    {
        throw std::runtime_error("wrong input file content: alpha_clust_hvp");
    }
    if (jsonObject.find("m_ratio_hvp") != jsonObject.end() && checkDimArray(jsonObject["m_ratio_hvp"]))
    {
        depended.m_ratio_hvp = jsonObject["m_ratio_hvp"].get<std::vector<std::vector<std::vector<double>>>>();
    }
    else
    {
        throw std::runtime_error("wrong input file content: m_ratio_hvp");
    }
    if (jsonObject.find("t_start_hvp") != jsonObject.end() && checkDimArray(jsonObject["t_start_hvp"]))
    {
        depended.t_start_hvp = jsonObject["t_start_hvp"].get<std::vector<std::vector<std::vector<double>>>>();
    }
    else
    {
        throw std::runtime_error("wrong input file content: t_start_hvp");
    }
}

void GunTable::assembleParams (nlohmann::json& jsonObject) const
{
    // fixed params assembling
    jsonObject["t_open"] = fixed.t_open;
    // depended params assembling
    jsonObject["A_hvp"] = nlohmann::json(depended.A_hvp);
    jsonObject["M_hc_hvp"] = nlohmann::json(depended.M_hc_hvp);
    jsonObject["alpha_hvp"] = nlohmann::json(depended.alpha_hvp);
    jsonObject["alpha_c_hvp"] = nlohmann::json(depended.alpha_c_hvp);
    jsonObject["alpha_mu_hvp"] = nlohmann::json(depended.alpha_mu_hvp);
    jsonObject["alpha_clust_hvp"] = nlohmann::json(depended.alpha_clust_hvp);
    jsonObject["by_hvp"] = nlohmann::json(depended.by_hvp);
    jsonObject["m_ratio_hvp"] = nlohmann::json(depended.m_ratio_hvp);
    jsonObject["t_start_hvp"] = nlohmann::json(depended.t_start_hvp);
}

bool GunTable::checkDimArray (const nlohmann::json& jsonObject) const
{
    if (!jsonObject.is_array() || jsonObject.size() != mesh.H_SIZE)
        return false;
    for (const auto& matrix: jsonObject)
    {
        if (!matrix.is_array() || matrix.size() != mesh.V_SIZE)
            return false;
        for (const auto& row: matrix)
            if (!row.is_array() || row.size() != mesh.P_SIZE)
                return false;
    }
    return true;
}

Params GunTable::interpolate (const gund_structs::Gun& gun) const
{
    Params result;

    result.t_open = fixed.t_open;

    double h = gun.z, v = gun.volume, p = gun.pressure;
    size_t i_h_left, i_h_right, i_v_left, i_v_right, i_p_left, i_p_right;
    double bar_h = (h - mesh.H_START) / static_cast<double>(mesh.H_STEP) - std::floor((h - mesh.H_START) / static_cast<double>(mesh.H_STEP));
    double bar_v = (v - mesh.V_START) / static_cast<double>(mesh.V_STEP) - std::floor((v - mesh.V_START) / static_cast<double>(mesh.V_STEP));
    double bar_p = (p - mesh.P_START) / static_cast<double>(mesh.P_STEP) - std::floor((p - mesh.P_START) / static_cast<double>(mesh.P_STEP));
    if (h < static_cast<double>(mesh.H_END) && h >= static_cast<double>(mesh.H_START))
    {
        i_h_left = static_cast<size_t>(std::floor((h - mesh.H_START) / static_cast<double>(mesh.H_STEP)));
        i_h_right = i_h_left + 1;
        // на всякий случай
        if (i_h_left < 0)
            i_h_left = 0;
        if (i_h_left >= mesh.H_SIZE)
            i_h_left = mesh.H_SIZE - 1;
        if (i_h_right < 0)
            i_h_right = 0;
        if (i_h_right >= mesh.H_SIZE)
            i_h_right = mesh.H_SIZE - 1;
    }
    else
    {
        if (h >= static_cast<double>(mesh.H_END))
        {
            i_h_left = mesh.H_SIZE - 1;
            i_h_right = mesh.H_SIZE - 1;
        }
        else
        {
            i_h_left = 0;
            i_h_right = 0;
        }
    }
    if (v < static_cast<double>(mesh.V_END) && v >= static_cast<double>(mesh.V_START))
    {
        i_v_left = static_cast<size_t>(std::floor((v - mesh.V_START) / static_cast<double>(mesh.V_STEP)));
        i_v_right = i_v_left + 1;
        // на всякий случай
        if (i_v_left < 0)
            i_v_left = 0;
        if (i_v_left >= mesh.V_SIZE)
            i_v_left = mesh.V_SIZE - 1;
        if (i_v_right < 0)
            i_v_right = 0;
        if (i_v_right >= mesh.V_SIZE)
            i_v_right = mesh.V_SIZE - 1;
    }
    else
    {
        if (v >= static_cast<double>(mesh.V_END))
        {
            i_v_left = mesh.V_SIZE - 1;
            i_v_right = mesh.V_SIZE - 1;
        }
        else
        {
            i_v_left = 0;
            i_v_right = 0;
        }
    }
    if (p < static_cast<double>(mesh.P_END) && p >= static_cast<double>(mesh.P_START))
    {
        i_p_left = static_cast<size_t>(std::floor((p - mesh.P_START) / static_cast<double>(mesh.P_STEP)));
        i_p_right = i_p_left + 1;
        // на всякий случай
        if (i_p_left < 0)
            i_p_left = 0;
        if (i_p_left >= mesh.P_SIZE)
            i_p_left = mesh.P_SIZE - 1;
        if (i_p_right < 0)
            i_p_right = 0;
        if (i_p_right >= mesh.P_SIZE)
            i_p_right = mesh.P_SIZE - 1;
    }
    else
    {
        if (p >= static_cast<double>(mesh.P_END))
        {
            i_p_left = mesh.P_SIZE - 1;
            i_p_right = mesh.P_SIZE - 1;
        }
        else
        {
            i_p_left = 0;
            i_p_right = 0;
        }
    }

    result.A = (1 - bar_h) * (1 - bar_v) * (1 - bar_p) * depended.A_hvp[i_h_left][i_v_left][i_p_left] + (1 - bar_h) * (1 - bar_v) * bar_p * depended.A_hvp[i_h_left][i_v_left][i_p_right] + (1 - bar_h) * bar_v * (1 - bar_p) * depended.A_hvp[i_h_left][i_v_right][i_p_left] + (1 - bar_h) * bar_v * bar_p * depended.A_hvp[i_h_left][i_v_right][i_p_right] + bar_h * (1 - bar_v) * (1 - bar_p) * depended.A_hvp[i_h_right][i_v_left][i_p_left] + bar_h * (1 - bar_v) * bar_p * depended.A_hvp[i_h_right][i_v_left][i_p_right] + bar_h * bar_v * (1 - bar_p) * depended.A_hvp[i_h_right][i_v_right][i_p_left] + bar_h * bar_v * bar_p * depended.A_hvp[i_h_right][i_v_right][i_p_right];
    result.M_hc = (1 - bar_h) * (1 - bar_v) * (1 - bar_p) * depended.M_hc_hvp[i_h_left][i_v_left][i_p_left] + (1 - bar_h) * (1 - bar_v) * bar_p * depended.M_hc_hvp[i_h_left][i_v_left][i_p_right] + (1 - bar_h) * bar_v * (1 - bar_p) * depended.M_hc_hvp[i_h_left][i_v_right][i_p_left] + (1 - bar_h) * bar_v * bar_p * depended.M_hc_hvp[i_h_left][i_v_right][i_p_right] + bar_h * (1 - bar_v) * (1 - bar_p) * depended.M_hc_hvp[i_h_right][i_v_left][i_p_left] + bar_h * (1 - bar_v) * bar_p * depended.M_hc_hvp[i_h_right][i_v_left][i_p_right] + bar_h * bar_v * (1 - bar_p) * depended.M_hc_hvp[i_h_right][i_v_right][i_p_left] + bar_h * bar_v * bar_p * depended.M_hc_hvp[i_h_right][i_v_right][i_p_right];
    result.alpha = (1 - bar_h) * (1 - bar_v) * (1 - bar_p) * depended.alpha_hvp[i_h_left][i_v_left][i_p_left] + (1 - bar_h) * (1 - bar_v) * bar_p * depended.alpha_hvp[i_h_left][i_v_left][i_p_right] + (1 - bar_h) * bar_v * (1 - bar_p) * depended.alpha_hvp[i_h_left][i_v_right][i_p_left] + (1 - bar_h) * bar_v * bar_p * depended.alpha_hvp[i_h_left][i_v_right][i_p_right] + bar_h * (1 - bar_v) * (1 - bar_p) * depended.alpha_hvp[i_h_right][i_v_left][i_p_left] + bar_h * (1 - bar_v) * bar_p * depended.alpha_hvp[i_h_right][i_v_left][i_p_right] + bar_h * bar_v * (1 - bar_p) * depended.alpha_hvp[i_h_right][i_v_right][i_p_left] + bar_h * bar_v * bar_p * depended.alpha_hvp[i_h_right][i_v_right][i_p_right];
    result.alpha_c = (1 - bar_h) * (1 - bar_v) * (1 - bar_p) * depended.alpha_c_hvp[i_h_left][i_v_left][i_p_left] + (1 - bar_h) * (1 - bar_v) * bar_p * depended.alpha_c_hvp[i_h_left][i_v_left][i_p_right] + (1 - bar_h) * bar_v * (1 - bar_p) * depended.alpha_c_hvp[i_h_left][i_v_right][i_p_left] + (1 - bar_h) * bar_v * bar_p * depended.alpha_c_hvp[i_h_left][i_v_right][i_p_right] + bar_h * (1 - bar_v) * (1 - bar_p) * depended.alpha_c_hvp[i_h_right][i_v_left][i_p_left] + bar_h * (1 - bar_v) * bar_p * depended.alpha_c_hvp[i_h_right][i_v_left][i_p_right] + bar_h * bar_v * (1 - bar_p) * depended.alpha_c_hvp[i_h_right][i_v_right][i_p_left] + bar_h * bar_v * bar_p * depended.alpha_c_hvp[i_h_right][i_v_right][i_p_right];
    result.alpha_mu = (1 - bar_h) * (1 - bar_v) * (1 - bar_p) * depended.alpha_mu_hvp[i_h_left][i_v_left][i_p_left] + (1 - bar_h) * (1 - bar_v) * bar_p * depended.alpha_mu_hvp[i_h_left][i_v_left][i_p_right] + (1 - bar_h) * bar_v * (1 - bar_p) * depended.alpha_mu_hvp[i_h_left][i_v_right][i_p_left] + (1 - bar_h) * bar_v * bar_p * depended.alpha_mu_hvp[i_h_left][i_v_right][i_p_right] + bar_h * (1 - bar_v) * (1 - bar_p) * depended.alpha_mu_hvp[i_h_right][i_v_left][i_p_left] + bar_h * (1 - bar_v) * bar_p * depended.alpha_mu_hvp[i_h_right][i_v_left][i_p_right] + bar_h * bar_v * (1 - bar_p) * depended.alpha_mu_hvp[i_h_right][i_v_right][i_p_left] + bar_h * bar_v * bar_p * depended.alpha_mu_hvp[i_h_right][i_v_right][i_p_right];
    result.alpha_clust = (1 - bar_h) * (1 - bar_v) * (1 - bar_p) * depended.alpha_clust_hvp[i_h_left][i_v_left][i_p_left] + (1 - bar_h) * (1 - bar_v) * bar_p * depended.alpha_clust_hvp[i_h_left][i_v_left][i_p_right] + (1 - bar_h) * bar_v * (1 - bar_p) * depended.alpha_clust_hvp[i_h_left][i_v_right][i_p_left] + (1 - bar_h) * bar_v * bar_p * depended.alpha_clust_hvp[i_h_left][i_v_right][i_p_right] + bar_h * (1 - bar_v) * (1 - bar_p) * depended.alpha_clust_hvp[i_h_right][i_v_left][i_p_left] + bar_h * (1 - bar_v) * bar_p * depended.alpha_clust_hvp[i_h_right][i_v_left][i_p_right] + bar_h * bar_v * (1 - bar_p) * depended.alpha_clust_hvp[i_h_right][i_v_right][i_p_left] + bar_h * bar_v * bar_p * depended.alpha_clust_hvp[i_h_right][i_v_right][i_p_right];
    result.by = (1 - bar_h) * (1 - bar_v) * (1 - bar_p) * depended.by_hvp[i_h_left][i_v_left][i_p_left] + (1 - bar_h) * (1 - bar_v) * bar_p * depended.by_hvp[i_h_left][i_v_left][i_p_right] + (1 - bar_h) * bar_v * (1 - bar_p) * depended.by_hvp[i_h_left][i_v_right][i_p_left] + (1 - bar_h) * bar_v * bar_p * depended.by_hvp[i_h_left][i_v_right][i_p_right] + bar_h * (1 - bar_v) * (1 - bar_p) * depended.by_hvp[i_h_right][i_v_left][i_p_left] + bar_h * (1 - bar_v) * bar_p * depended.by_hvp[i_h_right][i_v_left][i_p_right] + bar_h * bar_v * (1 - bar_p) * depended.by_hvp[i_h_right][i_v_right][i_p_left] + bar_h * bar_v * bar_p * depended.by_hvp[i_h_right][i_v_right][i_p_right];
    result.m_ratio = (1 - bar_h) * (1 - bar_v) * (1 - bar_p) * depended.m_ratio_hvp[i_h_left][i_v_left][i_p_left] + (1 - bar_h) * (1 - bar_v) * bar_p * depended.m_ratio_hvp[i_h_left][i_v_left][i_p_right] + (1 - bar_h) * bar_v * (1 - bar_p) * depended.m_ratio_hvp[i_h_left][i_v_right][i_p_left] + (1 - bar_h) * bar_v * bar_p * depended.m_ratio_hvp[i_h_left][i_v_right][i_p_right] + bar_h * (1 - bar_v) * (1 - bar_p) * depended.m_ratio_hvp[i_h_right][i_v_left][i_p_left] + bar_h * (1 - bar_v) * bar_p * depended.m_ratio_hvp[i_h_right][i_v_left][i_p_right] + bar_h * bar_v * (1 - bar_p) * depended.m_ratio_hvp[i_h_right][i_v_right][i_p_left] + bar_h * bar_v * bar_p * depended.m_ratio_hvp[i_h_right][i_v_right][i_p_right];
    result.t_start = (1 - bar_h) * (1 - bar_v) * (1 - bar_p) * depended.t_start_hvp[i_h_left][i_v_left][i_p_left] + (1 - bar_h) * (1 - bar_v) * bar_p * depended.t_start_hvp[i_h_left][i_v_left][i_p_right] + (1 - bar_h) * bar_v * (1 - bar_p) * depended.t_start_hvp[i_h_left][i_v_right][i_p_left] + (1 - bar_h) * bar_v * bar_p * depended.t_start_hvp[i_h_left][i_v_right][i_p_right] + bar_h * (1 - bar_v) * (1 - bar_p) * depended.t_start_hvp[i_h_right][i_v_left][i_p_left] + bar_h * (1 - bar_v) * bar_p * depended.t_start_hvp[i_h_right][i_v_left][i_p_right] + bar_h * bar_v * (1 - bar_p) * depended.t_start_hvp[i_h_right][i_v_right][i_p_left] + bar_h * bar_v * bar_p * depended.t_start_hvp[i_h_right][i_v_right][i_p_right];

    return result;
}

void GunMap::setDataPath ()
{
    try
    {
        dataPath = gund_utility::getExecutableDir() + "/data/gun_parameters_table";
    }
    catch (...)
    {
    }
}

void GunMap::parseInput (std::string inputFileName)
{
    if (!std::filesystem::exists(dataPath + "/" + inputFileName))
    {
        throw std::runtime_error("File could not be opened");
    }

    nlohmann::json jsonObject{nlohmann::json::parse(std::ifstream(dataPath + "/" + inputFileName))};
    if (jsonObject.is_array())
    {
        jsonObject = jsonObject[0];
    }

    parseMesh(jsonObject);
    parseData(jsonObject);
}

void GunMap::assembleOutput (std::string outputFileName) const
{
    nlohmann::json jsonObject;

    assembleMesh(jsonObject);
    assembleData(jsonObject);

    std::ofstream outputFile(dataPath + "/" + outputFileName);
    outputFile << std::setw(4) << std::setprecision(17) << jsonObject;
    outputFile.close();
}

void GunMap::parseMesh (const nlohmann::json& jsonObject)
{
    if (jsonObject.find("H_START") != jsonObject.end())
    {
        jsonObject["H_START"].get_to(mesh.H_START);
    }
    else
    {
        throw std::runtime_error("wrong input file content: H_START");
    }
    if (jsonObject.find("H_END") != jsonObject.end())
    {
        jsonObject["H_END"].get_to(mesh.H_END);
    }
    else
    {
        throw std::runtime_error("wrong input file content: H_END");
    }
    if (jsonObject.find("H_STEP") != jsonObject.end())
    {
        jsonObject["H_STEP"].get_to(mesh.H_STEP);
    }
    else
    {
        throw std::runtime_error("wrong input file content: H_STEP");
    }
    if ((mesh.H_END - mesh.H_START) % mesh.H_STEP != 0)
    {
        throw std::runtime_error("wrong input file content: H_STEP");
    }
    mesh.H_SIZE = static_cast<size_t>((mesh.H_END - mesh.H_START) / mesh.H_STEP + 1);

    if (jsonObject.find("V_START") != jsonObject.end())
    {
        jsonObject["V_START"].get_to(mesh.V_START);
    }
    else
    {
        throw std::runtime_error("wrong input file content: V_START");
    }
    if (jsonObject.find("V_END") != jsonObject.end())
    {
        jsonObject["V_END"].get_to(mesh.V_END);
    }
    else
    {
        throw std::runtime_error("wrong input file content: V_END");
    }
    if (jsonObject.find("V_STEP") != jsonObject.end())
    {
        jsonObject["V_STEP"].get_to(mesh.V_STEP);
    }
    else
    {
        throw std::runtime_error("wrong input file content: V_STEP");
    }
    if ((mesh.V_END - mesh.V_START) % mesh.V_STEP != 0)
    {
        throw std::runtime_error("wrong input file content: V_STEP");
    }
    mesh.V_SIZE = static_cast<size_t>((mesh.V_END - mesh.V_START) / mesh.V_STEP + 1);

    if (jsonObject.find("P_START") != jsonObject.end())
    {
        jsonObject["P_START"].get_to(mesh.P_START);
    }
    else
    {
        throw std::runtime_error("wrong input file content: P_START");
    }
    if (jsonObject.find("P_END") != jsonObject.end())
    {
        jsonObject["P_END"].get_to(mesh.P_END);
    }
    else
    {
        throw std::runtime_error("wrong input file content: P_END");
    }
    if (jsonObject.find("P_STEP") != jsonObject.end())
    {
        jsonObject["P_STEP"].get_to(mesh.P_STEP);
    }
    else
    {
        throw std::runtime_error("wrong input file content: P_STEP");
    }
    if ((mesh.P_END - mesh.P_START) % mesh.P_STEP != 0)
    {
        throw std::runtime_error("wrong input file content: P_STEP");
    }
    mesh.P_SIZE = static_cast<size_t>((mesh.P_END - mesh.P_START) / mesh.P_STEP + 1);
}

void GunMap::parseData (const nlohmann::json& jsonObject)
{
    std::vector<std::string> gunNames = {
        "600B",
        "800C",
        "1500C",
        "1500LL",
        "1900C",
        "1900D-DHS",
        "1900LLX",
        "1900LLXT",
        "2800",
        "2800LLX",
        "8500APG",
        "G-GUN",
        "G-GUNII",
        "GI-GUN",
        "Sleeve",
        "Sleevell"
    };
    for (const auto& name: gunNames)
        if (jsonObject.find(name) != jsonObject.end())
        {
            std::string name_cpy = name;
            auto gunType = gund_json_parser::convertGunType(name_cpy);
            data[gunType];
            data[gunType].setType(gunType);
            data[gunType].setMesh(mesh);
            data[gunType].parseParams(jsonObject[name]);
        }
}

void GunMap::assembleMesh (nlohmann::json& jsonObject) const
{
    jsonObject["H_START"] = mesh.H_START;
    jsonObject["H_END"] = mesh.H_END;
    jsonObject["H_STEP"] = mesh.H_STEP;
    jsonObject["V_START"] = mesh.V_START;
    jsonObject["V_END"] = mesh.V_END;
    jsonObject["V_STEP"] = mesh.V_STEP;
    jsonObject["P_START"] = mesh.P_START;
    jsonObject["P_END"] = mesh.P_END;
    jsonObject["P_STEP"] = mesh.P_STEP;
}

void GunMap::assembleData (nlohmann::json& jsonObject) const
{
    for (const auto& [guntype, table]: data)
        table.assembleParams(jsonObject[gund_json_parser::reconvertGunType(guntype)]);
}

Params GunMap::interpolate (const gund_structs::Gun& gun) const
{
    if (data.find(gun.type) != data.end())
        return data.at(gun.type).interpolate(gun);
    return Params(); // default gun if not match
}

}; // namespace gun_parameters_table
