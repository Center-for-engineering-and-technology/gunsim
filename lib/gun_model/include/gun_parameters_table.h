#pragma once

#include <filesystem>
#include <fstream>
#include <gund_json_parser.h>
#include <gund_structs.h>
#include <gund_utility.h>
#include <nlohmann/json.hpp>

namespace gun_parameters_table
{

struct Params
{
    double t_open = 0.01;
    double t_start = 0.;
    double M_hc = 1.2;
    double A = 0.005;
    double alpha = 1.6;
    double alpha_c = 0.3;
    double alpha_clust = 0;
    double alpha_mu = 0.4;
    double by = 5;
    double m_ratio = 1;

    Params() = default;
};

struct GunTable
{
    gund_structs::GunType type;
    struct Mesh
    {
        int H_START = 6;
        int H_END = 10;
        int H_STEP = 2;
        size_t H_SIZE = 3;
        int V_START = 100;
        int V_END = 200;
        int V_STEP = 50;
        size_t V_SIZE = 3;
        int P_START = 2000;
        int P_END = 3000;
        int P_STEP = 500;
        size_t P_SIZE = 3;
    } mesh;
    struct FixedParams
    {
        double t_open;
    } fixed;
    struct DependedParams
    {
        std::vector<std::vector<std::vector<double>>> A_hvp;
        std::vector<std::vector<std::vector<double>>> M_hc_hvp;
        std::vector<std::vector<std::vector<double>>> alpha_hvp;
        std::vector<std::vector<std::vector<double>>> alpha_c_hvp;
        std::vector<std::vector<std::vector<double>>> alpha_mu_hvp;
        std::vector<std::vector<std::vector<double>>> alpha_clust_hvp;
        std::vector<std::vector<std::vector<double>>> by_hvp;
        std::vector<std::vector<std::vector<double>>> m_ratio_hvp;
        std::vector<std::vector<std::vector<double>>> t_start_hvp;
    } depended;

    GunTable() = default;
    void setType (const gund_structs::GunType& type_)
    {
        type = type_;
    }
    void setMesh (const Mesh& mesh_)
    {
        mesh = mesh_;
    }
    void parseParams(const nlohmann::json& jsonObject);
    void assembleParams(nlohmann::json& jsonObject) const;
    Params interpolate(const gund_structs::Gun& gun) const;

private:

    bool checkDimArray(const nlohmann::json& jsonObject) const;
};

struct GunMap
{

    GunMap() = default;
    void setDataPath (const std::string& path)
    {
        dataPath = path;
    }
    void setDataPath();
    void parseInput(std::string inputFileName);
    void assembleOutput(std::string outputFileName) const;
    Params interpolate(const gund_structs::Gun& gun) const;

private:

    std::string dataPath;
    GunTable::Mesh mesh;
    std::map<gund_structs::GunType, GunTable> data;
    void parseMesh(const nlohmann::json& jsonObject);
    void parseData(const nlohmann::json& jsonObject);
    void assembleMesh(nlohmann::json& jsonObject) const;
    void assembleData(nlohmann::json& jsonObject) const;
};

}; // namespace gun_parameters_table