#pragma once

#include <algorithm>
#include <ctime>
#include <fstream>
#include <functional>
#include <gun_model.h>
#include <nlohmann/json.hpp>
#include <random>
#include <string>

namespace sig_optimizer
{

double randomDouble(double min, double max);
double norm_L2(const std::vector<double>& sig_1, const std::vector<double>& sig_2);
double norm_TECHNICAL(const std::vector<double>& sig_1, const std::vector<double>& sig_2);
double norm_TECHNICAL_L2(const std::vector<double>& sig_1, const std::vector<double>& sig_2);

struct Mesh
{
    int H_START = 2;
    int H_END = 20;
    int H_STEP = 2;
    int V_START = 50;
    int V_END = 500;
    int V_STEP = 50;
    int P_START = 500;
    int P_END = 5000;
    int P_STEP = 500;
    size_t SIGN_SIZE = 1000;
};

class Data
{
public:

    Mesh mesh;
    std::string dataPath;
    std::vector<std::vector<std::vector<std::vector<double>>>> data_h_v_p;
    std::vector<std::vector<std::vector<std::vector<double>>>> data_h_v_p_clust;

    void parseData(bool ref);
    void parseDataClust(bool ref);
    void setMesh (Mesh mesh_)
    {
        mesh = mesh_;
    }
    void setDataPath (std::string dataPath_)
    {
        dataPath = dataPath_;
    }
};

struct Options
{
    size_t POP_SIZE = 200;
    size_t POP_COUNT = 100;
    double MUTATION_RATE = 0.1;
    double MUTATION_SCALE = 0.5;
    double POPULATION_REJECTION = 0.05;
    size_t ELITE_GROUP_SIZE = 1;
};

struct Parameter
{
    double min_value;
    double max_value;
    double step;
    bool fixed;
    double value;

    Parameter (double min, double max, bool fix = false, double val = NAN)
        : min_value(min), max_value(max), step((max - min) / 100.), fixed(fix), value(val)
    {
    }
    void choose(); // value initialisation
    void mutation(const Options& options);
    Parameter operator+(const Parameter& other) const; // crossingover
};

struct Genome
{
    std::vector<Parameter> params;
    double fitness;

    void choose(); // initialisation
    void mutation(const Options& options);
    Genome operator+(const Genome& other) const; // crossingover
    template<typename Regime>
    void computeFitness(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm);
    template<typename Regime>
    void computeFitnessClust(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm);
    template<typename Regime>
    void echoSignature(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, double h, double v, double p) const;
};

struct Population
{
    std::vector<Genome> individs;

    template<typename Regime>
    void initialisation(size_t pop_size);
    void choose();
    void mutation(const Options& options);
    template<typename Regime>
    void computeFitness(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm);
    template<typename Regime>
    void computeFitnessClust(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm);
    void fitSort();
    const Genome& selection(const Options& options) const;
    template<typename Regime>
    Population geneticStep(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm, const Options& options);
    template<typename Regime>
    Population geneticStepClust(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm, const Options& options);
    void echoFit() const;
    void echoBest() const;
};

template<typename Regime>
class Algorithm
{
public:

    Options options;
    std::string dataPath;
    std::string gunName;
    Data data;
    gun_model::PhysParams PP;
    gun_model::MethodParams MP;
    std::function<double(const std::vector<double>&, const std::vector<double>&)> norm;
    Population population;

    Algorithm (Options options_, Mesh mesh_, std::string dataPath_, std::string gunName_, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm_, bool ref, bool clust)
        : options(options_), dataPath(dataPath_), gunName(gunName_), norm(norm_)
    {
        data.setDataPath(dataPath + gunName + "/");
        data.setMesh(mesh_);
        data.parseData(ref);
        if (clust)
            data.parseDataClust(ref);
        PP.ref = ref ? -1 : 0;
        population.initialisation<Regime>(options.POP_SIZE);
    }
    std::string genetic();
    std::string geneticClust();
    std::string gradient(Genome& start) const;
};

template<typename Derived>
class OptimisationRegime
{
public:

    static gun_model::Gun::GunType gun_type (const Genome& genome)
    {
        return Derived::gun_type_impl(genome);
    }
    static Genome genome ()
    {
        return Derived::genome_impl();
    }
};

class Default : public OptimisationRegime<Default>
{
    struct Limits
    {
        double A_MIN = 0;
        double A_MAX = 0.1;
        double T_OPEN_MIN = 0.005;
        double T_OPEN_MAX = 0.015;
        double T_START_MIN = 0;
        double T_START_MAX = 0.001;
        double M_HC_MIN = 0;
        double M_HC_MAX = 20;
        double ALPHA_MIN = 0;
        double ALPHA_MAX = 0.3;
        double ALPHA_C_MIN = 0;
        double ALPHA_C_MAX = 0.3;
        double ALPHA_CLUST_MIN = 0;
        double ALPHA_CLUST_MAX = 3;
        double ALPHA_MU_MIN = 0;
        double ALPHA_MU_MAX = 0.3;
        double BY_MIN = 0;
        double BY_MAX = 40;
        double M_RATIO_MIN = 0;
        double M_RATIO_MAX = 1.;
    };
    static const Limits& limits ()
    {
        static Limits lims;
        return lims;
    }

public:

    static gun_model::Gun::GunType gun_type_impl(const Genome& genome);
    static Genome genome_impl();
};

class TableMaker
{
public:

    Mesh mesh;
    std::string dataPath;
    std::string gunName;
    bool clust;
    std::vector<std::string> names;
    std::vector<std::vector<std::vector<std::vector<double>>>> optimized_params;

    TableMaker (Mesh mesh_, std::string dataPath_, std::string gunName_, bool clust_)
        : mesh(mesh_), dataPath(dataPath_), gunName(gunName_), clust(clust_)
    {
        initialisation();
    }
    void initialisation();
    void make_params();
    void assemble(nlohmann::json& jsonObject) const;
};

} // namespace sig_optimizer
