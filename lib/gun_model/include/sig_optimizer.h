#pragma once

#include <random>
#include <ctime>
#include <fstream>
#include <functional>
#include <algorithm>
#include <gun_model.h>

namespace sig_optimizer
{

double randomDouble(double min, double max);
double norm_L2(const std::vector<double> &sig_1, const std::vector<double> &sig_2);
double norm_PEAKS(const std::vector<double> &sig_1, const std::vector<double> &sig_2);
double norm_HYBRID(const std::vector<double> &sig_1, const std::vector<double> &sig_2);

struct Mesh {
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

class Data {
public:
    Mesh mesh;
    std::string dataPath;
    std::vector<std::vector<std::vector<std::vector<double>>>> data_h_v_p;

    void parseData();
    void setMesh(Mesh mesh_) { mesh = mesh_; }
    void setDataPath(std::string dataPath_) { dataPath = dataPath_; }
};

struct Options {
    size_t POP_SIZE = 10;
    size_t POP_COUNT = 10;
    double MUTATION_RATE = 0.1;
    double MUTATION_SCALE = 0.5;
    double POPULATION_REJECTION = 0.05;
    size_t ELITE_GROUP_SIZE = 1;
};

struct Parameter {
    double min_value;
    double max_value;
    bool fixed;
    double value;

    Parameter(double min, double max, bool fix = false, double val = NAN) : min_value(min), max_value(max), fixed(fix), value(val) { }
    void choose(); // value initialisation
    void mutation(const Options &options); 
    Parameter operator+(const Parameter& other) const; // crossingover
};

struct Genome {
    std::vector<Parameter> params;
    double fitness;

    void choose(); // initialisation
    void mutation(const Options &options); 
    Genome operator+(const Genome& other) const; // crossingover
    template <typename Regime>
    void computeFitness(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double> &, const std::vector<double> &)> &norm); 
    template <typename Regime>
    void echoSignature(const gun_model::PhysParams &PP, const gun_model::MethodParams &MP, const Data& data, double h, double v, double p) const;
};

struct Population {
    std::vector<Genome> individs;

    template <typename Regime>
    void initialisation(size_t pop_size);
    void choose();
    void mutation(const Options &options);
    template <typename Regime>
    void computeFitness(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double> &, const std::vector<double> &)> &norm);
    void fitSort();
    const Genome& selection(const Options &options) const;
    template <typename Regime>
    Population geneticStep(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double> &, const std::vector<double> &)> &norm, const Options &options);
    void echoFit() const;
    void echoBest() const;
};

template <typename Regime>
class Genetic {
public:
    Options options;
    std::string dataPath;
    std::string gunName;
    Data data;
    gun_model::PhysParams PP;
    gun_model::MethodParams MP;
    std::vector<std::function<double(const std::vector<double>&, const std::vector<double>&)>> norms;
    Population population;

    Genetic(Options options_, Mesh mesh_, std::string dataPath_, std::string gunName_, const std::vector<std::function<double(const std::vector<double>&, const std::vector<double>&)>> &norms_) 
    : options(options_), dataPath(dataPath_), gunName(gunName_), norms(norms_)
    { data.setDataPath(dataPath + gunName + "/"); data.setMesh(mesh_); data.parseData(); PP.ref = 0; population.initialisation<Regime>(options.POP_SIZE); }
    std::string solve();
};

template <typename Derived>
class OptimisationRegime {
public:
    static gun_model::Gun::GunType gun_type(const Genome& genome, double h, double v, double p) { return Derived::gun_type_impl(genome, h, v, p); }
    static Genome genome() { return Derived::genome_impl(); }
};

class Default : public OptimisationRegime<Default> {
    struct Limits {
        double A_MIN = 0.01;
        double A_MAX = 0.5;
        double T_OPEN_MIN = 0.001;
        double T_OPEN_MAX = 0.01;
        double M_HC_MIN = 0;
        double M_HC_MAX = 100;
        double ALPHA_MIN = 0;
        double ALPHA_MAX = 10;
        double ALPHA_C_MIN = 0;
        double ALPHA_C_MAX = 10;
        double BY_MIN = 0;
        double BY_MAX = 20;
        double M_R_MIN = 0.5;
        double M_R_MAX = 10;
        double I_R_V_MIN = 0.1;
        double I_R_V_MAX = 1;
        double I_M_V_MIN = 0.1;
        double I_M_V_MAX = 1;
        double I_M_P_MIN = 0.1;
        double I_M_P_MAX = 1;
    };
    static const Limits& limits() { static Limits lims; return lims; }
public:
    static gun_model::Gun::GunType gun_type_impl(const Genome& genome, double h, double v, double p);
    static Genome genome_impl();
};

class Linear : public OptimisationRegime<Linear> {
    struct Limits {
        double A_MIN = 0.01;
        double A_MAX = 0.5;
        double A_MIN_P = -0.00001;
        double A_MAX_P = 0.00001;
        double A_MIN_V = -0.0001;
        double A_MAX_V = 0.0001;
        double A_MIN_H = -0.001;
        double A_MAX_H = 0.001;
        double T_OPEN_MIN = 0.001;
        double T_OPEN_MAX = 0.01;
        double T_OPEN_MIN_P = -0.000001;
        double T_OPEN_MAX_P = 0.000001;
        double T_OPEN_MIN_V = -0.00001;
        double T_OPEN_MAX_V = 0.00001;
        double T_OPEN_MIN_H = -0.0001;
        double T_OPEN_MAX_H = 0.0001;
        double M_HC_MIN = 0;
        double M_HC_MAX = 100;
        double M_HC_MIN_P = -0.001;
        double M_HC_MAX_P = 0.001;
        double M_HC_MIN_V = -0.01;
        double M_HC_MAX_V = 0.01;
        double M_HC_MIN_H = -0.1;
        double M_HC_MAX_H = 0.1;
        double ALPHA_MIN = 0;
        double ALPHA_MAX = 10;
        double ALPHA_MIN_P = -0.001;
        double ALPHA_MAX_P = 0.001;
        double ALPHA_MIN_V = -0.01;
        double ALPHA_MAX_V = 0.01;
        double ALPHA_MIN_H = -0.1;
        double ALPHA_MAX_H = 0.1;
        double ALPHA_C_MIN = 0;
        double ALPHA_C_MAX = 10;
        double ALPHA_C_MIN_P = -0.001;
        double ALPHA_C_MAX_P = 0.001;
        double ALPHA_C_MIN_V = -0.01;
        double ALPHA_C_MAX_V = 0.01;
        double ALPHA_C_MIN_H = -0.1;
        double ALPHA_C_MAX_H = 0.1;
        double BY_MIN = 0;
        double BY_MAX = 20;
        double BY_MIN_P = -0.001;
        double BY_MAX_P = 0.001;
        double BY_MIN_V = -0.01;
        double BY_MAX_V = 0.01;
        double BY_MIN_H = -0.1;
        double BY_MAX_H = 0.1;
        double M_R_MIN = 0.5;
        double M_R_MAX = 10;
        double M_R_MIN_P = -0.001;
        double M_R_MAX_P = 0.001;
        double M_R_MIN_V = -0.01;
        double M_R_MAX_V = 0.01;
        double M_R_MIN_H = -0.1;
        double M_R_MAX_H = 0.1;
        double I_R_V_MIN = 0.1;
        double I_R_V_MAX = 1;
        double I_M_V_MIN = 0.1;
        double I_M_V_MAX = 1;
        double I_M_P_MIN = 0.1;
        double I_M_P_MAX = 1;
    };
    static const Limits& limits() { static Limits lims; return lims; }
public:
    static gun_model::Gun::GunType gun_type_impl(const Genome& genome, double h, double v, double p);
    static Genome genome_impl();
};

} // namespace sig_optimizer


