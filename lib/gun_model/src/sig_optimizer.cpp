#include <sig_optimizer.h>

namespace sig_optimizer {

// template specification
template void Genome::computeFitness<OptimisationRegime<Default>>(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double> &, const std::vector<double> &)> &norm);
template void Genome::echoSignature<OptimisationRegime<Default>>(const gun_model::PhysParams &PP, const gun_model::MethodParams &MP, const Data& data, double h, double v, double p) const;
template void Population::initialisation<OptimisationRegime<Default>>(size_t pop_size);
template void Population::computeFitness<OptimisationRegime<Default>>(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double> &, const std::vector<double> &)> &norm);
template Population Population::geneticStep<OptimisationRegime<Default>>(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double> &, const std::vector<double> &)> &norm, const Options &options);
template class Genetic<OptimisationRegime<Default>>;
template void Genome::computeFitness<OptimisationRegime<Linear>>(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double> &, const std::vector<double> &)> &norm);
template void Genome::echoSignature<OptimisationRegime<Linear>>(const gun_model::PhysParams &PP, const gun_model::MethodParams &MP, const Data& data, double h, double v, double p) const;
template void Population::initialisation<OptimisationRegime<Linear>>(size_t pop_size);
template void Population::computeFitness<OptimisationRegime<Linear>>(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double> &, const std::vector<double> &)> &norm);
template Population Population::geneticStep<OptimisationRegime<Linear>>(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double> &, const std::vector<double> &)> &norm, const Options &options);
template class Genetic<OptimisationRegime<Linear>>;

double randomDouble(double min, double max) {
    if (min >= max)
        return min;
    static std::mt19937 rng(std::time(nullptr));
    std::uniform_real_distribution<double> dist(min, max);
    return dist(rng);
}

double norm_L2(const std::vector<double> &sig_1, const std::vector<double> &sig_2) {
    if (sig_1.size() != sig_2.size())
        throw std::logic_error("wrong signature size");
    double ret = 0;
    for (size_t i = 0; i < sig_1.size(); ++i)
        ret += std::pow(sig_1[i] - sig_2[i], 2);
    return ret;
}

double norm_PEAKS(const std::vector<double> &sig_1, const std::vector<double> &sig_2) {
    if (sig_1.size() != sig_2.size())
        throw std::logic_error("wrong signature size");
    double ret = 0;
    std::vector<double> times_1, times_2, values_1, values_2;
    auto start_iterator_1 = sig_1.begin();
    auto start_iterator_2 = sig_2.begin();
    for (size_t i = 0; i < 5; ++i) {
        auto max_1 = std::max_element(start_iterator_1, sig_1.end());
        auto max_2 = std::max_element(start_iterator_2, sig_2.end());
        if (max_1 != sig_1.end()) {
            times_1.push_back(static_cast<double>(std::distance(sig_1.begin(), max_1)));
            values_1.push_back(*max_1);
        }
        if (max_2 != sig_2.end()) {
            times_2.push_back(static_cast<double>(std::distance(sig_2.begin(), max_2)));
            values_2.push_back(*max_2);
        }
        start_iterator_1 = std::min_element(max_1, sig_1.end());
        start_iterator_2 = std::min_element(max_2, sig_2.end());
    }
    for (size_t i = 0; i < std::min(times_1.size(), times_2.size()); ++i) {
        ret += std::pow((times_1[i] - times_2[i]) * 0.001, 2);
        ret += std::pow(values_1[i] - values_2[i], 2);
    }
    if (times_1.size() != times_2.size()) 
        ret += 5;
    return ret;
}

double norm_HYBRID(const std::vector<double> &sig_1, const std::vector<double> &sig_2) {
    return norm_L2(sig_1, sig_2) + 20 * norm_PEAKS(sig_1, sig_2);
}

void Data::parseData() {
    data_h_v_p.resize((mesh.H_END - mesh.H_START) / mesh.H_STEP + 1);
    for (size_t i_h = 0; i_h < data_h_v_p.size(); ++i_h) {
        int h = mesh.H_START + mesh.H_STEP * i_h;
        data_h_v_p[i_h].resize((mesh.V_END - mesh.V_START) / mesh.V_STEP + 1);
        for (size_t i_v = 0; i_v < data_h_v_p[i_h].size(); ++i_v) {
            int v = mesh.V_START + mesh.V_STEP * i_v;
            data_h_v_p[i_h][i_v].resize((mesh.P_END - mesh.P_START) / mesh.P_STEP + 1);
            for (size_t i_p = 0; i_p < data_h_v_p[i_h][i_v].size(); ++i_p) {
                int p = mesh.P_START + mesh.P_STEP * i_p;
                data_h_v_p[i_h][i_v][i_p].resize(mesh.SIGN_SIZE);
                std::string filename = dataPath + "h" + std::to_string(h) + "_v" + std::to_string(v) + "_p" + std::to_string(p) + ".nsg";
                std::ifstream input;
                input.open(filename);
                if (!input.is_open()) {
                    throw std::runtime_error("File " + filename + " could not be opened");
                }
                std::string line; double val;
                for (size_t i = 0; i < 5; ++i)
                    std::getline(input, line);
                for (size_t i = 0; i < mesh.SIGN_SIZE; ++i) 
                    input >> val >> val >> data_h_v_p[i_h][i_v][i_p][i];
            }
        }
    }
}

void Parameter::choose() {   
    if (!fixed)
        value = randomDouble(min_value, max_value);
}

void Parameter::mutation(const Options &options) {
    if (std::isnan(value))
        throw std::logic_error("wrong mutation context");
    if (!fixed && randomDouble(0, 1) < options.MUTATION_RATE)
        value = value * (1 - options.MUTATION_SCALE) + options.MUTATION_SCALE * randomDouble(min_value, max_value);
}

Parameter Parameter::operator+(const Parameter& other) const {
    if (min_value != other.min_value || max_value != other.max_value || std::isnan(value) || std::isnan(other.value) || fixed != other.fixed)
        throw std::logic_error("wrong crossingover context");
    Parameter ret(min_value, max_value, fixed);
    if (fixed) {
        ret.value = value;
    } else {
        auto alpha = randomDouble(0, 1);
        ret.value = alpha * value + (1 - alpha) * other.value;
    }
    return ret;
}

void Genome::choose() {
    for (size_t i = 0; i < params.size(); ++i)
        params[i].choose();
    fitness = NAN;
}

void Genome::mutation(const Options &options) {
    if (std::isnan(fitness))
        for (size_t i = 0; i < params.size(); ++i)
            params[i].mutation(options);
}

Genome Genome::operator+(const Genome& other) const {
    if (params.size() != other.params.size())
        throw std::logic_error("wrong crossingover context");
    Genome ret;
    for (size_t i = 0; i < params.size(); ++i)
        ret.params.push_back(params[i] + other.params[i]);
    ret.fitness = NAN;
    return ret;
}

template <typename Regime>
void Genome::computeFitness(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double> &, const std::vector<double> &)> &norm) {
    if (!std::isnan(fitness))
        return;
    if (data.mesh.SIGN_SIZE != static_cast<size_t>(MP.repr_N))
        throw std::logic_error("signature sizes not coinside");
    int counter = 0;
    fitness = 0;
    for (size_t i_h = 0; i_h < data.data_h_v_p.size(); ++i_h) { 
        double h = data.mesh.H_START + data.mesh.H_STEP * i_h;
        for (size_t i_v = 0; i_v < data.data_h_v_p[i_h].size(); ++i_v) { 
            double v = data.mesh.V_START + data.mesh.V_STEP * i_v;
            for (size_t i_p = 0; i_p < data.data_h_v_p[i_h][i_v].size(); ++i_p) {
                double p = data.mesh.P_START + data.mesh.P_STEP * i_p;
                gun_model::Gun::GunInputs gi(p * 6894.76, v * 0.000016387064, 0, 0, 0, h);
                gun_model::Gun::GunPhysics gp;
                gun_model::Gun::GunType gt = Regime::gun_type(*this, h, v * 0.000016387064, p * 6894.76);
                gun_model::GunArray GA;
                GA.emplace_back(gt, gp, gi);
                gun_model::GunArraySolver<> solver(PP, MP, GA, gun_model::GunArraySolver<>::Cached::NONE);
                auto echo = solver.solve();
                if (echo == "Solving done well") {
                    counter++;
                    auto result = solver.getResult();
                    auto &signature = result.signatures[0];
                    auto &data_signature = data.data_h_v_p[i_h][i_v][i_p];
                    fitness += norm(signature, data_signature);
                } else {
                    fitness = -INFINITY;
                    return;
                }
            }
        }
    }
    fitness /= static_cast<double>(counter);
    fitness = -std::sqrt(fitness);
}

template <typename Regime>
void Genome::echoSignature(const gun_model::PhysParams &PP, const gun_model::MethodParams &MP, const Data& data, double h, double v, double p) const {
    size_t i_h = static_cast<size_t>((h - data.mesh.H_START) / data.mesh.H_STEP);
    size_t i_v = static_cast<size_t>((v - data.mesh.V_START) / data.mesh.V_STEP);
    size_t i_p = static_cast<size_t>((p - data.mesh.P_START) / data.mesh.P_STEP);
    gun_model::Gun::GunInputs gi(p * 6894.76, v * 0.000016387064, 0, 0, 0, h);
    gun_model::Gun::GunPhysics gp;
    gun_model::Gun::GunType gt = Regime::gun_type(*this, h, v * 0.000016387064, p * 6894.76);
    gun_model::GunArray GA;
    GA.emplace_back(gt, gp, gi);
    gun_model::GunArraySolver<> solver(PP, MP, GA, gun_model::GunArraySolver<>::Cached::NONE);
    auto echo = solver.solve();
    if (echo == "Solving done well") {
        auto result = solver.getResult();
        const auto &signature = result.signatures[0];
        const auto &data_signature = data.data_h_v_p[i_h][i_v][i_p];
        std::cout << "solver's signature: ";
        for (size_t i = 0; i < data.mesh.SIGN_SIZE; ++i)
            std::cout << signature[i] << ", ";
        std::cout << "\ngundalf's signature: ";
        for (size_t i = 0; i < data.mesh.SIGN_SIZE; ++i)
            std::cout << data_signature[i] << ", ";
    }
}

template <typename Regime>
void Population::initialisation(size_t pop_size) {
    for (size_t i = 0; i < pop_size; ++i) 
        individs.push_back(Regime::genome());
}

void Population::choose() {
    for (auto &genome : individs) 
        genome.choose();
}

void Population::mutation(const Options &options) {
    for (auto &genome : individs) 
        genome.mutation(options);
}

template <typename Regime>
void Population::computeFitness(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double> &, const std::vector<double> &)> &norm) {
    for (auto &genome : individs) 
        genome.computeFitness<Regime>(PP, MP, data, norm);
}

void Population::fitSort() {
    std::sort(individs.begin(), individs.end(), [](const Genome &gen_1, const Genome &gen_2) { return gen_1.fitness < gen_2.fitness; });
}

const Genome& Population::selection(const Options &options) const {
    double total_rank = individs.size() * (individs.size() + 1) / 2;
    double rank = randomDouble(total_rank * options.POPULATION_REJECTION, total_rank);
    size_t i = 0;
    double sum = 1;
    while (sum < rank) { // селекция рейнджированием
        i++;
        sum += i + 1;
    }
    if (i < individs.size())
        return individs[i];
    return individs.back();
}

template <typename Regime>
Population Population::geneticStep(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double> &, const std::vector<double> &)> &norm, const Options &options) {
    Population next_pop;
    for (size_t i = 0; i < individs.size() - options.ELITE_GROUP_SIZE; ++i) {
        auto par_1 = selection(options);
        auto par_2 = selection(options);
        auto child = par_1 + par_2;
        next_pop.individs.push_back(child);
    }
    for (size_t i = 0; i < options.ELITE_GROUP_SIZE; ++i)
        next_pop.individs.push_back(individs[individs.size() - 1 - i]);
    next_pop.mutation(options);
    next_pop.computeFitness<Regime>(PP, MP, data, norm);
    next_pop.fitSort();
    return next_pop;
}

void Population::echoFit() const {
    std::cout << "fitness: ";
    for (const auto &ind : individs)
        std::cout << ind.fitness << " ";
    std::cout << "\n";
}

void Population::echoBest() const {
    std::cout << "best genome: ";
    for (const auto &param : individs.back().params)
        std::cout << param.value << " ";
    std::cout << "\n";
}

template <typename Regime>
std::string Genetic<Regime>::solve() {
    try {    
        std::cout << "population 1:\n";
        int norm_num = 0;
        population.choose();
        population.computeFitness<Regime>(PP, MP, data, norms[norm_num]);
        population.fitSort();
        population.echoFit();
        population.echoBest();
        for (size_t i = 2; i <= options.POP_COUNT; ++i) {
            std::cout << "population " << i << ":\n";
            norm_num = (norm_num + 1) % norms.size();
            population = population.geneticStep<Regime>(PP, MP, data, norms[norm_num], options);
            population.echoFit();
            population.echoBest();
        }
    } catch (const std::exception &e) {
        return e.what();
    } 
    return std::string();
}

// дефолтный оптимизатор (без зависимостей от входных параметров h, v, p)
gun_model::Gun::GunType Default::gun_type_impl(const Genome& genome, double /*h*/, double /*v*/, double /*p*/) {
    if (genome.params.size() != 10)
        throw std::logic_error("wrong input genome");
    for (size_t i = 0; i < genome.params.size(); ++i)
        if (std::isnan(genome.params[i].value))
            throw std::logic_error("wrong input genome");
    
    gun_model::Gun::GunType gt;
    gt.A = genome.params[0].value;
    gt.t_open = genome.params[1].value;
    gt.M_hc = genome.params[2].value;
    gt.alpha = genome.params[3].value;
    gt.alpha_c = genome.params[4].value;
    gt.by = genome.params[5].value;
    gt.M_r = genome.params[6].value;
    gt.i_r_v = genome.params[7].value;
    gt.i_m_v = genome.params[8].value;
    gt.i_m_p = genome.params[9].value;

    return gt;
}

Genome Default::genome_impl() {
    Genome ret;
    const auto& lims = Default::limits();
    ret.params.emplace_back(lims.A_MIN, lims.A_MAX);
    ret.params.emplace_back(lims.T_OPEN_MIN, lims.T_OPEN_MAX);
    ret.params.emplace_back(lims.M_HC_MIN, lims.M_HC_MAX);
    ret.params.emplace_back(lims.ALPHA_MIN, lims.ALPHA_MAX);
    ret.params.emplace_back(lims.ALPHA_C_MIN, lims.ALPHA_C_MAX);
    ret.params.emplace_back(lims.BY_MIN, lims.BY_MAX);
    ret.params.emplace_back(lims.M_R_MIN, lims.M_R_MAX);
    ret.params.emplace_back(lims.I_R_V_MIN, lims.I_R_V_MAX);
    ret.params.emplace_back(lims.I_M_V_MIN, lims.I_M_V_MAX);
    ret.params.emplace_back(lims.I_M_P_MIN, lims.I_M_P_MAX);
    ret.fitness = NAN;

    return ret;
}

// линейный оптимизатор 
gun_model::Gun::GunType Linear::gun_type_impl(const Genome& genome, double h, double v, double p) {
    if (genome.params.size() != 31)
        throw std::logic_error("wrong input genome");
    for (size_t i = 0; i < genome.params.size(); ++i)
        if (std::isnan(genome.params[i].value))
            throw std::logic_error("wrong input genome");
    
    const auto& lims = Linear::limits();
    gun_model::Gun::GunType gt;
    gt.A = genome.params[0].value + genome.params[1].value * p + genome.params[2].value * v + genome.params[3].value * h;
    if (gt.A < lims.A_MIN) gt.A = lims.A_MIN;
    if (gt.A > lims.A_MAX) gt.A = lims.A_MAX;

    gt.t_open = genome.params[4].value + genome.params[5].value * p + genome.params[6].value * v + genome.params[7].value * h;
    if (gt.t_open < lims.T_OPEN_MIN) gt.t_open = lims.T_OPEN_MIN;
    if (gt.t_open > lims.T_OPEN_MAX) gt.t_open = lims.T_OPEN_MAX;

    gt.M_hc = genome.params[8].value + genome.params[9].value * p + genome.params[10].value * v + genome.params[11].value * h;
    if (gt.M_hc < lims.M_HC_MIN) gt.M_hc = lims.M_HC_MIN;
    if (gt.M_hc > lims.M_HC_MAX) gt.M_hc = lims.M_HC_MAX;

    gt.alpha = genome.params[12].value + genome.params[13].value * p + genome.params[14].value * v + genome.params[15].value * h;
    if (gt.alpha < lims.ALPHA_MIN) gt.alpha = lims.ALPHA_MIN;
    if (gt.alpha > lims.ALPHA_MAX) gt.alpha = lims.ALPHA_MAX;

    gt.alpha_c = genome.params[16].value + genome.params[17].value * p + genome.params[18].value * v + genome.params[19].value * h;
    if (gt.alpha_c < lims.ALPHA_C_MIN) gt.alpha_c = lims.ALPHA_C_MIN;
    if (gt.alpha_c > lims.ALPHA_C_MAX) gt.alpha_c = lims.ALPHA_C_MAX;

    gt.by = genome.params[20].value + genome.params[21].value * p + genome.params[22].value * v + genome.params[23].value * h;
    if (gt.by < lims.BY_MIN) gt.by = lims.BY_MIN;
    if (gt.by > lims.BY_MAX) gt.by = lims.BY_MAX;

    gt.M_r = genome.params[24].value + genome.params[25].value * p + genome.params[26].value * v + genome.params[27].value * h;
    if (gt.M_r < lims.M_R_MIN) gt.M_r = lims.M_R_MIN;
    if (gt.M_r > lims.M_R_MAX) gt.M_r = lims.M_R_MAX;
    
    gt.i_r_v = genome.params[28].value;
    gt.i_m_v = genome.params[29].value;
    gt.i_m_p = genome.params[30].value;

    return gt;
}

Genome Linear::genome_impl() {
    Genome ret;
    const auto& lims = Linear::limits();
    ret.params.emplace_back(lims.A_MIN, lims.A_MAX);
    ret.params.emplace_back(lims.A_MIN_P, lims.A_MAX_P);
    ret.params.emplace_back(lims.A_MIN_V, lims.A_MAX_V);
    ret.params.emplace_back(lims.A_MIN_H, lims.A_MAX_H);

    ret.params.emplace_back(lims.T_OPEN_MIN, lims.T_OPEN_MAX);
    ret.params.emplace_back(lims.T_OPEN_MIN_P, lims.T_OPEN_MAX_P);
    ret.params.emplace_back(lims.T_OPEN_MIN_V, lims.T_OPEN_MAX_V);
    ret.params.emplace_back(lims.T_OPEN_MIN_H, lims.T_OPEN_MAX_H);

    ret.params.emplace_back(lims.M_HC_MIN, lims.M_HC_MAX);
    ret.params.emplace_back(lims.M_HC_MIN_P, lims.M_HC_MAX_P);
    ret.params.emplace_back(lims.M_HC_MIN_V, lims.M_HC_MAX_V);
    ret.params.emplace_back(lims.M_HC_MIN_H, lims.M_HC_MAX_H);

    ret.params.emplace_back(lims.ALPHA_MIN, lims.ALPHA_MAX);
    ret.params.emplace_back(lims.ALPHA_MIN_P, lims.ALPHA_MAX_P);
    ret.params.emplace_back(lims.ALPHA_MIN_V, lims.ALPHA_MAX_V);
    ret.params.emplace_back(lims.ALPHA_MIN_H, lims.ALPHA_MAX_H);

    ret.params.emplace_back(lims.ALPHA_C_MIN, lims.ALPHA_C_MAX);
    ret.params.emplace_back(lims.ALPHA_C_MIN_P, lims.ALPHA_C_MAX_P);
    ret.params.emplace_back(lims.ALPHA_C_MIN_V, lims.ALPHA_C_MAX_V);
    ret.params.emplace_back(lims.ALPHA_C_MIN_H, lims.ALPHA_C_MAX_H);

    ret.params.emplace_back(lims.BY_MIN, lims.BY_MAX);
    ret.params.emplace_back(lims.BY_MIN_P, lims.BY_MAX_P);
    ret.params.emplace_back(lims.BY_MIN_V, lims.BY_MAX_V);
    ret.params.emplace_back(lims.BY_MIN_H, lims.BY_MAX_H);

    ret.params.emplace_back(lims.M_R_MIN, lims.M_R_MAX);
    ret.params.emplace_back(lims.M_R_MIN_P, lims.M_R_MAX_P);
    ret.params.emplace_back(lims.M_R_MIN_V, lims.M_R_MAX_V);
    ret.params.emplace_back(lims.M_R_MIN_H, lims.M_R_MAX_H);

    ret.params.emplace_back(lims.I_R_V_MIN, lims.I_R_V_MAX);
    ret.params.emplace_back(lims.I_M_V_MIN, lims.I_M_V_MAX);
    ret.params.emplace_back(lims.I_M_P_MIN, lims.I_M_P_MAX);
    ret.fitness = NAN;

    return ret;
}

} // namespace sig_optimizer
