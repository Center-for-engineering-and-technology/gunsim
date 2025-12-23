#include <sig_optimizer.h>

namespace sig_optimizer
{

// template specification
template void Genome::computeFitness<OptimisationRegime<Default>>(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm);
template void Genome::computeFitnessClust<OptimisationRegime<Default>>(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm);
template void Genome::echoSignature<OptimisationRegime<Default>>(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, double h, double v, double p) const;
template void Population::initialisation<OptimisationRegime<Default>>(size_t pop_size);
template void Population::computeFitness<OptimisationRegime<Default>>(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm);
template void Population::computeFitnessClust<OptimisationRegime<Default>>(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm);
template Population Population::geneticStep<OptimisationRegime<Default>>(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm, const Options& options);
template Population Population::geneticStepClust<OptimisationRegime<Default>>(const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm, const Options& options);
template class Algorithm<OptimisationRegime<Default>>;

double randomDouble (double min, double max)
{
    if (min >= max)
        return min;
    static std::mt19937 rng(static_cast<int>(std::time(nullptr)));
    std::uniform_real_distribution<double> dist(min, max);
    return dist(rng);
}

double norm_L2 (const std::vector<double>& sig_1, const std::vector<double>& sig_2)
{
    double s = 0, s_1 = 0, s_2 = 0;
    for (size_t i = 0; i < sig_1.size(); ++i)
    {
        s += std::pow(sig_1[i] - sig_2[i], 2);
        s_1 += std::pow(sig_1[i], 2);
        s_2 += std::pow(sig_2[i], 2);
    }
    s = std::sqrt(s);
    s_1 = std::sqrt(s_1);
    s_2 = std::sqrt(s_2);
    return 2 * s / (s_1 + s_2);
}

double norm_TECHNICAL (const std::vector<double>& sig_1, const std::vector<double>& sig_2)
{
    auto tech = [] (const std::vector<double>& sig) -> std::array<double, 4>
    {
        auto [min, max] = std::minmax_element(sig.begin(), sig.end());
        auto bubble_max = std::max_element(min, sig.end());
        auto bubble_min = std::min_element(bubble_max, sig.end());
        return {*max - *min, *max, (*max - *min) / (*bubble_max - *bubble_min), static_cast<double>(distance(sig.begin(), bubble_max)) / static_cast<double>(sig.size())};
    };
    auto [peak_to_peak_1, zero_to_peak_1, primary_to_bubble_1, bubble_period_1] = tech(sig_1);
    auto [peak_to_peak_2, zero_to_peak_2, primary_to_bubble_2, bubble_period_2] = tech(sig_2);
    return std::sqrt(0.5 * (std::pow(peak_to_peak_1 - peak_to_peak_2, 2) / (std::pow(peak_to_peak_1, 2) + std::pow(peak_to_peak_2, 2)) + std::pow(zero_to_peak_1 - zero_to_peak_2, 2) / (std::pow(zero_to_peak_1, 2) + std::pow(zero_to_peak_2, 2)) + std::pow(primary_to_bubble_1 - primary_to_bubble_2, 2) / (std::pow(primary_to_bubble_1, 2) + std::pow(primary_to_bubble_2, 2)) + std::pow(bubble_period_1 - bubble_period_2, 2) / (std::pow(bubble_period_1, 2) + std::pow(bubble_period_2, 2))));
}

double norm_TECHNICAL_L2 (const std::vector<double>& sig_1, const std::vector<double>& sig_2)
{
    auto tech = [] (const std::vector<double>& sig) -> std::array<double, 4>
    {
        auto [min, max] = std::minmax_element(sig.begin(), sig.end());
        auto bubble_max = std::max_element(min, sig.end());
        auto bubble_min = std::min_element(bubble_max, sig.end());
        return {*max - *min, *max, (*max - *min) / (*bubble_max - *bubble_min), static_cast<double>(distance(sig.begin(), bubble_max))};
    };
    auto [peak_to_peak_1, zero_to_peak_1, primary_to_bubble_1, bubble_period_1] = tech(sig_1);
    auto [peak_to_peak_2, zero_to_peak_2, primary_to_bubble_2, bubble_period_2] = tech(sig_2);
    double s = 0, s_1 = 0, s_2 = 0;
    for (size_t i = 0; i < sig_1.size(); ++i)
    {
        s += std::pow(sig_1[i] - sig_2[i], 2);
        s_1 += std::pow(sig_1[i], 2);
        s_2 += std::pow(sig_2[i], 2);
    }
    s = std::sqrt(s);
    s_1 = std::sqrt(s_1);
    s_2 = std::sqrt(s_2);
    double q = s / (s_1 + s_2);
    return 0.4 * (std::abs(peak_to_peak_1 - peak_to_peak_2) / (peak_to_peak_1 + peak_to_peak_2) + std::abs(zero_to_peak_1 - zero_to_peak_2) / (zero_to_peak_1 + zero_to_peak_2) + std::abs(primary_to_bubble_1 - primary_to_bubble_2) / (primary_to_bubble_1 + primary_to_bubble_2) + std::abs(bubble_period_1 - bubble_period_2) / (bubble_period_1 + bubble_period_2) + 0.1 * q);
}

void Data::parseData (bool ref)
{
    data_h_v_p.resize((mesh.H_END - mesh.H_START) / mesh.H_STEP + 1);
    for (size_t i_h = 0; i_h < data_h_v_p.size(); ++i_h)
    {
        int h = static_cast<int>(mesh.H_START + mesh.H_STEP * i_h);
        data_h_v_p[i_h].resize((mesh.V_END - mesh.V_START) / mesh.V_STEP + 1);
        for (size_t i_v = 0; i_v < data_h_v_p[i_h].size(); ++i_v)
        {
            int v = static_cast<int>(mesh.V_START + mesh.V_STEP * i_v);
            data_h_v_p[i_h][i_v].resize((mesh.P_END - mesh.P_START) / mesh.P_STEP + 1);
            for (size_t i_p = 0; i_p < data_h_v_p[i_h][i_v].size(); ++i_p)
            {
                int p = static_cast<int>(mesh.P_START + mesh.P_STEP * i_p);
                data_h_v_p[i_h][i_v][i_p].resize(mesh.SIGN_SIZE);
                std::string filename = dataPath + "h" + std::to_string(h) + "_v" + std::to_string(v) + "_p" + std::to_string(p) + ".nsg";
                std::ifstream input;
                input.open(filename);
                if (!input.is_open())
                {
                    throw std::runtime_error("File " + filename + " could not be opened");
                }
                std::string line;
                double val;
                for (size_t i = 0; i < 5; ++i)
                    std::getline(input, line);
                for (size_t i = 0; i < mesh.SIGN_SIZE; ++i)
                    input >> val >> val >> data_h_v_p[i_h][i_v][i_p][i];
                if (ref)
                {
                    auto tmp_vec = data_h_v_p[i_h][i_v][i_p];
                    int j = static_cast<int>(std::floor(2 * h / 1496. / 0.0005));
                    for (int i = 0; i < static_cast<int>(mesh.SIGN_SIZE); ++i)
                    {
                        if (i - j >= 0)
                            data_h_v_p[i_h][i_v][i_p][i] -= tmp_vec[i - j];
                    }
                }
                input.close();
            }
        }
    }
}

void Data::parseDataClust (bool ref)
{
    data_h_v_p_clust.resize((mesh.H_END - mesh.H_START) / mesh.H_STEP + 1);
    for (size_t i_h = 0; i_h < data_h_v_p_clust.size(); ++i_h)
    {
        int h = static_cast<int>(mesh.H_START + mesh.H_STEP * i_h);
        data_h_v_p_clust[i_h].resize((mesh.V_END - mesh.V_START) / mesh.V_STEP + 1);
        for (size_t i_v = 0; i_v < data_h_v_p_clust[i_h].size(); ++i_v)
        {
            int v = static_cast<int>(mesh.V_START + mesh.V_STEP * i_v);
            data_h_v_p_clust[i_h][i_v].resize((mesh.P_END - mesh.P_START) / mesh.P_STEP + 1);
            for (size_t i_p = 0; i_p < data_h_v_p_clust[i_h][i_v].size(); ++i_p)
            {
                int p = static_cast<int>(mesh.P_START + mesh.P_STEP * i_p);
                data_h_v_p_clust[i_h][i_v][i_p].resize(mesh.SIGN_SIZE, 0);
                std::string filename_clust = dataPath + "h" + std::to_string(h) + "_v" + std::to_string(v) + "_p" + std::to_string(p) + "_clust.nsg";
                std::ifstream input_clust;
                input_clust.open(filename_clust);
                if (!input_clust.is_open())
                {
                    throw std::runtime_error("File " + filename_clust + " could not be opened");
                }
                std::string line;
                double val;
                for (size_t i = 0; i < 5; ++i)
                    std::getline(input_clust, line);
                for (size_t i = 0; i < mesh.SIGN_SIZE; ++i)
                {
                    input_clust >> val >> val >> val;
                    data_h_v_p_clust[i_h][i_v][i_p][i] += 2 * val;
                }
                if (ref)
                {
                    auto tmp_vec = data_h_v_p_clust[i_h][i_v][i_p];
                    int j = static_cast<int>(std::floor(2 * h / 1496. / 0.0005));
                    for (int i = 0; i < static_cast<int>(mesh.SIGN_SIZE); ++i)
                    {
                        if (i - j >= 0)
                            data_h_v_p_clust[i_h][i_v][i_p][i] -= tmp_vec[i - j];
                    }
                }
                input_clust.close();
            }
        }
    }
}

void Parameter::choose ()
{
    if (!fixed)
        value = randomDouble(min_value, max_value);
}

void Parameter::mutation (const Options& options)
{
    if (std::isnan(value))
        throw std::logic_error("wrong mutation context");
    if (!fixed && randomDouble(0, 1) < options.MUTATION_RATE)
        value = value * (1 - options.MUTATION_SCALE) + options.MUTATION_SCALE * randomDouble(min_value, max_value);
}

Parameter Parameter::operator+(const Parameter& other) const
{
    if (min_value != other.min_value || max_value != other.max_value || std::isnan(value) || std::isnan(other.value) || fixed != other.fixed)
        throw std::logic_error("wrong crossingover context");
    Parameter ret(min_value, max_value, fixed);
    if (fixed)
    {
        ret.value = value;
    }
    else
    {
        auto alpha = randomDouble(0, 1);
        ret.value = alpha * value + (1 - alpha) * other.value;
    }
    return ret;
}

void Genome::choose ()
{
    for (size_t i = 0; i < params.size(); ++i)
        params[i].choose();
    fitness = NAN;
}

void Genome::mutation (const Options& options)
{
    if (std::isnan(fitness))
        for (size_t i = 0; i < params.size(); ++i)
            params[i].mutation(options);
}

Genome Genome::operator+(const Genome& other) const
{
    if (params.size() != other.params.size())
        throw std::logic_error("wrong crossingover context");
    Genome ret;
    for (size_t i = 0; i < params.size(); ++i)
        ret.params.push_back(params[i] + other.params[i]);
    ret.fitness = NAN;
    return ret;
}

template<typename Regime>
void Genome::computeFitness (const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm)
{
    if (!std::isnan(fitness))
        return;
    int counter = 0;
    fitness = 0;
    for (size_t i_h = 0; i_h < data.data_h_v_p.size(); ++i_h)
    {
        double h = static_cast<double>(data.mesh.H_START + data.mesh.H_STEP * i_h);
        for (size_t i_v = 0; i_v < data.data_h_v_p[i_h].size(); ++i_v)
        {
            double v = static_cast<double>(data.mesh.V_START + data.mesh.V_STEP * i_v);
            for (size_t i_p = 0; i_p < data.data_h_v_p[i_h][i_v].size(); ++i_p)
            {
                double p = static_cast<double>(data.mesh.P_START + data.mesh.P_STEP * i_p);
                gun_model::Gun::GunInputs gi(p * 6894.76, v * 0.000016387064, 0, 0, 0, h);
                gun_model::Gun::GunPhysics gp;
                gun_model::Gun::GunType gt = Regime::gun_type(*this);
                gun_model::GunArray GA;
                GA.emplace_back(gt, gp, gi);
                gun_model::GunArraySolver<> solver(PP, MP, GA, gun_model::GunArraySolver<>::Cached::NONE, gun_model::GunArraySolver<>::InteractionModel::DEFAULT_INTERACTION);
                auto echo = solver.solve();
                if (echo == "Solving done well")
                {
                    counter++;
                    auto result = solver.getResult();
                    auto signature = result.signatures[0];
                    int j = static_cast<int>(std::floor(2 * h / PP.c / MP.sampleMax * MP.repr_N));
                    for (int i = 0; i < static_cast<int>(signature.size()); ++i)
                    {
                        if (i - j >= 0)
                            signature[i] += PP.ref * result.signatures[0][i - j];
                    }
                    auto& data_signature = data.data_h_v_p[i_h][i_v][i_p];
                    fitness += norm(signature, data_signature);
                }
                else
                {
                    fitness = -INFINITY;
                    return;
                }
            }
        }
    }
    fitness /= static_cast<double>(counter);
    fitness *= -1;
}

template<typename Regime>
void Genome::computeFitnessClust (const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm)
{
    if (!std::isnan(fitness))
        return;
    int counter = 0;
    fitness = 0;
    for (size_t i_h = 0; i_h < data.data_h_v_p.size(); ++i_h)
    {
        double h = static_cast<double>(data.mesh.H_START + data.mesh.H_STEP * i_h);
        for (size_t i_v = 0; i_v < data.data_h_v_p[i_h].size(); ++i_v)
        {
            double v = static_cast<double>(data.mesh.V_START + data.mesh.V_STEP * i_v);
            for (size_t i_p = 0; i_p < data.data_h_v_p[i_h][i_v].size(); ++i_p)
            {
                double p = static_cast<double>(data.mesh.P_START + data.mesh.P_STEP * i_p);
                gun_model::Gun::GunInputs gi(p * 6894.76, v * 0.000016387064, 0, 0, 0, h);
                gun_model::Gun::GunPhysics gp;
                gun_model::Gun::GunType gt = Regime::gun_type(*this);
                gun_model::GunArray GA;
                GA.emplace_back(gt, gp, gi);
                gun_model::GunArraySolver<> solver(PP, MP, GA, gun_model::GunArraySolver<>::Cached::NONE, gun_model::GunArraySolver<>::InteractionModel::DEFAULT_INTERACTION);
                auto echo = solver.solve();
                if (echo == "Solving done well")
                {
                    counter++;
                    auto result = solver.getResult();
                    auto signature = result.signatures[0];
                    int j = static_cast<int>(std::floor(2 * h / PP.c / MP.sampleMax * MP.repr_N));
                    for (int i = 0; i < static_cast<int>(signature.size()); ++i)
                    {
                        if (i - j >= 0)
                            signature[i] += PP.ref * result.signatures[0][i - j];
                    }
                    auto& data_signature = data.data_h_v_p[i_h][i_v][i_p];
                    fitness += norm(signature, data_signature);
                }
                else
                {
                    fitness = -INFINITY;
                    return;
                }
                gi.x = 1.;
                GA.emplace_back(gt, gp, gi);
                gun_model::GunArraySolver<> solver_clust(PP, MP, GA, gun_model::GunArraySolver<>::Cached::NONE, gun_model::GunArraySolver<>::InteractionModel::DEFAULT_INTERACTION);
                echo = solver_clust.solve();
                if (echo == "Solving done well")
                {
                    counter++;
                    auto result = solver_clust.getResult();
                    std::vector<double> signature(result.signatures[0].size(), 0);
                    int j = static_cast<int>(std::floor(2 * h / PP.c / MP.sampleMax * MP.repr_N));
                    for (size_t k = 0; k < 2; ++k)
                        for (int i = 0; i < static_cast<int>(signature.size()); ++i)
                        {
                            signature[i] += result.signatures[k][i];
                            if (i - j >= 0)
                                signature[i] += PP.ref * result.signatures[k][i - j];
                        }
                    auto& data_signature = data.data_h_v_p_clust[i_h][i_v][i_p];
                    fitness += norm(signature, data_signature);
                }
                else
                {
                    fitness = -INFINITY;
                    return;
                }
            }
        }
    }
    fitness /= static_cast<double>(counter);
    fitness *= -1;
}

template<typename Regime>
void Genome::echoSignature (const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, double h, double v, double p) const
{
    size_t i_h = static_cast<size_t>((h - data.mesh.H_START) / data.mesh.H_STEP);
    size_t i_v = static_cast<size_t>((v - data.mesh.V_START) / data.mesh.V_STEP);
    size_t i_p = static_cast<size_t>((p - data.mesh.P_START) / data.mesh.P_STEP);
    gun_model::Gun::GunInputs gi(p * 6894.76, v * 0.000016387064, 0, 0, 0, h);
    gun_model::Gun::GunPhysics gp;
    gun_model::Gun::GunType gt = Regime::gun_type(*this);
    gun_model::GunArray GA;
    GA.emplace_back(gt, gp, gi);
    gun_model::GunArraySolver<> solver(PP, MP, GA, gun_model::GunArraySolver<>::Cached::NONE, gun_model::GunArraySolver<>::InteractionModel::DEFAULT_INTERACTION);
    auto echo = solver.solve();
    if (echo == "Solving done well")
    {
        auto result = solver.getResult();
        auto signature = result.signatures[0];
        int j = static_cast<int>(std::floor(2 * h / PP.c / MP.sampleMax * MP.repr_N));
        for (int i = 0; i < static_cast<int>(signature.size()); ++i)
        {
            if (i - j >= 0)
                signature[i] += PP.ref * result.signatures[0][i - j];
        }
        auto& data_signature = data.data_h_v_p[i_h][i_v][i_p];
        std::cout << "solver's signature: ";
        for (size_t i = 0; i < signature.size(); ++i)
            std::cout << signature[i] << ", ";
        std::cout << "\n\ngundalf's signature: ";
        for (size_t i = 0; i < data.mesh.SIGN_SIZE; ++i)
            std::cout << data_signature[i] << ", ";
    }
}

template<typename Regime>
void Population::initialisation (size_t pop_size)
{
    for (size_t i = 0; i < pop_size; ++i)
        individs.push_back(Regime::genome());
}

void Population::choose ()
{
    for (auto& genome: individs)
        genome.choose();
}

void Population::mutation (const Options& options)
{
    for (auto& genome: individs)
        genome.mutation(options);
}

template<typename Regime>
void Population::computeFitness (const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm)
{
    // #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(individs.size()); ++i)
        individs[i].computeFitness<Regime>(PP, MP, data, norm);
}

template<typename Regime>
void Population::computeFitnessClust (const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm)
{
    // #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(individs.size()); ++i)
        individs[i].computeFitnessClust<Regime>(PP, MP, data, norm);
}

void Population::fitSort ()
{
    std::sort(individs.begin(), individs.end(), [] (const Genome& gen_1, const Genome& gen_2)
              {
                  return gen_1.fitness < gen_2.fitness;
              });
}

const Genome& Population::selection (const Options& options) const
{
    double total_rank = static_cast<double>(individs.size() * (individs.size() + 1) / 2);
    double rank = randomDouble(total_rank * options.POPULATION_REJECTION, total_rank);
    size_t i = 0;
    double sum = 1;
    while (sum < rank)
    { // селекция рейнджированием
        i++;
        sum += i + 1;
    }
    if (i < individs.size())
        return individs[i];
    return individs.back();
}

template<typename Regime>
Population Population::geneticStep (const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm, const Options& options)
{
    Population next_pop;
    for (size_t i = 0; i < individs.size() - options.ELITE_GROUP_SIZE; ++i)
    {
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

template<typename Regime>
Population Population::geneticStepClust (const gun_model::PhysParams& PP, const gun_model::MethodParams& MP, const Data& data, const std::function<double(const std::vector<double>&, const std::vector<double>&)>& norm, const Options& options)
{
    Population next_pop;
    for (size_t i = 0; i < individs.size() - options.ELITE_GROUP_SIZE; ++i)
    {
        auto par_1 = selection(options);
        auto par_2 = selection(options);
        auto child = par_1 + par_2;
        next_pop.individs.push_back(child);
    }
    for (size_t i = 0; i < options.ELITE_GROUP_SIZE; ++i)
        next_pop.individs.push_back(individs[individs.size() - 1 - i]);
    next_pop.mutation(options);
    next_pop.computeFitnessClust<Regime>(PP, MP, data, norm);
    next_pop.fitSort();
    return next_pop;
}

void Population::echoFit () const
{
    std::cout << "fitness: ";
    for (auto& ind: individs)
        std::cout << ind.fitness << " ";
    std::cout << "\n";
}

void Population::echoBest () const
{
    std::cout << "best genome: ";
    for (auto& param: individs.back().params)
        std::cout << param.value << " ";
    std::cout << "\n";
}

template<typename Regime>
std::string Algorithm<Regime>::genetic ()
{
    try
    {
        // std::cout << "population 1:\n";
        population.choose();
        population.computeFitness<Regime>(PP, MP, data, norm);
        population.fitSort();
        // population.echoFit();
        // population.echoBest();
        for (size_t i = 2; i <= options.POP_COUNT; ++i)
        {
            // std::cout << "population " << i << ":\n";
            population = population.geneticStep<Regime>(PP, MP, data, norm, options);
            // population.echoFit();
            // population.echoBest();
            // gradient(population.individs.back());
        }
    }
    catch (const std::exception& e)
    {
        return e.what();
    }
    return std::string();
}

template<typename Regime>
std::string Algorithm<Regime>::geneticClust ()
{
    try
    {
        // std::cout << "population 1:\n";
        population.choose();
        population.computeFitnessClust<Regime>(PP, MP, data, norm);
        population.fitSort();
        // population.echoFit();
        // population.echoBest();
        for (size_t i = 2; i <= options.POP_COUNT; ++i)
        {
            // std::cout << "population " << i << ":\n";
            population = population.geneticStepClust<Regime>(PP, MP, data, norm, options);
            // population.echoFit();
            // population.echoBest();
            // gradient(population.individs.back());
        }
    }
    catch (const std::exception& e)
    {
        return e.what();
    }
    return std::string();
}

template<typename Regime>
std::string Algorithm<Regime>::gradient (Genome& start) const
{
    try
    {
        start.computeFitness<Regime>(PP, MP, data, norm);
        for (size_t k = 0; k < 20; ++k)
        {
            double h = 0.01;
            std::vector<double> grad(start.params.size());
            for (size_t i = 0; i < grad.size(); ++i)
            {
                Genome moved_forward = start;
                Genome moved_back = start;
                moved_forward.fitness = NAN;
                moved_back.fitness = NAN;
                moved_forward.params[i].value += start.params[i].step;
                moved_back.params[i].value -= start.params[i].step;
                moved_forward.computeFitness<Regime>(PP, MP, data, norm);
                moved_back.computeFitness<Regime>(PP, MP, data, norm);
                grad[i] = (moved_forward.fitness - moved_back.fitness) / start.params[i].step / 2.;
            }

            Genome next;

            for (size_t j = 0; j < 20; ++j)
            {
                next = start;
                next.fitness = NAN;
                for (size_t i = 0; i < grad.size(); ++i)
                {
                    next.params[i].value += h * grad[i];
                }
                next.computeFitness<Regime>(PP, MP, data, norm);
                if (next.fitness > start.fitness)
                {
                    // std::cout << "gradient iteartions: " << j + 1 << "\n";
                    break;
                }
                else
                    h *= 0.5;
            }
            if (next.fitness > start.fitness)
            {
                start = next;
                // std::cout << "new best fitness: " << start.fitness << "\n";
            }
            else
                break;
        }
    }
    catch (const std::exception& e)
    {
        return e.what();
    }
    return std::string();
}

// дефолтный оптимизатор
gun_model::Gun::GunType Default::gun_type_impl (const Genome& genome)
{
    if (genome.params.size() != 10)
        throw std::logic_error("wrong input genome");
    for (size_t i = 0; i < genome.params.size(); ++i)
        if (std::isnan(genome.params[i].value))
            throw std::logic_error("wrong input genome");

    gun_model::Gun::GunType gt;

    gt.bubbleModel = gun_model::Gun::GunType::BubbleModel::KELLER;
    gt.massModel = gun_model::Gun::GunType::MassModel::EMPIRICAL;
    gt.decayModel = gun_model::Gun::GunType::DecayModel::VISCOSITY;
    gt.byoancyModel = gun_model::Gun::GunType::ByoancyModel::LINEAR;
    gt.A = genome.params[0].value;
    gt.M_hc = genome.params[1].value;
    gt.alpha = genome.params[2].value;
    gt.alpha_c = genome.params[3].value;
    gt.alpha_mu = genome.params[4].value;
    gt.alpha_clust = genome.params[5].value;
    gt.by = genome.params[6].value;
    gt.m_ratio = genome.params[7].value;
    gt.t_start = genome.params[8].value;
    gt.t_open = genome.params[9].value;

    return gt;
}

Genome Default::genome_impl ()
{
    Genome ret;
    const auto& lims = Default::limits();
    ret.params.emplace_back(lims.A_MIN, lims.A_MAX);
    ret.params.emplace_back(lims.M_HC_MIN, lims.M_HC_MAX);
    ret.params.emplace_back(lims.ALPHA_MIN, lims.ALPHA_MAX);
    ret.params.emplace_back(lims.ALPHA_C_MIN, lims.ALPHA_C_MAX);
    ret.params.emplace_back(lims.ALPHA_MU_MIN, lims.ALPHA_MU_MAX);
    ret.params.emplace_back(lims.ALPHA_CLUST_MIN, lims.ALPHA_CLUST_MAX, true, 2.7);
    ret.params.emplace_back(lims.BY_MIN, lims.BY_MAX);
    ret.params.emplace_back(lims.M_RATIO_MIN, lims.M_RATIO_MAX);
    ret.params.emplace_back(lims.T_START_MIN, lims.T_START_MAX, true, 0);
    ret.params.emplace_back(lims.T_OPEN_MIN, lims.T_OPEN_MAX, true, 0.015);
    ret.fitness = NAN;

    return ret;
}

void TableMaker::initialisation ()
{
    names = {"A_hvp", "M_hc_hvp", "alpha_hvp", "alpha_c_hvp", "alpha_mu_hvp", "alpha_clust_hvp", "by_hvp", "m_ratio_hvp", "t_start_hvp", "error_hvp"};
    optimized_params.resize(names.size());
    for (size_t i = 0; i < optimized_params.size(); ++i)
    {
        optimized_params[i].resize((mesh.H_END - mesh.H_START) / mesh.H_STEP + 1);
        for (size_t i_h = 0; i_h < optimized_params[i].size(); ++i_h)
        {
            optimized_params[i][i_h].resize((mesh.V_END - mesh.V_START) / mesh.V_STEP + 1);
            for (size_t i_v = 0; i_v < optimized_params[i][i_h].size(); ++i_v)
            {
                optimized_params[i][i_h][i_v].resize((mesh.P_END - mesh.P_START) / mesh.P_STEP + 1);
            }
        }
    }
}

void TableMaker::make_params ()
{
    for (size_t i_h = 0; i_h < optimized_params[0].size(); ++i_h)
    {
        int h = static_cast<int>(mesh.H_START + mesh.H_STEP * i_h);
        for (size_t i_v = 0; i_v < optimized_params[0][i_h].size(); ++i_v)
        {
            int v = static_cast<int>(mesh.V_START + mesh.V_STEP * i_v);
            for (size_t i_p = 0; i_p < optimized_params[0][i_h][i_v].size(); ++i_p)
            {
                int p = static_cast<int>(mesh.P_START + mesh.P_STEP * i_p);
                Options options;
                Mesh local_mesh = mesh;
                std::cout << "h " << h << " v " << v << " p " << p << "\n";
                local_mesh.H_START = h;
                local_mesh.H_END = h;
                local_mesh.V_START = v;
                local_mesh.V_END = v;
                local_mesh.P_START = p;
                local_mesh.P_END = p;
                Algorithm<OptimisationRegime<Default>> algorithm(options, local_mesh, dataPath, gunName, sig_optimizer::norm_TECHNICAL, true, clust);
                if (!clust)
                    algorithm.genetic();
                else
                    algorithm.geneticClust();
                for (size_t i = 0; i < names.size() - 1; ++i)
                {
                    optimized_params[i][i_h][i_v][i_p] = algorithm.population.individs.back().params[i].value;
                }
                optimized_params[names.size() - 1][i_h][i_v][i_p] = algorithm.population.individs.back().fitness;
            }
        }
    }
}

void TableMaker::assemble (nlohmann::json& jsonObject) const
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
    for (size_t i = 0; i < names.size(); ++i)
    {
        jsonObject[gunName][names[i]] = nlohmann::json(optimized_params[i]);
    }
}

} // namespace sig_optimizer