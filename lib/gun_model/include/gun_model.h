#pragma once

#define _USE_MATH_DEFINES

#include <array>
#include <cmath>
#include <functional>
#include <gun_parameters_table.h>
#include <gund_structs.h>
#include <iostream>
#include <omp.h>
#include <spectrum_solver_structs.h>
#include <string>
#include <vector>

namespace gun_model
{

template<size_t dim>
class RungeKutta
{
public:

    explicit RungeKutta (const std::function<std::array<double, dim>(double, std::array<double, dim>)>& RHS_)
        : RHS(RHS_)
    {
    }
    virtual std::array<double, dim> solve(double t, double h, std::array<double, dim> y) = 0;

protected:

    std::function<std::array<double, dim>(double, std::array<double, dim>)> RHS;
};

template<size_t dim>
class RungeKuttaImplicit final : public RungeKutta<dim>
{
public:

    using RungeKutta<dim>::RungeKutta;
    std::array<double, dim> solve(double t, double h, std::array<double, dim> y) override;
};

template<size_t dim>
class RungeKuttaExplicit final : public RungeKutta<dim>
{
public:

    using RungeKutta<dim>::RungeKutta;
    std::array<double, dim> solve(double t, double h, std::array<double, dim> y) override;
};

struct PhysParams
{
    double c = 1496;    // [m/s] salt water speed of sound
    double rho = 1020.; // [kg/m^3] salt water density
    double g = 9.806;   // [m/s^2] gravitation acceleration
    double ref = -1.;   // [] reflection coefficient
    double T_sea = 283; // [K] water temperature
    double mu = 1.3e-3; // [kg/m/s] dynamical viscosity

    PhysParams() = default;
    PhysParams(const spectrum_solver_structs::SpectrumSolverParams& params, const spectrum_solver_structs::SpectrumSolverOptions& options);
};

struct MethodParams
{
    size_t N = 4000; // mesh size
    double sampleMax = 0.5;
    double dt = sampleMax / static_cast<double>(N); // [s] time step
    int repr_N = 4000;
    int repr_step = static_cast<int>(N) / repr_N; // outputing step
    int numtreads = 12;                           // the number of processors

    double INTERACTION_DIST = 20; // [m] distance between sources above that interaction neglected
    size_t REFINEMENT = 40;
    double NEAR_INTERACTION_RATIO = 0.6;

    MethodParams() = default;
    MethodParams(const spectrum_solver_structs::SpectrumSolverParams& params, const spectrum_solver_structs::SpectrumSolverOptions& options);
};

class Gun
{
public:

    struct GunType
    {                          // leaning parameters
        double A = 0.005;      // [m^2] cross-sectional gun area
        double t_open = 0.008; // [s] time that port is open for
        double t_start = 0.002;
        double M_hc = 1.2;        // magnification of heat conduction
        double alpha = 1.6;       // model damping parameter
        double alpha_c = 0.3;     // model damping parameter
        double alpha_clust = 0.3; // model damping parameter
        double alpha_mu = 0.1;    // model damping parameter
        double by = 5;            // [m/s] buoyancy speed
        double m_ratio = 1;

        enum struct MassModel
        {
            ADIABATIC,
            EMPIRICAL
        } massModel = MassModel::EMPIRICAL;
        enum struct BubbleModel
        {
            RELEIGH,
            KELLER
        } bubbleModel = BubbleModel::KELLER;
        enum struct DecayModel
        {
            LINEAR,
            LINEAR_R,
            VISCOSITY
        } decayModel = DecayModel::VISCOSITY;
        enum struct ByoancyModel
        {
            LINEAR,
            HERRING
        } byoancyModel = ByoancyModel::LINEAR;

        GunType ()
        {
        }
        GunType(const gun_parameters_table::GunMap& gun_map, const gund_structs::Gun& gun, MassModel m_model, BubbleModel b_model, DecayModel d_model, ByoancyModel by_model);
    } GT;

    struct GunPhysics
    {
        double kappa = 4000;    // [J/(m^2*s*K)] gas/water heat conduction coefficient
        double cv = 717.258;    // [J/(kg*K)] gas specific heat capacity at constant volume
        double cp = 1006.313;   // [J/(kg*K)] gas specific heat capacity at constant pressure
        double gamma = cp / cv; // [] gas adiabatic index
        double Q = cp - cv;     //[J/(kg*K)] gas specific constant
    } GP;

    struct GunInputs
    {
        double P;       // [Pa] gun pressure
        double Vol;     // [m^3] gun volume
        double delay;   // [s] time-delay of firing
        double x, y, z; // [m] gun position

        GunInputs ()
        {
        }
        GunInputs (double P_, double Vol_, double delay_, double x_, double y_, double z_)
            : P(P_), Vol(Vol_), delay(delay_), x(x_), y(y_), z(z_)
        {
        }
        explicit GunInputs(const gund_structs::Gun& gun);
    } GI;

    Gun ()
    {
    }
    Gun (const GunType& GT_, const GunPhysics& GP_, const GunInputs& GI_)
        : GT(GT_), GP(GP_), GI(GI_)
    {
    }
    Gun (const gun_parameters_table::GunMap& gun_map, const gund_structs::Gun& gun)
        : GT(gun_map, gun, GunType::MassModel::EMPIRICAL, GunType::BubbleModel::KELLER, GunType::DecayModel::VISCOSITY, GunType::ByoancyModel::LINEAR), GP(), GI(gun)
    {
    }

    double ambientWaterPressure(double t, const PhysParams& PP) const;
    double ambientWaterPressureVelocity(double t, const PhysParams& PP) const;
    double gunPressure(double tg, double mg) const;
    double bubblePressure(double t, double r, double tb, double mb, const PhysParams& PP) const;
    double bubblePressureVelocity(double t, double r, double v, double tb, double mb, double mbv, double tbv) const;

    double initialBubbleRadius() const;
    double initialBubbleVolume() const;
    double initialWallVelocity() const;
    double initialBubbleTemperature(const PhysParams& PP) const;
    double initialBubblePressure(const PhysParams& PP) const;
    double initialBubbleMass(const PhysParams& PP) const;
    double initialGunTemperature(const PhysParams& PP) const;
    double initialGunMass(const PhysParams& PP) const;

    double bubbleWallVelocity(double v) const;
    double bubbleDecay(double t, double r, double v, double re) const;
    double bubbleDecayVelocityWithWallAccelerationMultiplier(double t, double r, double v, double u_other) const;
    double bubbleDecayVelocityOtherTerm(double t, double r, double v, double u_other, double u_other_velo) const;
    double bubbleWallAcceleration(double t, double r, double v, std::array<double, 4> apv, const PhysParams& PP, double pb, double pbv) const;
    double bubbleTemperatureVelocity(double t, double r, double v, double tb, double mb, const PhysParams& PP, double pb, double mbv) const;
    double bubbleMassVelocity(double t, double tb, double mg, const PhysParams& PP, double pb, double pg) const;
    double gunTemperatureVelocity(double t, double tg, double mg, double mgv) const;
    double gunMassVelocity(double mbv) const;
};

class GunArray
{
public:

    GunArray() = default;
    GunArray(const gun_parameters_table::GunMap& gun_map, const spectrum_solver_structs::SpectrumSolverParams& params, const spectrum_solver_structs::SpectrumSolverOptions& options);

    size_t size() const;
    void resize(size_t size_val);
    Gun& operator[](size_t g);
    const Gun& operator[](size_t g) const;
    void emplace_back(const Gun::GunType& GT_, const Gun::GunPhysics& GP_, const Gun::GunInputs& GI_);
    void emplace_back(const gun_parameters_table::GunMap& gun_map, const gund_structs::Gun& gun);

private:

    std::vector<Gun> G;
};

struct GunModelResult
{
    std::vector<std::vector<double>> signatures;
};

template<typename MRK = RungeKuttaExplicit<6>>
class GunArraySolver
{
public:

    enum Cached
    {
        NONE,
        BUBBLE_RADIUS,
        WALL_VELOCITY,
        BUBBLE_TEMPERATURE,
        BUBBLE_MASS,
        BUBBLE_PRESSURE,
        GUN_TEMPERATURE,
        GUN_MASS,
        GUN_PRESSURE
    } cache_name;
    enum InteractionModel
    {
        DEFAULT_INTERACTION,
        CLOSE_INTERACTION
    } interaction_model;

    GunArraySolver(const PhysParams& PP_, const MethodParams& MP_, const GunArray& GA_, Cached cache_name_, InteractionModel interaction_model_);
    GunArraySolver(const gun_parameters_table::GunMap& gun_map, const spectrum_solver_structs::SpectrumSolverParams& params, const spectrum_solver_structs::SpectrumSolverOptions& options, InteractionModel interaction_model_);
    void initialisation();

    std::string solve();
    std::vector<double> getSignal(double x, double y, double z, bool print) const;
    std::vector<double> getSignatureInfinity(bool print) const;
    std::vector<double> getCache(size_t g, bool print) const;

    GunModelResult getResult() const;

private:

    PhysParams PP;
    MethodParams MP;
    GunArray GA;

    std::vector<double> R;  // bubble radiuses
    std::vector<double> V;  // bubble wall velocities
    std::vector<double> TB; // bubble temperatures
    std::vector<double> MB; // bubble masses
    std::vector<double> TG; // gun temperatures
    std::vector<double> MG; // gun masses

    std::vector<std::vector<double>> cache;   // some computed physical values
    std::vector<MRK> RK;                      // Runge-Kutta solver structs
    std::vector<std::vector<double>> dist;    // [m] distances between sourses and ghosts
    std::vector<std::vector<double>> S_r;     // [m] bubble radiuses
    std::vector<std::vector<double>> S_f;     // [m^3 * s] far field signature component
    std::vector<std::vector<double>> S_der_f; // [m^3 * s^2] far field signature compunent

    void signatureComputation(size_t g, size_t i);
    std::array<double, 4> acousticPressureAndVelocityComputation(size_t g, double t) const;
    void setInitialData();
    void oneTimeStep(size_t i, std::vector<std::exception_ptr>& exceptions);
};

void CheckUpGun(const gun_parameters_table::GunMap& table, gund_structs::GunType gun_type, double min_p, double max_p, double min_v, double max_v, double min_z, double max_z, size_t mesh_size, std::vector<std::array<double, 3>>& errors);

}; // namespace gun_model
