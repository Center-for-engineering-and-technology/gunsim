#pragma once

#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <functional>
#include <cmath>
#include <omp.h>

namespace gun_model {

template <size_t dim>
class RungeKutta {
public:
	explicit RungeKutta(const std::function<std::array<double, dim>(double, std::array<double, dim>)>& RHS_) : RHS(RHS_) { }
	virtual std::array<double, dim> solve(double t, double h, std::array<double, dim> y) = 0;

protected:
	std::function<std::array<double, dim>(double, std::array<double, dim>)> RHS;
};

template <size_t dim>
class RungeKuttaImplicit final : public RungeKutta<dim> {
public:
	using RungeKutta<dim>::RungeKutta;
	std::array<double, dim> solve(double t, double h, std::array<double, dim> y) override;
};

template <size_t dim>
class RungeKuttaExplicit final : public RungeKutta<dim> {
public:
	using RungeKutta<dim>::RungeKutta;
	std::array<double, dim> solve(double t, double h, std::array<double, dim> y) override;
};

struct PhysParams {
	double c = 1496; // [m/s] salt water speed of sound
	double rho = 1020.; // [kg/m^3] salt water density
	double g = 9.806; // [m/s^2] gravitation acceleration
	double ref = -1.; // [] reflection coefficient
	double T_sea = 282; // [K] water temperature

	PhysParams() = default;
};

struct MethodParams {
	size_t N = 2000; // mesh size
	double sampleMax = 0.5;
	double dt = sampleMax / static_cast<double>(N); // [s] time step
	int repr_N = 1000;
	int repr_step = static_cast<int>(N) / repr_N; // outputing step
    int numtreads = 12; // the number of processors
	
	double INTERACTION_DIST = 12; // [m] distance between sources above that interaction neglected
	size_t REFINEMENT = 20;
	double NEAR_INTERACTION_RATIO = 0.6;

	MethodParams() = default;
};

class Gun {
public:
	struct GunType { // leaning parameters
		double A = 0.005; // [m^2] cross-sectional gun area
		double t_open = 0.008; // [s] time that port is open for
		double M_hc = 1.2; // magnification of heat conduction
		double alpha = 1.6; // [m/s^2] model damping parameter
		double alpha_c = 0.3; // [m/s] model damping parameter
		double by = 5; // [m/s] buoyancy speed
		double M_r = 0.7; // magnification of initial bubble radius
		double i_r_v = 0.3333; // volume index of initial radius
		double i_m_v = 0; // volume index of mass flow
		double i_m_p = 1; // pressure index of mass flow

		GunType() { }
	} GT;

	struct GunPhysics {
		double kappa = 4000; // [J/(m^2*s*K)] gas/water heat conduction coefficient
		double cv = 717.258; // [J/(kg*K)] gas specific heat capacity at constant volume
		double cp = 1006.313; // [J/(kg*K)] gas specific heat capacity at constant pressure
		double gamma = cp / cv; // [] gas adiabatic index
		double Q = cp - cv; //[J/(kg*K)] gas specific constant
	} GP;

	struct GunInputs {
		double P; // [Pa] gun pressure
		double Vol; // [m^3] gun volume
		double delay; // [s] time-delay of firing
		double x, y, z; // [m] gun position

		GunInputs(double P_, double Vol_, double delay_, double x_, double y_, double z_) :
			 P(P_), Vol(Vol_), delay(delay_), x(x_), y(y_), z(z_) { }
	} GI;

	Gun(double P_, double Vol_, double delay_, double x_, double y_, double z_) :
		GT(), GP(), GI(P_, Vol_, delay_, x_, y_, z_) { }
	Gun(const GunType& GT_, const GunPhysics& GP_, const GunInputs& GI_) :
		GT(GT_), GP(GP_), GI(GI_) { }

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
	double bubbleWallAcceleration(double t, double r, double v, std::pair<double, double> ap, const PhysParams& PP, double pb, double pbv) const;
	double bubbleTemperatureVelocity(double t, double r, double v, double tb, double mb, double tg, const PhysParams& PP, double pb, double mbv) const;
	double bubbleMassVelocity(double t, double tg, double mg, double pb, double pg) const;
	double gunTemperatureVelocity(double t, double tg, double mg, double mgv) const;
	double gunMassVelocity(double mbv) const;
};

class GunArray {
public:
	GunArray() = default;

	size_t size() const;
	Gun& operator[](size_t g);
	const Gun& operator[](size_t g) const;
	void emplace_back(double P_, double Vol_, double delay_, double x_, double y_, double z_);
	void emplace_back(const Gun::GunType& GT_, const Gun::GunPhysics& GP_, const Gun::GunInputs& GI_);

private:
	std::vector<Gun> G;
};

struct GunModelResult {
	std::vector<std::vector<double>> signatures;
};

template<typename MRK = RungeKuttaExplicit<6>>
class GunArraySolver {
public:
	enum Cached { NONE, BUBBLE_RADIUS, WALL_VELOCITY, BUBBLE_TEMPERATURE, BUBBLE_MASS, BUBBLE_PRESSURE, GUN_TEMPERATURE, GUN_MASS, GUN_PRESSURE } cache_name;
	
	GunArraySolver(const PhysParams& PP_, const MethodParams& MP_, const GunArray& GA_, Cached cache_name_);

	void initialisation();

    std::string solve();
	void getSignal(double x, double y, double z) const;
	void getSignatureInfinity() const;
	void getCache(size_t g) const;
	
	GunModelResult getResult() const;

private:
	PhysParams PP;
	MethodParams MP;
	GunArray GA;

	std::vector<double> R; // bubble radiuses
	std::vector<double> V; // bubble wall velocities
	std::vector<double> TB; // bubble temperatures
	std::vector<double> MB; // bubble masses
	std::vector<double> TG; // gun temperatures
	std::vector<double> MG; // gun masses

	std::vector<std::vector<double>> cache; // some computed physical values
    std::vector<MRK> RK; // Runge-Kutta solver structs
	std::vector<std::vector<double>> dist; // [m] distances between sourses and ghosts
	std::vector<std::vector<double>> S_r; // [m] bubble radiuses
	std::vector<std::vector<double>> S_der_f; // [Pa * m] far field signature 

	void signatureComputation(size_t g, size_t i);
	std::pair<double, double> acousticPressureComputation(size_t g, double t) const;
	void setInitialData();
	void oneTimeStep(size_t i, std::vector<std::exception_ptr> &exceptions);

};

}; // namespace gun_model
