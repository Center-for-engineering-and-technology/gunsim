#include <gun_model.h>

namespace gun_model
{

// template specialization
template class RungeKuttaExplicit<6>;
template class RungeKuttaImplicit<6>;
template class GunArraySolver<RungeKuttaExplicit<6>>;
template class GunArraySolver<RungeKuttaImplicit<6>>;

namespace vector_operations {

template <size_t dim>
std::array<double, dim> operator+(const std::array<double, dim>& v1, const std::array<double, dim>& v2) {
	std::array<double, dim> ret;
	for (size_t i = 0; i < dim; ++i)
		ret[i] = v1[i] + v2[i];
	return ret;
}

template <size_t dim>
std::array<double, dim> operator-(const std::array<double, dim>& v1, const std::array<double, dim>& v2) {
	std::array<double, dim> ret;
	for (size_t i = 0; i < dim; ++i)
		ret[i] = v1[i] - v2[i];
	return ret;
}

template <size_t dim>
std::array<double, dim> operator-(const std::array<double, dim>& v) {
	std::array<double, dim> ret;
	for (size_t i = 0; i < dim; ++i)
		ret[i] = -v[i];
	return ret;
}

template <size_t dim>
std::array<double, dim> operator*(double a, const std::array<double, dim>& v) {
	std::array<double, dim> ret;
	for (size_t i = 0; i < dim; ++i)
		ret[i] = a * v[i];
	return ret;
}

template <size_t dim>
double abs(const std::array<double, dim>& v) {
	double ret = 0;
	for (size_t i = 0; i < dim; ++i)
		ret += std::abs(v[i]);
	return ret;
}

}; // namespace vector_operations

template <size_t dim>
std::array<double, dim> RungeKuttaImplicit<dim>::solve(double t, double h, std::array<double, dim> y) {
	using namespace vector_operations;

	auto k1 = RungeKutta<dim>::RHS(t + h / 4., y);
	auto k1_next = RungeKutta<dim>::RHS(t + h / 4., y + h / 4. * k1);
	int counter = 0;
	while (abs(k1 - k1_next) > 1e-8) {
		k1 = RungeKutta<dim>::RHS(t + h / 4., y + h / 4. * k1_next);
		k1_next = RungeKutta<dim>::RHS(t + h / 4., y + h / 4. * k1);
		counter++;
		if (counter > 1000)
			throw std::logic_error("too many simple iterations in runge-kutta solver");
	}
	counter = 0;
	auto k2 = RungeKutta<dim>::RHS(t + 3 * h / 4., y + h / 2. * k1_next);
	auto k2_next = RungeKutta<dim>::RHS(t + 3 * h / 4., y + h / 2. * k1_next + h / 4. * k2);
	while (abs(k2 - k2_next) > 1e-8) {
		k2 = RungeKutta<dim>::RHS(t + 3 * h / 4., y + h / 2. * k1_next + h / 4. * k2_next);
		k2_next = RungeKutta<dim>::RHS(t + 3 * h / 4., y + h / 2. * k1_next + h / 4. * k2);
		counter++;
		if (counter > 1000)
			throw std::logic_error("too many simple iterations in runge-kutta solver");
	}
	return y + h / 2. * (k1_next + k2_next);
}

template <size_t dim>
std::array<double, dim> RungeKuttaExplicit<dim>::solve(double t, double h, std::array<double, dim> y) {
	using namespace vector_operations;

	auto k1 = RungeKutta<dim>::RHS(t, y);
	auto k2 = RungeKutta<dim>::RHS(t + h / 2., y + h / 2. * k1);
	auto k3 = RungeKutta<dim>::RHS(t + h / 2., y + h / 2. * k2);
	auto k4 = RungeKutta<dim>::RHS(t + h, y + h * k3);
	return y + h / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
}

double Gun::ambientWaterPressure(double t, const PhysParams& PP) const {
	/* gidrostatic pressure */
	return 98700 + PP.rho * PP.g * std::max(GI.z - GT.by * std::max(t - GI.delay, 0.), 0.);
}

double Gun::ambientWaterPressureVelocity(double t, const PhysParams& PP) const {
	double z = GI.z - GT.by * (t - GI.delay);
	if (z > 0 && t > GI.delay)
		return -PP.rho * PP.g * GT.by;
	return 0;
}

double Gun::gunPressure(double tg, double mg) const {
	/* using perfect gas equation */
	return mg * GP.Q * tg / GI.Vol;
}

double Gun::bubblePressure(double t, double r, double tb, double mb, const PhysParams& PP) const {
	/* using perfect gas equation */
	if (t < GI.delay)
		return initialBubblePressure(PP);
	if (r < 1e-10)
		throw std::logic_error("too small bubble radius occurs: r=" + std::to_string(r));
	double ret = 3 * mb * GP.Q * tb / 4. / M_PI / std::pow(r, 3);
	if (std::isnan(ret) || std::isinf(ret))
		throw std::logic_error("nan or inf in bubble pressure");
	return ret;
}

double Gun::bubblePressureVelocity(double t, double r, double v, double tb, double mb, double mbv, double tbv) const {
	if (t < GI.delay)
		return 0;
	if (r < 1e-10)
		throw std::logic_error("too small bubble radius occurs: r=" + std::to_string(r));
	double ret = 3 * GP.Q / 4. / M_PI * (mbv * tb + mb * tbv - 3 * mb * tb * v / r) / std::pow(r, 3);
	if (std::isnan(ret) || std::isinf(ret))
		throw std::logic_error("nan or inf in bubble pressure velocity");
	return ret;
}

double Gun::initialBubbleRadius() const {
	return GT.M_r * std::pow(GI.Vol, GT.i_r_v);
}

double Gun::initialBubbleVolume() const {
	return 4 * M_PI / 3. * std::pow(initialBubbleRadius(), 3);
}

double Gun::initialWallVelocity() const {
	/* assuming initial wall speed to be zero */
	return 0;
}

double Gun::initialBubblePressure(const PhysParams& PP) const {
	/* assuming initial bubbble pressure equals to ambient pressure */
	return ambientWaterPressure(0., PP);
}

double Gun::initialBubbleMass(const PhysParams& PP) const {
	/* using perfect gas equation */
	return initialBubblePressure(PP) * initialBubbleVolume() / GP.Q / initialBubbleTemperature(PP);
}

double Gun::initialGunTemperature(const PhysParams& PP) const {
	/* assuming initial gun camera temperature coinside with water */
	return PP.T_sea;
}

double Gun::initialBubbleTemperature(const PhysParams& PP) const {
	/* assuming initial bubble temperature coinside with water */
	return PP.T_sea;
}

double Gun::initialGunMass(const PhysParams& PP) const {
	/* using perfect gas equation */
	return GI.P * GI.Vol / GP.Q / initialGunTemperature(PP);
}

double Gun::bubbleWallVelocity(double v) const {
	/* introducing a new variable */
	return v;
}

double Gun::bubbleWallAcceleration(double t, double r, double v, std::pair<double, double> ap, const PhysParams& PP, double pb, double pbv) const {
	/* using Rayleigh spherical bubble model with dumping */
	if (t < GI.delay)
		return 0;
	if (r < 1e-10)
		throw std::logic_error("too small bubble radius occurs: r=" + std::to_string(r));
	double h = (pb - ambientWaterPressure(t, PP) - ap.first) / PP.rho;
	double der_h = (pbv - ambientWaterPressureVelocity(t, PP) - ap.second) / PP.rho;
	double ret = (h - 3 / 2. * std::pow(v, 2) - (t * GT.alpha + GT.alpha_c) * v + (std::pow(v, 3) - 2 * v * h - r * der_h) / (v - PP.c)) / r;
	if (std::isnan(ret) || std::isinf(ret))
		throw std::logic_error("nan or inf in bubble wall acceleration");
	return ret;
}

double Gun::bubbleTemperatureVelocity(double t, double r, double v, double tb, double mb, double tg, const PhysParams& PP, double pb, double mbv) const {
	/* using first thermodynamics law */
	if (t < GI.delay)
		return 0;
	if (mb < 1e-10)
		throw std::logic_error("too small bubble mass occurs: m=" + std::to_string(mb));
	double ret = ((GP.cp * tg - GP.cv * tb) * mbv - 4 * M_PI * std::pow(r, 2) * (GT.M_hc * GP.kappa * (tb - PP.T_sea) + pb * v)) / GP.cv / mb;
	if (std::isnan(ret) || std::isinf(ret))
		throw std::logic_error("nan or inf in bubble temperature velocity");
	/*if (mbv > 1e-5 && ret > 5000)
		return 5000;*/
	return ret;
}

double Gun::bubbleMassVelocity(double t, double tg, double mg, double pb, double pg) const {
	/* using Bernoulli equation */
	double ret = 0;
	if (t >= GT.t_open + GI.delay || t <= GI.delay)
		return ret;
	if (pg < pb || mg < 1e-3)
		return ret;
	if (tg < 1e-10)
		throw std::logic_error("too small gun temperature occurs: tg=" + std::to_string(tg));
	if (pg > pb * std::pow((GP.gamma + 1) / 2., GP.gamma / (GP.gamma - 1.))) {
		ret = pg * GT.A * std::pow(GI.Vol, GT.i_m_v) * std::sqrt(GP.gamma / GP.Q / tg * std::pow(2 / (GP.gamma + 1.), GP.gamma * GT.i_m_p / (GP.gamma - 1.)));
		if (std::isnan(ret) || std::isinf(ret))
			throw std::logic_error("nan or inf in supersonic mass velocity");
		return ret;
	}
	if (pb < 1e-10 || pg < 1e-10)
		throw std::logic_error("negative pressures occur: pb=" + std::to_string(pb) + ", pg=" + std::to_string(pg));
	ret = pg * GT.A * std::pow(GI.Vol, GT.i_m_v) * std::sqrt((GP.gamma / GP.Q / tg) * (2 / (GP.gamma - 1.)) * std::pow(pb / pg, GT.i_m_p) * (std::pow(pg / pb, (GP.gamma - 1) / GP.gamma) - 1));
	if (std::isnan(ret) || std::isinf(ret))
		throw std::logic_error("nan or inf in undersonic mass velocity");
	return ret;
}

double Gun::gunTemperatureVelocity(double t, double tg, double mg, double mgv) const {
	/* using first thermodynamics law */
	if (t < GI.delay)
		return 0;
	if (mg < 1e-9)
		throw std::logic_error("too small gun mass occurs: mg=" + std::to_string(mg));
	double ret = (GP.Q / GP.cv) * (tg / mg) * mgv;
	if (std::isnan(ret) || std::isinf(ret))
		throw std::logic_error("nan or inf in gun temperature velocity");
	return ret;
}

double Gun::gunMassVelocity(double mbv) const {
	/* opposite to bubble mass velocity */
	return -mbv;
}

size_t GunArray::size() const {
	return G.size();
}

Gun& GunArray::operator[](size_t g) {
	return G[g];
}

const Gun& GunArray::operator[](size_t g) const {
	return G[g];
}

void GunArray::emplace_back(double P_, double Vol_, double delay_, double x_, double y_, double z_) {
	G.emplace_back(P_, Vol_, delay_, x_, y_, z_);
}

void GunArray::emplace_back(const Gun::GunType& GT_, const Gun::GunPhysics& GP_, const Gun::GunInputs& GI_) {
	G.emplace_back(GT_, GP_, GI_);
}

template <typename MRK>
void GunArraySolver<MRK>::initialisation() {

	R.resize(GA.size());
	V.resize(GA.size());
	TB.resize(GA.size());
	MB.resize(GA.size());
	TG.resize(GA.size());
	MG.resize(GA.size());

	S_r.resize(GA.size());
	for (size_t g = 0; g < S_r.size(); ++g) {
		S_r[g].resize(MP.N);
		for (size_t i = 0; i < S_r[g].size(); ++i)
			S_r[g][i] = NAN;
	}
	S_der_f.resize(GA.size());
	for (size_t g = 0; g < S_der_f.size(); ++g) {
		S_der_f[g].resize(MP.N);
		for (size_t i = 0; i < S_der_f[g].size(); ++i)
			S_der_f[g][i] = NAN;
	}

	cache.resize(GA.size());

	for (size_t g = 0; g < GA.size(); ++g) {
		std::function<std::array<double, 6>(double, std::array<double, 6>)> RHS;
		RHS = [g, this](double t, std::array<double, 6> r_v_tb_mb_tg_mg)->std::array<double, 6> {
			std::array<double, 6> ret;
			double r = r_v_tb_mb_tg_mg[0]; // bubble raius
			double v = r_v_tb_mb_tg_mg[1]; // wall velocity
			double tb = r_v_tb_mb_tg_mg[2]; // bubble temperature
			double mb = r_v_tb_mb_tg_mg[3]; // bubble mass
			double tg = r_v_tb_mb_tg_mg[4]; // gun temperature
			double mg = r_v_tb_mb_tg_mg[5]; // gun mass
				
			ret[0] = GA[g].bubbleWallVelocity(v); 
			double pb = GA[g].bubblePressure(t, r, tb, mb, PP); // bubble pressure
			double pg = GA[g].gunPressure(tg, mg); // gun pressure
			double mbv = ret[3] = GA[g].bubbleMassVelocity(t, tg, mg, pb, pg); // bubble mass velocity
			double mgv = ret[5] = GA[g].gunMassVelocity(mbv); // gun mass velocity
			double tbv = ret[2] = GA[g].bubbleTemperatureVelocity(t, r, v, tb, mb, tg, PP, pb, mbv); // bubble temperature velocity
			double pbv = GA[g].bubblePressureVelocity(t, r, v, tb, mb, mbv, tbv); // bubble pressure velocity
			ret[1] = GA[g].bubbleWallAcceleration(t, r, v, acousticPressureComputation(g, t), PP, pb, pbv); 
			ret[4] = GA[g].gunTemperatureVelocity(t, tg, mg, mgv);

			return ret;
		};
		RK.emplace_back(RHS);
	}

	dist.resize(GA.size());
	for (size_t g = 0; g < GA.size(); ++g) {
		dist[g].resize(2 * GA.size());
		if (GA[g].GI.z < 0)
			throw std::logic_error("sources must lie below water");
		for (size_t h = 0; h < GA.size(); ++h) {
			dist[g][h] = std::sqrt(std::pow(GA[g].GI.x - GA[h].GI.x, 2) + std::pow(GA[g].GI.y - GA[h].GI.y, 2) + std::pow(GA[g].GI.z - GA[h].GI.z, 2));
			if (g != h && dist[g][h] < 0.099)
				throw std::logic_error("too close configuration of sources");
			if (g != h && 2 * MP.dt > dist[g][h] / PP.c)
				throw std::logic_error("time step dt is too big for the source configuration");
			dist[g][GA.size() + h] = std::sqrt(std::pow(GA[g].GI.x - GA[h].GI.x, 2) + std::pow(GA[g].GI.y - GA[h].GI.y, 2) + std::pow(GA[g].GI.z + GA[h].GI.z, 2));
			if (dist[g][GA.size() + h] < 0.099)
				throw std::logic_error("too close configuration of ghosts");
			if (2 * MP.dt > dist[g][GA.size() + h] / PP.c)
				throw std::logic_error("time step dt is too big for the source configuration");
		}
	}
}

template <typename MRK>
GunArraySolver<MRK>::GunArraySolver(const PhysParams& PP_, const MethodParams& MP_, const GunArray& GA_, Cached cache_name_) : cache_name(cache_name_), PP(PP_), MP(MP_), GA(GA_) {
	initialisation();
}

template <typename MRK>
std::string GunArraySolver<MRK>::solve() {
	omp_set_num_threads(MP.numtreads);
	std::vector<std::exception_ptr> exceptions(MP.numtreads, nullptr);
	try {
		setInitialData();
	} catch (const std::logic_error& e) {
		return e.what();
	}
	for (size_t i = 1; i < MP.N; ++i)
		try { 
			oneTimeStep(i, exceptions);
		} catch (const std::logic_error& e) {
			return e.what();
		}
	return "Solving done well";
}

template <typename MRK>
void GunArraySolver<MRK>::getSignal(double x, double y, double z) const {
	/* gives the signal [bar] in itput point */
	if (z <= 0)
		throw std::logic_error("the reciever must lay below water");
	std::vector<double> length(2 * GA.size());

	for (size_t g = 0; g < GA.size(); ++g) {
		length[g] = std::sqrt(std::pow(x - GA[g].GI.x, 2) + std::pow(y - GA[g].GI.y, 2) + std::pow(z - GA[g].GI.z, 2));
		length[GA.size() + g] = std::sqrt(std::pow(x - GA[g].GI.x, 2) + std::pow(y - GA[g].GI.y, 2) + std::pow(z + GA[g].GI.z, 2));
		if (length[g] < 1e-1 || length[GA.size() + g] < 1e-1)
			throw std::logic_error("the reciever is too close to a source");
	}

	for (size_t i = 0; i < MP.N; ++i) {
		double signal = 0;
		for (size_t g = 0; g < length.size(); ++g) {
			int j = static_cast<int>(std::floor(i - length[g] / PP.c / MP.dt));
			double bar = i - length[g] / PP.c / MP.dt - j;
			double left_sig = 0;
			if (j >= 0 && j < static_cast<int>(S_der_f[g % GA.size()].size()))
				left_sig = PP.rho * (-S_der_f[g % GA.size()][j] / length[g]);
			double right_sig = 0;
			if (j + 1 >= 0 && j + 1 < static_cast<int>(S_der_f[g % GA.size()].size()))
				right_sig = PP.rho * (-S_der_f[g % GA.size()][j + 1] / length[g]);
			if (g < GA.size())
				signal += bar * right_sig + (1 - bar) * left_sig;
			else
				signal += PP.ref * (bar * right_sig + (1 - bar) * left_sig);
		}
		if (i % MP.repr_step == 0)
			std::cout << signal / 1e5 << ",  "; // [bar]
	}
	std::cout << "\n";
}

template <typename MRK>
void GunArraySolver<MRK>::getSignatureInfinity() const {
	/* gives the signature [bar * m] in infinite point on z axis without damping effect */
	double z_max = 0;
	for (size_t g = 0; g < GA.size(); ++g)
		if (GA[g].GI.z > z_max)
			z_max = GA[g].GI.z;

	std::vector<double> z_delta(2 * GA.size());
	for (size_t g = 0; g < GA.size(); ++g) {
		z_delta[g] = z_max - GA[g].GI.z;
		z_delta[g + GA.size()] = z_max + GA[g].GI.z;
	}

	for (size_t i = 0; i < MP.N; ++i) {
		double signature = 0;
		for (size_t g = 0; g < z_delta.size(); ++g) {
			int j = static_cast<int>(std::floor(i - z_delta[g] / PP.c / MP.dt));
			double bar = i - z_delta[g] / PP.c / MP.dt - j;
			double left_sig = 0;
			if (j >= 0 && j < static_cast<int>(S_der_f[g % GA.size()].size()))
				left_sig = -PP.rho * S_der_f[g % GA.size()][j];
			double right_sig = 0;
			if (j + 1 >= 0 && j + 1 < static_cast<int>(S_der_f[g % GA.size()].size()))
				right_sig = -PP.rho * S_der_f[g % GA.size()][j + 1];
			if (g < GA.size())
				signature += bar * right_sig + (1 - bar) * left_sig;
			else
				signature += PP.ref * (bar * right_sig + (1 - bar) * left_sig);
		}
		if (i % MP.repr_step == 0)
			std::cout << signature / 1e5 << ",  "; // [bar * m]
	}
	std::cout << "\n";
}

template <typename MRK>
GunModelResult GunArraySolver<MRK>::getResult() const {
	GunModelResult result;
	result.signatures.resize(GA.size());
	for (size_t g = 0; g < GA.size(); ++g) {
		result.signatures[g].resize(MP.repr_N);
		for (size_t i = 0; i < static_cast<size_t>(MP.repr_N); ++i) {
			int j_start = static_cast<int>(std::floor(GA[g].GI.delay / MP.dt));
			result.signatures[g][i] = 0;
			if (j_start + i * MP.repr_step >= 0 && j_start + i * MP.repr_step < S_der_f[g].size())
				result.signatures[g][i] = -PP.rho * S_der_f[g][j_start + i * MP.repr_step] / 1e5; // Pa * m to bar * m
		}
	}
	return result;
}

template <typename MRK>
void GunArraySolver<MRK>::getCache(size_t g) const {
	for (size_t i = 0; i < cache[g].size(); ++i)
		if (i % MP.repr_step == 0)
			std::cout << cache[g][i] << ", ";
	std::cout << "\n";
}

template <typename MRK>
void GunArraySolver<MRK>::signatureComputation(size_t g, size_t i) {
		
	if (std::isnan(S_r[g][i]) && std::isnan(S_der_f[g][i])) {
		S_r[g][i] = R[g];
		double h = (GA[g].bubblePressure(i * MP.dt, R[g], TB[g], MB[g], PP) - GA[g].ambientWaterPressure(i * MP.dt, PP) - acousticPressureComputation(g, i * MP.dt).first) / PP.rho;
		S_der_f[g][i] = -R[g] * (h + 0.5 * std::pow(V[g], 2));
	}
	
    if (std::isnan(S_r[g][i]) || std::isinf(S_r[g][i]) || std::isnan(S_der_f[g][i]) || std::isinf(S_der_f[g][i]))
		throw std::logic_error("nan or inf signature on gun " + std::to_string(g) + " and step " + std::to_string(i) + " with parameters: r=" + std::to_string(R[g]) + ", v=" + std::to_string(V[g]) + ", tb=" + std::to_string(TB[g]) + ", mb=" + std::to_string(MB[g]) + ", tg=" + std::to_string(TG[g]) + ", mg=" + std::to_string(MG[g]));
}

template <typename MRK>
std::pair<double, double> GunArraySolver<MRK>::acousticPressureComputation(size_t g, double t) const {
	/* returns acoustic pressure and acoustic pressure velocity of gun g at moment t */
	std::pair<double, double> ret = { 0, 0 };
	std::pair<double, double> ret_covered = { 0, 0 };
	int covered_count = 0;
	for (size_t h = 0; h < 2 * GA.size(); h++)
		if (h != g && dist[g][h] < MP.INTERACTION_DIST) {
			bool covered = false;
			double distance = dist[g][h];
			int j = static_cast<int>(std::floor((t - distance / PP.c) / MP.dt));
			double bar = (t - distance / PP.c) / MP.dt - j;

			double left_h_rad = 0, left_g_rad = 0;
			if (j >= 0 && j < static_cast<int>(S_r[h % GA.size()].size())) {
				left_h_rad = S_r[h % GA.size()][j];
				left_g_rad = S_r[g][j];
			}
			double right_h_rad = 0, right_g_rad = 0;
			if (j + 1 >= 0 && j + 1 < static_cast<int>(S_r[h % GA.size()].size())) {
				right_h_rad = S_r[h % GA.size()][j + 1];
				right_g_rad = S_r[g][j + 1];
			}
			double h_rad = bar * right_h_rad + (1 - bar) * left_h_rad;
			double g_rad = bar * right_g_rad + (1 - bar) * left_g_rad;
			if (g_rad + h_rad > MP.NEAR_INTERACTION_RATIO * distance) {
				covered_count += 1;
				covered = true;
				distance = std::max(distance, g_rad + h_rad);
			}

			double left_sig = 0;
			if (j >= 0 && j < static_cast<int>(S_der_f[h % GA.size()].size()))
				left_sig = PP.rho * (-S_der_f[h % GA.size()][j] / distance);
			double right_sig = 0;
			if (j + 1 >= 0 && j + 1 < static_cast<int>(S_der_f[h % GA.size()].size()))
				right_sig = PP.rho * (-S_der_f[h % GA.size()][j + 1] / distance);
			if (h >= GA.size()) {
				left_sig *= PP.ref;
				right_sig *= PP.ref;
			}

			if (!covered) {
				ret.first += bar * right_sig + (1 - bar) * left_sig;
				ret.second += (right_sig - left_sig) / MP.dt;
			}
			else {
				ret_covered.first += bar * right_sig + (1 - bar) * left_sig;
				ret_covered.second += (right_sig - left_sig) / MP.dt;
			}
		}

	if (covered_count == 0)
		return ret;
	else {
		ret_covered.first += ret.first;
		ret_covered.first /= static_cast<double>(covered_count);
		ret_covered.second += ret.second;
		ret_covered.second /= static_cast<double>(covered_count);
	}
	return ret_covered;
}

template <typename MRK>
void GunArraySolver<MRK>::setInitialData() {

	for (size_t g = 0; g < GA.size(); ++g) {
		R[g] = GA[g].initialBubbleRadius();
		V[g] = GA[g].initialWallVelocity();
		TB[g] = GA[g].initialBubbleTemperature(PP);
		MB[g] = GA[g].initialBubbleMass(PP);
		TG[g] = GA[g].initialGunTemperature(PP);
		MG[g] = GA[g].initialGunMass(PP);
			
		signatureComputation(g, 0);
	}
}

template <typename MRK>
void GunArraySolver<MRK>::oneTimeStep(size_t i, std::vector<std::exception_ptr> &exceptions) {
    /* moving from (i - 1) * MP.dt to i * MP.dt time step */

	#pragma omp parallel for
	for (int g = 0; g < static_cast<int>(GA.size()); ++g) {
		if (cache_name != Cached::NONE) 
			switch (cache_name) {
				case Cached::NONE: break;
				case Cached::BUBBLE_RADIUS: cache[g].push_back(R[g]); break;
				case Cached::WALL_VELOCITY: cache[g].push_back(V[g]); break;
				case Cached::BUBBLE_TEMPERATURE: cache[g].push_back(TB[g]); break;
				case Cached::BUBBLE_MASS: cache[g].push_back(MB[g]); break;
				case Cached::BUBBLE_PRESSURE: cache[g].push_back(GA[g].bubblePressure((i - 1) * MP.dt, R[g], TB[g], MB[g], PP)); break;
				case Cached::GUN_TEMPERATURE: cache[g].push_back(TG[g]); break;
				case Cached::GUN_MASS: cache[g].push_back(MG[g]); break;
				case Cached::GUN_PRESSURE: cache[g].push_back(GA[g].gunPressure(TG[g], MG[g])); break;
			}
		try {
			if ((i - 1) * MP.dt <= GA[g].GT.t_open + GA[g].GI.delay && i * MP.dt >= GA[g].GI.delay) {
				std::array<double, 6> r_v_tb_mb_tg_mg = RK[g].solve((i - 1) * MP.dt, MP.dt / static_cast<double>(MP.REFINEMENT), { R[g], V[g], TB[g], MB[g], TG[g], MG[g] });
				for (size_t j = 1; j < MP.REFINEMENT; ++j)
					r_v_tb_mb_tg_mg = RK[g].solve((i - 1) * MP.dt + j * MP.dt / static_cast<double>(MP.REFINEMENT), MP.dt / static_cast<double>(MP.REFINEMENT), r_v_tb_mb_tg_mg);
				
				R[g] = r_v_tb_mb_tg_mg[0];
				V[g] = r_v_tb_mb_tg_mg[1];
				TB[g] = r_v_tb_mb_tg_mg[2];
				MB[g] = r_v_tb_mb_tg_mg[3];
				TG[g] = r_v_tb_mb_tg_mg[4];
				MG[g] = r_v_tb_mb_tg_mg[5];
			} else {
				std::array<double, 6> r_v_tb_mb_tg_mg = RK[g].solve((i - 1) * MP.dt, MP.dt, { R[g], V[g], TB[g], MB[g], TG[g], MG[g] });
				
				R[g] = r_v_tb_mb_tg_mg[0];
				V[g] = r_v_tb_mb_tg_mg[1];
				TB[g] = r_v_tb_mb_tg_mg[2];
				MB[g] = r_v_tb_mb_tg_mg[3];
				TG[g] = r_v_tb_mb_tg_mg[4];
				MG[g] = r_v_tb_mb_tg_mg[5];
			}
			
			signatureComputation(g, i);
		} catch (...) {
			exceptions[omp_get_thread_num()] = std::current_exception();
		}
	}

	/* перехват исключения из параллельной секции */
	for (const auto& ex : exceptions) {
		if (ex) {
			try {
				std::rethrow_exception(ex);
			} catch (const std::logic_error &e) {
				throw;
			}
		}
	}
}

}; // namespace gun_model
