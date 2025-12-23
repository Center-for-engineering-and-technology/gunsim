#include <energy_comp.h>

namespace energy_comp
{
void Energy::prepareSignalsAndLengths ()
{
    individual_signals.resize(energyParams.gunArray.size());
    lengths.resize(energyParams.gunArray.size());
    for (size_t g = 0; g < energyParams.gunArray.size(); g++)
    {
        computeSignalFromFile(g);
    }
    for (size_t g = 0; g < energyParams.gunArray.size(); g++)
    {
        lengths[g].resize(2 * energyParams.gunArray.size());
        for (size_t ag = 0; ag < 2 * energyParams.gunArray.size(); ag++)
        {
            lengths[g][ag] = std::sqrt(
                std::pow(energyParams.gunArray[g].x - energyParams.gunArray[ag % energyParams.gunArray.size()].x, 2) + std::pow(energyParams.gunArray[g].y - energyParams.gunArray[ag % energyParams.gunArray.size()].y, 2) + std::pow(energyParams.gunArray[g].z - std::pow(-1, ag / energyParams.gunArray.size()) * energyParams.gunArray[ag % energyParams.gunArray.size()].z, 2)
            );
        }
    }
}

void Energy::prepareSignalsAndLengths (const gun_model::GunModelResult& data)
{
    individual_signals.resize(energyParams.gunArray.size());
    lengths.resize(energyParams.gunArray.size());
    for (size_t g = 0; g < energyParams.gunArray.size(); g++)
    {
        computeSignalFromGunModel(g, data);
    }
    for (size_t g = 0; g < energyParams.gunArray.size(); g++)
    {
        lengths[g].resize(2 * energyParams.gunArray.size());
        for (size_t ag = 0; ag < 2 * energyParams.gunArray.size(); ag++)
        {
            lengths[g][ag] = std::sqrt(
                std::pow(energyParams.gunArray[g].x - energyParams.gunArray[ag % energyParams.gunArray.size()].x, 2) + std::pow(energyParams.gunArray[g].y - energyParams.gunArray[ag % energyParams.gunArray.size()].y, 2) + std::pow(energyParams.gunArray[g].z - std::pow(-1, ag / energyParams.gunArray.size()) * energyParams.gunArray[ag % energyParams.gunArray.size()].z, 2)
            );
        }
    }
}

void Energy::computeSignalFromFile (size_t g)
{
    auto& gun = energyParams.gunArray[g];
    std::string signalFilename = gund_json_parser::reconvertGunType(gun.type) + "_";
    signalFilename += (std::to_string(( int )gun.z) + "m_");
    signalFilename += ("V" + std::to_string(( int )gun.volume) + "_");
    signalFilename += ("P" + std::to_string(( int )gun.pressure));
    signalFilename += ".sig";
    const auto file = dataPath + signalFilename;

    gund_format_parser::GundalfOutputParser parser(file);
    const std::vector<double>& parsed_signal = parser.getData();
    std::vector<double> signal(static_cast<size_t>(energyOptions.sigParams.sampleNum));
    signal = parsed_signal; // geting signal in Bars

    if (energyOptions.filter.bandpass.mode == gund_structs::BandpassFilter::EXTERNAL)
    {
        // применение внешнего bandpass фильтра
        fourier_solver::ConcreteFourierSolver solver(energyOptions.sigParams);
        solver.solve(signal);
        const std::string fltFile = energyOptions.filter.bandpass.filename;
        if (!fltFile.empty())
        {
            gund_format_parser::GundalfOutputParser fltParser(dataPath + fltFile);
            energyOptions.filter.bandpass.filterData = fltParser.getData();
            energyOptions.filter.bandpass.filterSigParams = fltParser.getSigParams();
        }
        fourier_solver::FilterFourierSolver filteredSolver(solver, energyOptions.filter.bandpass.filterSigParams, energyOptions.filter.bandpass.filterData);
        filteredSolver.applyFilter();

        const std::vector<double>& filtered_signal = filteredSolver.getSignal();
        individual_signals[g].resize(filtered_signal.size());
        individual_signals[g] = filtered_signal;
    }
    else
    {
        individual_signals[g].resize(signal.size());
        individual_signals[g] = signal;
    }
}

void Energy::computeSignalFromGunModel (size_t g, const gun_model::GunModelResult& data)
{
    std::vector<double> signal = data.signatures[g]; // geting signal in Bars

    if (energyOptions.filter.bandpass.mode == gund_structs::BandpassFilter::EXTERNAL)
    {
        // применение внешнего bandpass фильтра
        fourier_solver::ConcreteFourierSolver solver(energyOptions.sigParams);
        solver.solve(signal);
        const std::string fltFile = energyOptions.filter.bandpass.filename;
        if (!fltFile.empty())
        {
            gund_format_parser::GundalfOutputParser fltParser(dataPath + fltFile);
            energyOptions.filter.bandpass.filterData = fltParser.getData();
            energyOptions.filter.bandpass.filterSigParams = fltParser.getSigParams();
        }
        fourier_solver::FilterFourierSolver filteredSolver(solver, energyOptions.filter.bandpass.filterSigParams, energyOptions.filter.bandpass.filterData);
        filteredSolver.applyFilter();

        const std::vector<double>& filtered_signal = filteredSolver.getSignal();
        individual_signals[g].resize(filtered_signal.size());
        individual_signals[g] = filtered_signal;
    }
    else
    {
        individual_signals[g].resize(signal.size());
        individual_signals[g] = signal;
    }
}

void Energy::solve ()
{
    energyResult.energy.resize(energyParams.gunArray.size(), 0.);
    energyResult.totalAcousticEnergy = 0.;
    energyResult.totalPotentialEnergy = 0.;
    energyResult.energyCenter.x = 0.;
    energyResult.energyCenter.y = 0.;
    energyResult.energyCenter.z = 0.;
    for (size_t g = 0; g < energyParams.gunArray.size(); ++g)
    {
        energyResult.totalPotentialEnergy += energyParams.gunArray[g].pressure * energyParams.gunArray[g].volume * 0.0000163871 * 6894.75728; // inner gun's energy in Joules
        double antider_i = 0.;
        for (size_t i = 0; i < 2 * energyOptions.sigParams.sampleNum; ++i)
        {
            double sig_i = interpolation(i * energyOptions.sigParams.sampleInterval - energyParams.gunArray[g].delay, g);
            antider_i += sig_i * energyOptions.sigParams.sampleInterval;
            energyResult.energy[g] += std::pow(sig_i, 2) * energyOptions.sigParams.sampleInterval / energyParams.physParams.soundVelocity;
            for (size_t ag = 0; ag < 2 * energyParams.gunArray.size(); ++ag)
                if (ag != g)
                {
                    double asig_i = interpolation(i * energyOptions.sigParams.sampleInterval - energyParams.gunArray[ag % energyParams.gunArray.size()].delay - lengths[g][ag] / energyParams.physParams.soundVelocity, ag % energyParams.gunArray.size()) / lengths[g][ag];
                    if (ag / energyParams.gunArray.size() == 1)
                        asig_i *= energyOptions.reflection.refCoef;
                    energyResult.energy[g] += asig_i * antider_i * energyOptions.sigParams.sampleInterval;
                }
        }
        energyResult.energy[g] *= 4. * math.pi / energyParams.physParams.density * 1e10; // geting energy in Joules [1 Bar = 10^5 Pa]
        energyResult.totalAcousticEnergy += energyResult.energy[g];
    }
    energyResult.acousticEnergyEffectiveness = energyResult.totalAcousticEnergy / energyResult.totalPotentialEnergy * 100.; // %-value of efficiency

    for (size_t g = 0; g < energyParams.gunArray.size(); ++g)
    {
        energyResult.energyCenter.x += energyResult.energy[g] * energyParams.gunArray[g].x / energyResult.totalAcousticEnergy;
        energyResult.energyCenter.y += energyResult.energy[g] * energyParams.gunArray[g].y / energyResult.totalAcousticEnergy;
        energyResult.energyCenter.z += energyResult.energy[g] * energyParams.gunArray[g].z / energyResult.totalAcousticEnergy;
    }
}

std::pair<int, double> Energy::lower_tick (double time, double dt) const
{
    int i = ( int )std::floor(time / dt);
    double bar = time / dt - i;
    return {i, bar};
}

double Energy::interpolation (double time, size_t g) const
{
    auto lt = lower_tick(time, energyOptions.sigParams.sampleInterval);
    int i = lt.first;
    double bar = lt.second;
    double sig_left;
    double sig_right;
    if (i < 0 || i >= ( int )individual_signals[g].size())
        sig_left = 0.;
    else
        sig_left = individual_signals[g][i];
    if (i + 1 < 0 || i + 1 >= ( int )individual_signals[g].size())
        sig_right = 0.;
    else
        sig_right = individual_signals[g][i + 1];
    return (1. - bar) * sig_left + bar * sig_right;
}

} // namespace energy_comp
