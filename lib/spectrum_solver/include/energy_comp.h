#pragma once

#include <fstream>
#include <vector>
#include <array>
#include <complex>
#include <cmath>
#include <string>
#include <iostream>

#include <gund_structs.h>
#include <gund_format_parser.h>
#include <gund_json_parser.h>
#include <nlohmann/json.hpp>
#include <fourier_solver.h>

namespace energy_comp {
    struct DelayedGun : gund_structs::Gun {
        double delay;
    };

    struct ExtendedPhysicalParameters : gund_structs::PhysicalParameters {
        double density = 1020.;
    };

    struct MathConstants {
        const double pi = std::acos(-1);
        const complex_type I = {0, 1};
    };

    struct EnergyParams {
        std::vector<DelayedGun> gunArray;
        ExtendedPhysicalParameters physParams;
    };

    struct EnergyOptions {
        gund_structs::Reflection reflection;
        gund_structs::Filter filter;
        gund_structs::SignatureParameters sigParams;
    };

    struct EnergyResult {
        std::vector<double> energy;
        struct { double x, y, z; } energyCenter;
        double totalAcousticEnergy;
        double totalPotentialEnergy;
        double acousticEnergyEffectiveness;
    };

    class Energy {
        public:
        Energy(EnergyParams params, EnergyOptions options) :    energyParams(params), energyOptions(options) { }
        void setDataPath(const std::string& path) { dataPath = path; }
        void prepareSignalsAndLengths();
        void solve();
        EnergyResult getResult() const { return energyResult; }
        
        protected:
        // подготовка сигнатур
        void computeSignalFromFile(size_t i);
        std::pair<int, double> lower_tick(double time, double dt) const;
        double interpolation(double time, size_t i) const;

        private:
        MathConstants math;
        EnergyParams energyParams;
        EnergyOptions energyOptions;
        std::string dataPath;

        std::vector<std::vector<double>> individual_signals;
        std::vector<std::vector<double>> lengths;
        EnergyResult energyResult;
    };
} // namespace energy_comp
