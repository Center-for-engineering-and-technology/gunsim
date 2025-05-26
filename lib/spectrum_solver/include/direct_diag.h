#pragma once

#include <fstream>
#include <vector>
#include <array>
#include <complex>
#include <cmath>
#include <string>

#include <gund_structs.h>
#include <gund_format_parser.h>
#include <gund_json_parser.h>
#include <nlohmann/json.hpp>
#include <fourier_solver.h>

namespace direct_diag {
    struct DelayedGun : gund_structs::Gun {
        double delay;
    };

    struct MathConstants {
        const double pi = std::acos(-1);
        const complex_type I = {0, 1};
    };

    struct DiagramParams {
        std::vector<DelayedGun> gunArray;
        gund_structs::PhysicalParameters physParams;
    };

    struct DiagramOptions {
        gund_structs::Reflection reflection;
        gund_structs::Filter filter;
        gund_structs::SignatureParameters sigParams;
    };

    struct OutputDiagramOptions {
        std::string outputFileName;
        double azimut; // degrees
        double maxDipAngle = 90; // degrees
        double dipIncr = 1; // degrees
        double lowerDB = 100; // dB
        double higherDB = 210; // dB 
    };

    class Diagram {
        public:
        Diagram(DiagramParams params, DiagramOptions options) : diagramParams(params), diagramOptions(options), individual_specs(params.gunArray.size()) {}
        void setDataPath(const std::string& path) { dataPath = path; }
        void prepareSpecs();
        void outputDiagram(OutputDiagramOptions &options);

        protected:
        // подготовка i-го спектра
        void computeSpectrumFromFile(size_t i);

        private:
        MathConstants math;
        DiagramParams diagramParams;
        DiagramOptions diagramOptions;
        std::string dataPath;
        std::vector<std::vector<complex_type>> individual_specs;
    };
} // namespace direct_diag
