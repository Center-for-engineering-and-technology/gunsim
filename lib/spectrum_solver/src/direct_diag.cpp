#include <direct_diag.h>

namespace direct_diag
{
    void Diagram::prepareSpecs() {
        for (size_t i = 0; i < diagramParams.gunArray.size(); i++) {
            computeSpectrumFromFile(i);
        }
    }

    void Diagram::computeSpectrumFromFile(size_t i) {
        
        auto &gun = diagramParams.gunArray[i];
        std::string signalFilename = gund_json_parser::reconvertGunType(gun.type) + "_";
        signalFilename += (std::to_string((int)gun.z) + "m_");
        signalFilename += ("V" + std::to_string((int)gun.volume) + "_");
        signalFilename += ("P" + std::to_string((int)gun.pressure));
        signalFilename += ".sig";
        const auto file = dataPath + signalFilename;

        gund_format_parser::GundalfOutputParser parser(file);
        const std::vector<double>& parsed_signal = parser.getData();
        std::vector<double> signal(diagramOptions.sigParams.sampleNum);
        signal = parsed_signal;

        fourier_solver::ConcreteFourierSolver solver(diagramOptions.sigParams);
        solver.solve(signal);

        const std::vector<complex_type>& computed_spectrum = solver.getSpectrum();
        individual_specs[i].resize(computed_spectrum.size());
        individual_specs[i] = computed_spectrum;
    }

    void Diagram::outputDiagram(OutputDiagramOptions &options) {
        
        size_t angleMeshSize = static_cast<size_t>(2 * options.maxDipAngle / options.dipIncr); 
        double angleMeshStep = math.pi / 180. * options.dipIncr; // radians
        double freqMeshStep = 1. / diagramOptions.sigParams.sampleInterval / static_cast<double>(gund_utility::findNextPowerOfTwo(static_cast<size_t>(diagramOptions.sigParams.sampleNum))); // Hz
        size_t freqMeshSize = static_cast<size_t>(250. / freqMeshStep);
        std::vector<double> angles((angleMeshSize + 1) * (freqMeshSize + 1)),
                            frequences((angleMeshSize + 1) * (freqMeshSize + 1)),
                            values((angleMeshSize + 1) * (freqMeshSize + 1));
        double T = diagramOptions.sigParams.sampleInterval * diagramOptions.sigParams.sampleNum;

        double azimutAngle = -options.azimut * math.pi / 180.; // radians

        for (size_t k = 0; k < freqMeshSize + 1; ++k) {
            auto pol_ang = -options.maxDipAngle * math.pi / 180.;
            for (size_t q = 0; q < angleMeshSize + 1; ++q) {
                complex_type J_kq = {0., 0.};
                for (size_t g = 0; g < diagramParams.gunArray.size(); ++g) {
                    double s =  diagramParams.gunArray[g].x * std::cos(azimutAngle) * std::sin(pol_ang) +
                                diagramParams.gunArray[g].y * std::sin(azimutAngle) * std::sin(pol_ang) +
                                diagramParams.gunArray[g].z * std::cos(pol_ang);
                    J_kq += individual_specs[g][k] * std::exp(-2. * math.pi * math.I * (diagramParams.gunArray[g].delay - s / diagramParams.physParams.soundVelocity) * (double)k / T);
                    if (diagramOptions.reflection.refCoef != 0.) {
                        double s_ref =  diagramParams.gunArray[g].x * std::cos(azimutAngle) * std::sin(pol_ang) +
                                        diagramParams.gunArray[g].y * std::sin(azimutAngle) * std::sin(pol_ang) -
                                        diagramParams.gunArray[g].z * std::cos(pol_ang);
                        J_kq += diagramOptions.reflection.refCoef * individual_specs[g][k] * std::exp(-2. * math.pi * math.I * (diagramParams.gunArray[g].delay - s_ref / diagramParams.physParams.soundVelocity) * (double)k / T);
                    }
                }
                angles[k * (angleMeshSize + 1) + q] = pol_ang / math.pi * 180; // degrees
                frequences[k * (angleMeshSize + 1) + q] = k * freqMeshStep;
                values[k * (angleMeshSize + 1) + q] = std::pow(std::abs(J_kq), 2);
                values[k * (angleMeshSize + 1) + q] = 10. * std::log10(values[k * (angleMeshSize + 1) + q] + 1e-20);
                pol_ang += angleMeshStep;
            }
        }

        nlohmann::json j;
        j["diagram description"]["x-count"] = angleMeshSize + 1;
        j["diagram description"]["x-step"] = options.dipIncr;
        j["diagram description"]["x-start"] = -options.maxDipAngle;
        j["diagram description"]["y-count"] = freqMeshSize + 1;
        j["diagram description"]["y-step"] = freqMeshStep;
        j["diagram description"]["y-start"] = 0;
        j["diagram description"]["values diapazone"] = {options.lowerDB, options.higherDB};
        j["x-axis"] = nlohmann::json(angles);
        j["y-axis"] = nlohmann::json(frequences);
        j["values"] = nlohmann::json(values);

        std::ofstream output_json;
        output_json.open(dataPath + options.outputFileName);
        output_json << std::setw(4) << j;
        output_json.close();
    }

} // namespace direct_diag
