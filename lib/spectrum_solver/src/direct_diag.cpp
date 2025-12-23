#include <direct_diag.h>

namespace direct_diag
{
void Diagram::prepareSpecs ()
{
    for (size_t i = 0; i < diagramParams.gunArray.size(); i++)
    {
        computeSpectrumFromFile(i);
    }
}

void Diagram::prepareSpecs (const gun_model::GunModelResult& data)
{
    for (size_t i = 0; i < diagramParams.gunArray.size(); i++)
    {
        computeSpectrumFromGunModel(i, data);
    }
}

void Diagram::computeSpectrumFromFile (size_t i)
{

    auto& gun = diagramParams.gunArray[i];
    std::string signalFilename = gund_json_parser::reconvertGunType(gun.type) + "_";
    signalFilename += (std::to_string(( int )gun.z) + "m_");
    signalFilename += ("V" + std::to_string(( int )gun.volume) + "_");
    signalFilename += ("P" + std::to_string(( int )gun.pressure));
    signalFilename += ".sig";
    const auto file = dataPath + signalFilename;

    gund_format_parser::GundalfOutputParser parser(file);
    const std::vector<double>& parsed_signal = parser.getData();
    std::vector<double> signal(static_cast<size_t>(diagramOptions.sigParams.sampleNum));
    signal = parsed_signal;

    fourier_solver::ConcreteFourierSolver solver(diagramOptions.sigParams);
    solver.solve(signal);

    if (diagramOptions.filter.bandpass.mode == gund_structs::BandpassFilter::EXTERNAL)
    {
        // применение внешнего bandpass фильтра
        const std::string fltFile = diagramOptions.filter.bandpass.filename;
        if (!fltFile.empty())
        {
            gund_format_parser::GundalfOutputParser fltParser(dataPath + fltFile);
            diagramOptions.filter.bandpass.filterData = fltParser.getData();
            diagramOptions.filter.bandpass.filterSigParams = fltParser.getSigParams();
        }
        fourier_solver::FilterFourierSolver filteredSolver(solver, diagramOptions.filter.bandpass.filterSigParams, diagramOptions.filter.bandpass.filterData);
        filteredSolver.applyFilter();

        const std::vector<complex_type>& computed_spectrum = filteredSolver.getSpectrum();
        individual_specs[i].resize(computed_spectrum.size());
        individual_specs[i] = computed_spectrum;
    }
    else
    {
        const std::vector<complex_type>& computed_spectrum = solver.getSpectrum();
        individual_specs[i].resize(computed_spectrum.size());
        individual_specs[i] = computed_spectrum;
    }
}

void Diagram::computeSpectrumFromGunModel (size_t i, const gun_model::GunModelResult& data)
{

    std::vector<double> signal = data.signatures[i];

    fourier_solver::ConcreteFourierSolver solver(diagramOptions.sigParams);
    solver.solve(signal);

    if (diagramOptions.filter.bandpass.mode == gund_structs::BandpassFilter::EXTERNAL)
    {
        // применение внешнего bandpass фильтра
        const std::string fltFile = diagramOptions.filter.bandpass.filename;
        if (!fltFile.empty())
        {
            gund_format_parser::GundalfOutputParser fltParser(dataPath + fltFile);
            diagramOptions.filter.bandpass.filterData = fltParser.getData();
            diagramOptions.filter.bandpass.filterSigParams = fltParser.getSigParams();
        }
        fourier_solver::FilterFourierSolver filteredSolver(solver, diagramOptions.filter.bandpass.filterSigParams, diagramOptions.filter.bandpass.filterData);
        filteredSolver.applyFilter();

        const std::vector<complex_type>& computed_spectrum = filteredSolver.getSpectrum();
        individual_specs[i].resize(computed_spectrum.size());
        individual_specs[i] = computed_spectrum;
    }
    else
    {
        const std::vector<complex_type>& computed_spectrum = solver.getSpectrum();
        individual_specs[i].resize(computed_spectrum.size());
        individual_specs[i] = computed_spectrum;
    }
}

void Diagram::outputDiagram (OutputDiagramOptions& options)
{

    size_t angleMeshSize = static_cast<size_t>(2 * options.maxDipAngle / options.dipIncr);
    double angleMeshStep = math.pi / 180. * options.dipIncr;                                                                                                                             // radians
    double freqMeshStep = 1. / diagramOptions.sigParams.sampleInterval / static_cast<double>(gund_utility::findNextPowerOfTwo(static_cast<size_t>(diagramOptions.sigParams.sampleNum))); // Hz
    size_t freqMeshSize = static_cast<size_t>(250. / freqMeshStep);
    std::vector<double> angles((angleMeshSize + 1) * (freqMeshSize + 1)),
        frequences((angleMeshSize + 1) * (freqMeshSize + 1)),
        values((angleMeshSize + 1) * (freqMeshSize + 1));
    double T = diagramOptions.sigParams.sampleInterval * diagramOptions.sigParams.sampleNum;

    double azimutAngle = -options.azimut * math.pi / 180.; // radians

    for (size_t k = 0; k < freqMeshSize + 1; ++k)
    {
        auto pol_ang = -options.maxDipAngle * math.pi / 180.;
        for (size_t q = 0; q < angleMeshSize + 1; ++q)
        {
            complex_type J_kq = {0., 0.};
            for (size_t g = 0; g < diagramParams.gunArray.size(); ++g)
            {
                double s = diagramParams.gunArray[g].x * std::cos(azimutAngle) * std::sin(pol_ang) + diagramParams.gunArray[g].y * std::sin(azimutAngle) * std::sin(pol_ang) + diagramParams.gunArray[g].z * std::cos(pol_ang);
                J_kq += individual_specs[g][k] * std::exp(-2. * math.pi * math.I * (diagramParams.gunArray[g].delay - s / diagramParams.physParams.soundVelocity) * ( double )k / T);
                if (diagramOptions.reflection.refCoef != 0.)
                {
                    double s_ref = diagramParams.gunArray[g].x * std::cos(azimutAngle) * std::sin(pol_ang) + diagramParams.gunArray[g].y * std::sin(azimutAngle) * std::sin(pol_ang) - diagramParams.gunArray[g].z * std::cos(pol_ang);
                    J_kq += diagramOptions.reflection.refCoef * individual_specs[g][k] * std::exp(-2. * math.pi * math.I * (diagramParams.gunArray[g].delay - s_ref / diagramParams.physParams.soundVelocity) * ( double )k / T);
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

void Diagram::outputSignatureDiagram (OutputSignatureDiagramOptions& options)
{
    std::vector<double> dipAngles{-2. * options.dipIncrForSignRepr / 180. * math.pi, -options.dipIncrForSignRepr / 180. * math.pi, 0, options.dipIncrForSignRepr / 180. * math.pi, 2 * options.dipIncrForSignRepr / 180. * math.pi};
    std::vector<std::vector<complex_type>> spectrums;
    std::vector<std::vector<double>> signatures, amplitudes;
    spectrums.resize(5);
    signatures.resize(5);
    amplitudes.resize(5);

    size_t specSize = 0;
    double freqMeshStep = 1. / diagramOptions.sigParams.sampleInterval / static_cast<double>(gund_utility::findNextPowerOfTwo(static_cast<size_t>(diagramOptions.sigParams.sampleNum)));
    if (individual_specs.size() > 0)
        specSize = individual_specs[0].size();
    double azimutAngle = -options.azimut * math.pi / 180.; // radians
    double T = diagramOptions.sigParams.sampleInterval * diagramOptions.sigParams.sampleNum;

    fourier_solver::ConcreteFourierSolver solver(diagramOptions.sigParams);
    double max_depth = 0;
    for (size_t g = 0; g < diagramParams.gunArray.size(); ++g)
        if (diagramParams.gunArray[g].z > max_depth)
            max_depth = diagramParams.gunArray[g].z;

    for (size_t i = 0; i < 5; ++i)
    {
        spectrums[i].resize(specSize);
        for (size_t k = 0; k < specSize; ++k)
        {
            spectrums[i][k] = {0., 0.};
            for (size_t g = 0; g < diagramParams.gunArray.size(); ++g)
            {
                double s = diagramParams.gunArray[g].x * std::cos(azimutAngle) * std::sin(dipAngles[i]) + diagramParams.gunArray[g].y * std::sin(azimutAngle) * std::sin(dipAngles[i]) + diagramParams.gunArray[g].z * std::cos(dipAngles[i]);
                spectrums[i][k] += individual_specs[g][k] * std::exp(-2. * math.pi * math.I * (diagramParams.gunArray[g].delay - s / diagramParams.physParams.soundVelocity) * ( double )k / T);
                if (diagramOptions.reflection.refCoef != 0.)
                {
                    double s_ref = diagramParams.gunArray[g].x * std::cos(azimutAngle) * std::sin(dipAngles[i]) + diagramParams.gunArray[g].y * std::sin(azimutAngle) * std::sin(dipAngles[i]) - diagramParams.gunArray[g].z * std::cos(dipAngles[i]);
                    spectrums[i][k] += diagramOptions.reflection.refCoef * individual_specs[g][k] * std::exp(-2. * math.pi * math.I * (diagramParams.gunArray[g].delay - s_ref / diagramParams.physParams.soundVelocity) * ( double )k / T);
                }
            }
            spectrums[i][k] *= std::exp(-2. * math.pi * math.I * max_depth / diagramParams.physParams.soundVelocity * ( double )k / T);
        }
        std::vector<complex_type> tmp_spec(gund_utility::findNextPowerOfTwo(specSize));
        for (size_t k = 0; k < tmp_spec.size(); ++k)
            if (k < tmp_spec.size() / 2)
                tmp_spec[k] = spectrums[i][k];
            else
                tmp_spec[k] = std::conj(spectrums[i][tmp_spec.size() - k]);
        signatures[i] = solver.inverseFastFourierTransform(tmp_spec);
        signatures[i].erase(signatures[i].begin() + static_cast<size_t>(diagramOptions.sigParams.sampleNum), signatures[i].end());

        amplitudes[i].resize(specSize);
        std::transform(spectrums[i].begin(), spectrums[i].end(), amplitudes[i].begin(), [] (const complex_type& z)
                       {
                           return 10. * std::log10(std::pow(std::abs(z), 2) + 1e-20);
                       });
    }

    nlohmann::json j_spec, j_sign;
    j_spec["mesh size"] = specSize;
    j_spec["mesh step"] = freqMeshStep;
    for (size_t i = 0; i < 5; ++i)
    {
        j_spec["dip_" + std::to_string(i + 1)]["dip"] = dipAngles[i] / math.pi * 180;
        j_spec["dip_" + std::to_string(i + 1)]["values"] = nlohmann::json(amplitudes[i]);
    }
    j_sign["mesh size"] = static_cast<size_t>(diagramOptions.sigParams.sampleNum);
    j_sign["mesh step"] = diagramOptions.sigParams.sampleInterval;
    for (size_t i = 0; i < 5; ++i)
    {
        j_sign["dip_" + std::to_string(i + 1)]["dip"] = dipAngles[i] / math.pi * 180;
        j_sign["dip_" + std::to_string(i + 1)]["values"] = nlohmann::json(signatures[i]);
    }

    std::ofstream output_spec_json, output_sign_json;
    output_spec_json.open(dataPath + "spec_repr_" + options.outputFileName);
    output_sign_json.open(dataPath + "sign_repr_" + options.outputFileName);
    output_spec_json << std::setw(4) << std::setprecision(17) << j_spec;
    output_sign_json << std::setw(4) << std::setprecision(17) << j_sign;
    output_spec_json.close();
    output_sign_json.close();
}

void Diagram::outputAzimutalDiagram (OutputAzimutalDiagramOptions& options)
{
    size_t angleMeshSize = 64;
    double angleMeshStep = math.pi / static_cast<double>(angleMeshSize);
    double freqMeshStep = 1. / diagramOptions.sigParams.sampleInterval / static_cast<double>(gund_utility::findNextPowerOfTwo(static_cast<size_t>(diagramOptions.sigParams.sampleNum))); // Hz
    std::vector<double> graph_x(angleMeshSize * angleMeshSize),
        graph_y(angleMeshSize * angleMeshSize),
        values(angleMeshSize * angleMeshSize);
    double T = diagramOptions.sigParams.sampleInterval * diagramOptions.sigParams.sampleNum;
    int k_f = static_cast<int>(std::floor(options.frequency / freqMeshStep));
    double bar = options.frequency / freqMeshStep - k_f;

    for (size_t k_x = 0; k_x < angleMeshSize; ++k_x)
    {
        double ang_x = -math.pi / 2. + k_x * angleMeshStep;
        for (size_t k_y = 0; k_y < angleMeshSize; ++k_y)
        {
            double ang_y = -math.pi / 2. + k_y * angleMeshStep;
            complex_type J = {0., 0.};
            double theta = std::sqrt(std::pow(ang_x, 2) + std::pow(ang_y, 2));
            if (theta < math.pi / 2.)
            {
                double cos_phi = 0.;
                double sin_phi = 1.;
                if (theta > 1e-6)
                {
                    cos_phi = ang_y / theta;
                    sin_phi = ang_x / theta;
                }
                for (size_t g = 0; g < diagramParams.gunArray.size(); ++g)
                {
                    double s = diagramParams.gunArray[g].x * cos_phi * std::sin(theta) + diagramParams.gunArray[g].y * sin_phi * std::sin(theta) + diagramParams.gunArray[g].z * std::cos(theta);
                    if (k_f >= 0 && k_f < static_cast<int>(individual_specs[g].size()))
                        J += (1 - bar) * individual_specs[g][k_f] * std::exp(-2. * math.pi * math.I * (diagramParams.gunArray[g].delay - s / diagramParams.physParams.soundVelocity) * static_cast<double>(k_f) / T);
                    if (k_f + 1 >= 0 && k_f + 1 < static_cast<int>(individual_specs[g].size()))
                        J += bar * individual_specs[g][k_f + 1] * std::exp(-2. * math.pi * math.I * (diagramParams.gunArray[g].delay - s / diagramParams.physParams.soundVelocity) * static_cast<double>(k_f + 1) / T);
                    if (diagramOptions.reflection.refCoef != 0.)
                    {
                        double s_ref = diagramParams.gunArray[g].x * cos_phi * std::sin(theta) + diagramParams.gunArray[g].y * sin_phi * std::sin(theta) - diagramParams.gunArray[g].z * std::cos(theta);
                        if (k_f >= 0 && k_f < static_cast<int>(individual_specs[g].size()))
                            J += (1 - bar) * diagramOptions.reflection.refCoef * individual_specs[g][k_f] * std::exp(-2. * math.pi * math.I * (diagramParams.gunArray[g].delay - s_ref / diagramParams.physParams.soundVelocity) * static_cast<double>(k_f) / T);
                        if (k_f + 1 >= 0 && k_f + 1 < static_cast<int>(individual_specs[g].size()))
                            J += bar * diagramOptions.reflection.refCoef * individual_specs[g][k_f + 1] * std::exp(-2. * math.pi * math.I * (diagramParams.gunArray[g].delay - s_ref / diagramParams.physParams.soundVelocity) * static_cast<double>(k_f + 1) / T);
                    }
                }
            }
            graph_x[k_x * angleMeshSize + k_y] = ang_x / math.pi * 180;
            graph_y[k_x * angleMeshSize + k_y] = ang_y / math.pi * 180;
            values[k_x * angleMeshSize + k_y] = 10. * std::log10(std::pow(std::abs(J), 2) + 1e-20);
        }
    }

    nlohmann::json j;
    j["diagram description"]["x-count"] = angleMeshSize;
    j["diagram description"]["x-step"] = angleMeshStep / math.pi * 180;
    j["diagram description"]["x-start"] = -90;
    j["diagram description"]["y-count"] = angleMeshSize;
    j["diagram description"]["y-step"] = angleMeshStep / math.pi * 180;
    j["diagram description"]["y-start"] = -90;
    j["diagram description"]["values diapazone"] = {options.lowerDB, options.higherDB};
    j["x-axis"] = nlohmann::json(graph_x);
    j["y-axis"] = nlohmann::json(graph_y);
    j["values"] = nlohmann::json(values);

    std::ofstream output_json;
    output_json.open(dataPath + options.outputFileName);
    output_json << std::setw(4) << j;
    output_json.close();
}
} // namespace direct_diag
