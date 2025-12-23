#include <gund_json_parser.h>

using gund_utility::toUpperCase;

namespace gund_json_parser
{

using namespace spectrum_solver_structs;

bool parseInitJson (std::string const& inputFilename, SpectrumSolverParams& params, SpectrumSolverOptions& options)
{
    if (!std::filesystem::exists(inputFilename))
    {
        throw std::runtime_error("File could not be opened");
    }

    nlohmann::json jsonObject{nlohmann::json::parse(std::ifstream(inputFilename))};
    if (jsonObject.is_array())
    {
        jsonObject = jsonObject[0];
    }

    if (jsonObject.find("sigParams") != jsonObject.end())
    {
        jsonObject["sigParams"]["sampleInterval"].get_to(options.sigParams.sampleInterval);
        jsonObject["sigParams"]["sampleNum"].get_to(options.sigParams.sampleNum);
        options.sigParams.sampleMax = options.sigParams.sampleNum * options.sigParams.sampleInterval;
    }
    else
    {
        return false;
    }

    if (jsonObject.find("observationPoint") != jsonObject.end())
    {
        bool infiniteType = jsonObject["observationPoint"]["infiniteType"].get<bool>();
        options.obsPoint.observationType = (infiniteType) ? gund_structs::ObservationPoint::ObservationPointType::INFINITE : gund_structs::ObservationPoint::ObservationPointType::SPECIFY;
        if (!infiniteType)
        {
            jsonObject["observationPoint"]["coordinate"]["x"].get_to(options.obsPoint.x);
            jsonObject["observationPoint"]["coordinate"]["y"].get_to(options.obsPoint.y);
            jsonObject["observationPoint"]["coordinate"]["z"].get_to(options.obsPoint.z);
        }
    }
    else
    {
        return false;
    }

    if (jsonObject.find("reflection") != jsonObject.end())
    {
        jsonObject["reflection"]["coeff"].get_to(options.reflection.refCoef);
        if (jsonObject["reflection"].contains("cable_reflection"))
        {
            jsonObject["reflection"]["cable_reflection"]["cable1_depth"].get_to(options.reflection.firstCableDepth);
            jsonObject["reflection"]["cable_reflection"]["cable2_depth"].get_to(options.reflection.secondCableDepth);
        }
    }
    else
    {
        return false;
    }

    if (jsonObject.find("physicalParameters") != jsonObject.end())
    {
        jsonObject["physicalParameters"]["temperature"].get_to(params.physParams.seaTemp);
        jsonObject["physicalParameters"]["soundVelocity"].get_to(params.physParams.soundVelocity);
    }
    else
    {
        return false;
    }

    if (jsonObject.find("filter") != jsonObject.end())
    {
        std::string mode = jsonObject["filter"]["accessMode"].get<std::string>();
        toUpperCase(mode);
        // сменить запись в структуры (при наличии нескольких типов фильтра)
        if (mode == "OFF")
        {
            options.filter.bandpass.mode = gund_structs::BandpassFilter::OFF;
        }
        else if (mode == "EXTERNAL")
        {
            options.filter.bandpass.mode = gund_structs::BandpassFilter::EXTERNAL;
            jsonObject["filter"]["filename"].get_to(options.filter.bandpass.filename);
        }
        else if (mode == "INTERNAL")
        {
            // добавить параметры внутреннего фильтра
            options.filter.bandpass.mode = gund_structs::BandpassFilter::INTERNAL;
        }
    }

    if (jsonObject.find("gunArray") != jsonObject.end())
    {
        for (auto& el: jsonObject["gunArray"])
        {
            gund_structs::Gun gun;
            std::string gun_type = el["type"].get<std::string>();
            gun.type = convertGunType(gun_type);

            el["coordinate"]["x"].get_to(gun.x);
            el["coordinate"]["y"].get_to(gun.y);
            el["coordinate"]["z"].get_to(gun.z);
            el["delay"].get_to(gun.delay);
            el["volume"].get_to(gun.volume);
            el["pressure"].get_to(gun.pressure);
            el["shapeRatio"].get_to(gun.shapeRatio);
            el["temperature"].get_to(gun.temperature);

            params.gunArray.push_back(gun);
        }
    }
    else
    {
        return false;
    }

    if (jsonObject.find("directivity") != jsonObject.end())
    {
        jsonObject["directivity"]["inlineAzimutAngle"].get_to(options.dirParams.inlineAzimutAngle);
        jsonObject["directivity"]["crosslineAzimutAngle"].get_to(options.dirParams.crosslineAzimutAngle);
        jsonObject["directivity"]["maxDipAngle"].get_to(options.dirParams.maxDipAngle);
        jsonObject["directivity"]["lowerDB"].get_to(options.dirParams.lowerDB);
        jsonObject["directivity"]["higherDB"].get_to(options.dirParams.higherDB);
        jsonObject["directivity"]["dipIncr"].get_to(options.dirParams.dipIncr);
        jsonObject["directivity"]["dipIncrForSignRepr"].get_to(options.dirParams.dipIncrForSignRepr);
        options.dirParams.azimutal = jsonObject["directivity"]["azimutalMode"].get<bool>();
        if (options.dirParams.azimutal)
        {
            jsonObject["directivity"]["azimutalFreq1"].get_to(options.dirParams.azimutalFreq1);
            jsonObject["directivity"]["azimutalFreq2"].get_to(options.dirParams.azimutalFreq2);
            jsonObject["directivity"]["azimutalFreq3"].get_to(options.dirParams.azimutalFreq3);
            jsonObject["directivity"]["azimutalFreq4"].get_to(options.dirParams.azimutalFreq4);
            jsonObject["directivity"]["azimutalLowerDB"].get_to(options.dirParams.azimutalLowerDB);
            jsonObject["directivity"]["azimutalHigherDB"].get_to(options.dirParams.azimutalHigherDB);
        }
    }

    if (jsonObject.find("diffModel") != jsonObject.end())
        options.diffModel = jsonObject["diffModel"].get<bool>();

    return true;
}

gund_structs::GunType convertGunType (std::string& type)
{
    toUpperCase(type);

    gund_structs::GunType gunType = gund_structs::GunType::C1500;
    if (type == "600B")
    {
        gunType = gund_structs::GunType::B600;
    }
    else if (type == "800C")
    {
        gunType = gund_structs::GunType::C800;
    }
    else if (type == "1500C")
    {
        gunType = gund_structs::GunType::C1500;
    }
    else if (type == "1500LL")
    {
        gunType = gund_structs::GunType::LL1500;
    }
    else if (type == "1900C")
    {
        gunType = gund_structs::GunType::C1900;
    }
    else if (type == "1900D-DHS")
    {
        gunType = gund_structs::GunType::DDHS1900;
    }
    else if (type == "1900LLX")
    {
        gunType = gund_structs::GunType::LLX1900;
    }
    else if (type == "1900LLXT")
    {
        gunType = gund_structs::GunType::LLXT1900;
    }
    else if (type == "2800")
    {
        gunType = gund_structs::GunType::GUN2800;
    }
    else if (type == "2800LLX")
    {
        gunType = gund_structs::GunType::LLX2800;
    }
    else if (type == "8500APG")
    {
        gunType = gund_structs::GunType::APG8500;
    }
    else if (type == "G-GUN")
    {
        gunType = gund_structs::GunType::GGUN;
    }
    else if (type == "G-GUNII")
    {
        gunType = gund_structs::GunType::GGUNII;
    }
    else if (type == "GI-GUN")
    {
        gunType = gund_structs::GunType::GIGUN;
    }
    else if (type == "Sleeve")
    {
        gunType = gund_structs::GunType::SLEEVE;
    }
    else if (type == "Sleevell")
    {
        gunType = gund_structs::GunType::SLEEVELL;
    }
    // To Do: поменять на variant с ошибкой и в else сделать ожидание ошибки.
    return gunType;
}

std::string reconvertGunType (gund_structs::GunType const& gunType)
{
    std::string type;
    if (gunType == gund_structs::GunType::B600)
    {
        type = "600B";
    }
    else if (gunType == gund_structs::GunType::C800)
    {
        type = "800C";
    }
    else if (gunType == gund_structs::GunType::C1500)
    {
        type = "1500C";
    }
    else if (gunType == gund_structs::GunType::LL1500)
    {
        type = "1500LL";
    }
    else if (gunType == gund_structs::GunType::C1900)
    {
        type = "1900C";
    }
    else if (gunType == gund_structs::GunType::DDHS1900)
    {
        type = "1900D-DHS";
    }
    else if (gunType == gund_structs::GunType::LLX1900)
    {
        type = "1900LLX";
    }
    else if (gunType == gund_structs::GunType::LLXT1900)
    {
        type = "1900LLXT";
    }
    else if (gunType == gund_structs::GunType::GUN2800)
    {
        type = "2800";
    }
    else if (gunType == gund_structs::GunType::LLX2800)
    {
        type = "2800LLX";
    }
    else if (gunType == gund_structs::GunType::APG8500)
    {
        type = "8500APG";
    }
    else if (gunType == gund_structs::GunType::GGUN)
    {
        type = "G-GUN";
    }
    else if (gunType == gund_structs::GunType::GGUNII)
    {
        type = "G-GUNII";
    }
    else if (gunType == gund_structs::GunType::GIGUN)
    {
        type = "GI-GUN";
    }
    else if (gunType == gund_structs::GunType::SLEEVE)
    {
        type = "Sleeve";
    }
    else if (gunType == gund_structs::GunType::SLEEVELL)
    {
        type = "Sleevell";
    }

    return type;
}

void writeJson (nlohmann::json& jsonObject, SpectrumSolverParams const& params, SpectrumSolverOptions const& options)
{
    // sigParams
    jsonObject["sigParams"]["sampleInterval"] = options.sigParams.sampleInterval;
    jsonObject["sigParams"]["sampleNum"] = options.sigParams.sampleNum;
    // observationPoint
    bool infiniteType = (options.obsPoint.observationType == gund_structs::ObservationPoint::ObservationPointType::INFINITE) ? 1 : 0;
    jsonObject["diffModel"] = options.diffModel;
    jsonObject["observationPoint"]["infiniteType"] = infiniteType;
    jsonObject["observationPoint"]["coordinate"]["x"] = options.obsPoint.x;
    jsonObject["observationPoint"]["coordinate"]["y"] = options.obsPoint.y;
    jsonObject["observationPoint"]["coordinate"]["z"] = options.obsPoint.z;
    // reflection
    jsonObject["reflection"]["coeff"] = options.reflection.refCoef;
    // physicalParameters
    jsonObject["physicalParameters"]["temperature"] = params.physParams.seaTemp;
    jsonObject["physicalParameters"]["soundVelocity"] = params.physParams.soundVelocity;
    // filter
    jsonObject["filter"]["accessMode"] = options.filter.bandpass.mode;
    jsonObject["filter"]["filename"] = options.filter.bandpass.filename;
    // gunArray
    for (size_t i = 0; i < params.gunArray.size(); i++)
    {
        nlohmann::json tmp_json;

        tmp_json["type"] = reconvertGunType(params.gunArray[i].type);

        tmp_json["coordinate"]["x"] = params.gunArray[i].x;
        tmp_json["coordinate"]["y"] = params.gunArray[i].y;
        tmp_json["coordinate"]["z"] = params.gunArray[i].z;
        tmp_json["delay"] = params.gunArray[i].delay;
        tmp_json["volume"] = params.gunArray[i].volume;
        tmp_json["pressure"] = params.gunArray[i].pressure;
        tmp_json["shapeRatio"] = params.gunArray[i].shapeRatio;
        tmp_json["temperature"] = params.gunArray[i].temperature;

        jsonObject["gunArray"].push_back(tmp_json);
    }
}

using namespace gund_structs;

std::unordered_map<GunType, ErrorData> parseErrorsJson (std::string const& inputFilename)
{
    if (!std::filesystem::exists(inputFilename))
    {
        throw std::runtime_error("File could not be opened");
    }

    nlohmann::json jsonObject = nlohmann::json::parse(std::ifstream(inputFilename));
    if (jsonObject.is_array())
    {
        jsonObject = jsonObject[0];
    }

    std::unordered_map<GunType, ErrorData> hashRelErrors;

    if (jsonObject.find("gun_rel_errors") != jsonObject.end())
    {
        for (auto& el: jsonObject["gun_rel_errors"])
        {
            std::string gun_type = el["gun_type"].get<std::string>();
            GunType type = convertGunType(gun_type);

            auto peakToPeakErr = el["peak_to_peak"].get<double>();
            auto bubblePeakErr = el["bubble_peak"].get<double>();
            auto bubblePeriodErr = el["bubble_period"].get<double>();

            ErrorData errData{peakToPeakErr, bubblePeakErr, bubblePeriodErr};
            hashRelErrors[type] = errData;
        }
    }
    return hashRelErrors;
}
} // namespace gund_json_parser
