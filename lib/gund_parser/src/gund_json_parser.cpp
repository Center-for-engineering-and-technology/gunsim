#include <gund_json_parser.h>

using gund_utility::toUpperCase;

namespace gund_json_parser {

using namespace spectrum_solver_structs;

    bool parseInitJson(std::string const& inputFilename, SpectrumSolverParams& params, SpectrumSolverOptions& options) {
        std::string jsonString;
        std::ifstream file;

        file.open(inputFilename);
        if (!file.is_open()) {
            throw std::runtime_error("File could not be opened");
        }
        // To Do: проверить тип файла

        std::getline(file, jsonString, '\0');

        file.close();

        nlohmann::json jsonObject{ nlohmann::json::parse(jsonString) };

        if (jsonObject.find("sigParams") != jsonObject.end()) {
            jsonObject["sigParams"]["sampleInterval"].get_to(options.sigParams.sampleInterval);
            jsonObject["sigParams"]["sampleNum"].get_to(options.sigParams.sampleNum);
            options.sigParams.sampleMax = options.sigParams.sampleNum * options.sigParams.sampleInterval;
        }
        else {
            return false;
        }

        if (jsonObject.find("observationPoint") != jsonObject.end()) {
            bool infiniteType = jsonObject["observationPoint"]["infiniteType"].get<bool>();
            options.obsPoint.observationType = (infiniteType) ? 
                gund_structs::ObservationPoint::ObservationPointType::INFINITE : gund_structs::ObservationPoint::ObservationPointType::SPECIFY;
            if (!infiniteType) {
                jsonObject["observationPoint"]["coordinate"]["x"].get_to(options.obsPoint.x);
                jsonObject["observationPoint"]["coordinate"]["y"].get_to(options.obsPoint.y);
                jsonObject["observationPoint"]["coordinate"]["z"].get_to(options.obsPoint.z);
            }
        }
        else {
            return false;
        }

        if (jsonObject.find("reflection") != jsonObject.end()) {
            jsonObject["reflection"]["coeff"].get_to(options.reflection.refCoef);
        }
        else {
            return false;
        }

        if (jsonObject.find("physicalParameters") != jsonObject.end()) {
            jsonObject["physicalParameters"]["temperature"].get_to(params.physParams.seaTemp);
            jsonObject["physicalParameters"]["soundVelocity"].get_to(params.physParams.soundVelocity);
        }
        else {
            return false;
        }

        if (jsonObject.find("filter") != jsonObject.end()) {
            std::string mode = jsonObject["filter"]["accessMode"].get<std::string>();
            toUpperCase(mode);
            // сменить запись в структуры (при наличии нескольких типов фильтра)
            if (mode == "OFF") { options.filter.bandpass.mode = gund_structs::BandpassFilter::OFF; }
            else if (mode == "EXTERNAL") {
                options.filter.bandpass.mode = gund_structs::BandpassFilter::EXTERNAL;
                jsonObject["filter"]["filename"].get_to(options.filter.bandpass.filename);
            }
            else if (mode == "INTERNAL") {
                // добавить параметры внутреннего фильтра
                options.filter.bandpass.mode = gund_structs::BandpassFilter::INTERNAL;
            }
        }

        if (jsonObject.find("gunArray") != jsonObject.end()) {
            for (auto& el : jsonObject["gunArray"]) {
                gund_structs::Gun gun;
                std::string gun_type = el["type"].get<std::string>();
                gun.type = convertGunType(gun_type);

                el["coordinate"]["x"].get_to(gun.x);
                el["coordinate"]["y"].get_to(gun.y);
                el["coordinate"]["z"].get_to(gun.z);

                el["volume"].get_to(gun.volume);
                el["pressure"].get_to(gun.pressure);
                el["shapeRatio"].get_to(gun.shapeRatio);
                el["temperature"].get_to(gun.temperature);

                params.gunArray.push_back(gun);
            }
        }
        else {
            return false;
        }

        return true;
    }

    gund_structs::GunType convertGunType(std::string& type) {
        toUpperCase(type);

        gund_structs::GunType gunType;
        if (type == "1500C") { gunType = gund_structs::GunType::C1500; }
        else if (type == "1900C") { gunType = gund_structs::GunType::C1900; }

        return gunType;
    }

    std::string reconvertGunType(gund_structs::GunType const& gunType) {
        std::string type;
        if (gunType == gund_structs::GunType::C1500) { type = "1500C"; }
        else if (gunType == gund_structs::GunType::C1900) { type = "1900C"; }

        return type;
    }

    void writeJson(nlohmann::json& jsonObject, SpectrumSolverParams const& params, SpectrumSolverOptions const& options) {
        // sigParams
        jsonObject["sigParams"]["sampleInterval"] = options.sigParams.sampleInterval;
        jsonObject["sigParams"]["sampleNum"] = options.sigParams.sampleNum;
        // observationPoint
        bool infiniteType = (options.obsPoint.observationType == gund_structs::ObservationPoint::ObservationPointType::INFINITE) ? 1 : 0;
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
        for (size_t i = 0; i < params.gunArray.size(); i++) {
            nlohmann::json tmp_json;

            tmp_json["type"] = reconvertGunType(params.gunArray[i].type);

            tmp_json["coordinate"]["x"] = params.gunArray[i].x;
            tmp_json["coordinate"]["y"] = params.gunArray[i].y;
            tmp_json["coordinate"]["z"] = params.gunArray[i].z;

            tmp_json["volume"] = params.gunArray[i].volume;
            tmp_json["pressure"] = params.gunArray[i].pressure;
            tmp_json["shapeRatio"] = params.gunArray[i].shapeRatio;
            tmp_json["temperature"] = params.gunArray[i].temperature;

            jsonObject["gunArray"].push_back(tmp_json);
        }
    }


} // gund_json_parser

