#pragma once

namespace gund_structs {

    struct SignatureParameters {
        SignatureParameters() = default;
        SignatureParameters(double sampleInterval_, double sampleNum_, double sampleMax_) :
            sampleInterval(sampleInterval_), sampleNum(sampleNum_), sampleMax(sampleMax_) {}
        // величина интервала по времени (с)
        double sampleInterval = 0.0005;
        // количество параметров в структуре сигнала
        double sampleNum = 1000;
        // максимальное значение для отображения данных (с)
        double sampleMax = 0.5;
        // To Do: проверить, что все файлы фильтров начинаются с нулевого значения
    };

    struct ObservationPoint {
        enum class ObservationPointType {
            INFINITE,
            SPECIFY
        };
        // тип точки наблюдения
        ObservationPointType observationType = ObservationPointType::INFINITE;
        // координаты точки наблюдения (м)
        double x;
        double y;
        double z;
    };

    struct Reflection {
        Reflection() = default;
        Reflection(double refCoef_) : refCoef(refCoef_) {};
        // коэффициент отражения от поверхности
        double refCoef = 0.;
    };

    struct PhysicalParameters {
        PhysicalParameters() = default;
        PhysicalParameters(double seaTemp_, double soundVelocity_) :
            seaTemp(seaTemp_), soundVelocity(soundVelocity_) {};
        // температура морской воды (С)
        double seaTemp = 10.;
        // скорость звука в воде (м/с)
        double soundVelocity = 1496.;
    };

    enum GunType {
        C1500, // 1500C
        C1900 // 1900C
    };

    struct Gun {
        
        GunType type = GunType::C1500;
        // координаты пушки (м)
        double x;
        double y;
        double z;
        // объем пушки (cuin)
        double volume;
        // отношение объема камеры к объему пушки [0., 1.]
        double shapeRatio = 1.;
        // давление пушки (psi)
        double pressure;
        // температура пушки (С)
        double temperature = 0.;
    };

    // полосовой фильтр, сгенерированный внутри
    struct BandpassFilter {

        enum AccessMode {
            OFF, // моделирование без фильтра
            EXTERNAL, // характеристика фильтра задана из файла
            INTERNAL // сгенерированный внутри фильтр
        };

        AccessMode mode = AccessMode::OFF;

        std::string filename;

        struct InternalFilter {
            // нижняя частота (Гц)
            double lowCutOffFrec;
            // верхняя частота (Гц)
            double highCutOffFrec;
            // наклон нижних частот (дБ/октаву)
            double lowCutOffSlope;
            // наклон вехних частот (дБ/октаву)
            double highCutOffSlope;
        };
    };

    struct Filter {

        enum FilterType {
            BANDPASS
        };

        FilterType filterType = FilterType::BANDPASS;

        // To Do: в дальнейшем использовать std::variant для выбора типа фильтра
        BandpassFilter bandpass;
    };
} // namespace gund_structs
