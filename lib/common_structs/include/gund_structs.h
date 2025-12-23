#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

#include <string>
#include <vector>

#include <Eigen/Core>

namespace gund_structs
{

using NewCoordinate = Eigen::Vector3d;

struct Coordinate
{
    Coordinate() = default;
    Coordinate (double x_, double y_, double z_)
        : x(x_), y(y_), z(z_)
    {
    }
    Coordinate (const NewCoordinate& other)
        : x(other.x()), y(other.y()), z(other.z())
    {
    }

    double x{0};
    double y{0};
    double z{0};
};

struct SignatureParameters
{
    SignatureParameters() = default;
    SignatureParameters (double sampleInterval_, double sampleNum_, double sampleMax_)
        : sampleInterval(sampleInterval_), sampleNum(sampleNum_), sampleMax(sampleMax_)
    {
    }
    // величина интервала по времени (с или Гц)
    double sampleInterval = 0.0005;
    // количество параметров в структуре сигнала
    double sampleNum = 1000;
    // максимальное значение для отображения данных (с или Гц)
    double sampleMax = 0.5;
    // сдвиг 0.8 максимума первого пика сигнала от пушки вправо по времени
};

struct ObservationPoint
{
    enum class ObservationPointType
    {
        INFINITE,
        SPECIFY
    };
    // тип точки наблюдения
    ObservationPointType observationType = ObservationPointType::INFINITE;
    // координаты точки наблюдения (м)
    double x{0.};
    double y{0.};
    double z{0.};

    operator Eigen::Vector3d() const
    {
        return {x, y, z};
    }
};

struct Reflection
{
    struct Jovanovich
    {
        double meanWaveHeight{std::numeric_limits<double>::quiet_NaN()};
        double dominantFrequency{std::numeric_limits<double>::quiet_NaN()};
    };

    Reflection() = default;
    Reflection(double refCoef_)
        : refCoef(refCoef_) {};
    Reflection (Jovanovich config, double soundSpeed)
    {
        // TODO: в будущем хорошо бы вынести в расчетную часть,
        // т.к. в такой моедли коэффициент зависит от частоты и угла к нормали пришедшего сигнала
        // Пока частоту принимаем за доминирующую частоту вейвлета, а падение считаем нормальным
        // refco = - exp (- 2 * omega ^ 2 * sigma ^ 2 *cos^2 (theta) / soundSpeed^2)
        // sigma - средняя высота волны
        refCoef = -std::exp(
            -4.
            * M_PI
            * config.dominantFrequency
            * config.dominantFrequency
            * config.meanWaveHeight
            * config.meanWaveHeight
            / (soundSpeed
               * soundSpeed)
        );
    };

    // коэффициент отражения от поверхности
    double refCoef = 0.;
    double firstCableDepth = 0;
    double secondCableDepth = 0;
};

struct PhysicalParameters
{
    PhysicalParameters() = default;
    PhysicalParameters(double seaTemp_, double soundVelocity_)
        : seaTemp(seaTemp_), soundVelocity(soundVelocity_) {};
    // температура морской воды (С)
    double seaTemp = 10.;
    // скорость звука в воде (м/с)
    double soundVelocity = 1496.;
};

enum GunType
{
    B600,     // 600B
    C800,     // 800C
    C1500,    // 1500C
    LL1500,   // 1500LL
    C1900,    // 1900C
    DDHS1900, // 1900D-DHS
    LLX1900,  // 1900LLX
    LLXT1900, // 1900LLXT
    GUN2800,  // 2800
    LLX2800,  // 2800LLX
    APG8500,  // 8500APG
    GGUN,     // G-GUN
    GGUNII,   // G-GUNII
    GIGUN,    // GI-GUN
    SLEEVE,   // Sleeve
    SLEEVELL  // Sleevell
};

struct Gun
{

    std::string name;
    GunType type = GunType::C1500;
    // координаты пушки (м)
    double x = 0;
    double y = 0;
    double z = 0;
    // начало испуска сигнала (с)
    double delay = 0.;
    // объем пушки (cuin)
    double volume = 0;
    // отношение объема камеры к объему пушки [0., 1.]
    double shapeRatio = 1.;
    // давление пушки (psi)
    double pressure = 0;
    // температура пушки (С)
    double temperature = 0.;
};

// полосовой фильтр, сгенерированный внутри
struct BandpassFilter
{

    enum AccessMode
    {
        OFF,      // моделирование без фильтра
        EXTERNAL, // характеристика фильтра задана из файла
        INTERNAL  // сгенерированный внутри фильтр
    };

    AccessMode mode = AccessMode::OFF;

    std::string filename;
    // характеристика фильтра из файла
    std::vector<double> filterData;

    SignatureParameters filterSigParams;

    struct InternalFilter
    {
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

struct Filter
{

    enum FilterType
    {
        BANDPASS
    };

    FilterType filterType = FilterType::BANDPASS;

    // To Do: в дальнейшем использовать std::variant для выбора типа фильтра
    BandpassFilter bandpass;
};
// свойства спектральной полосы пропускания
struct Passband
{
    // нижняя граница полосы пропускания (Гц)
    double lowFrequency = 0.;
    // верхняя граница полосы пропускания (Гц)
    double highFrequency = 0.;
};

struct ErrorData
{
    // относительная ошибка параметра peak-to-peak для конкретной пушки
    double peakToPeakErr;
    // относительная ошибка нахождения bubble peak для конкретной пушки
    double bubblePeakErr;
    // относительная ошибка параметра bubble period для конкретной пушки
    double bubblePeriodErr;
};

struct DirectivityParameters
{
    // настройки стандартных диаграмм направленности
    double inlineAzimutAngle = 0;
    double crosslineAzimutAngle = 90;
    double maxDipAngle = 90;
    double lowerDB = 100;
    double higherDB = 210;
    double dipIncr = 1;
    double dipIncrForSignRepr = 5;
    // настройки азимутальных диаграмм направленности
    bool azimutal = false;
    double azimutalFreq1 = 30;
    double azimutalFreq2 = 60;
    double azimutalFreq3 = 90;
    double azimutalFreq4 = 120;
    double azimutalLowerDB = 60;
    double azimutalHigherDB = 210;
};

struct DropOutParameters
{
    enum DropOutMode
    {
        OFF = 0,
        SINGLE = 1,
        DOUBLE = 2,
        TRIPLE = 3
    };
    DropOutMode mode = DropOutMode::OFF;
    // максимально допустимое процентное расхождение при отсеве
    double maxPeakToPeakDrop = 0;
    // минимально допустимое значение primary-to-bubble
    double minPrimaryToBubble = 0;
    // максимально допустимое процентное расхождение величины primary-to-bubble
    double maxPrimaryToBubbleDrop = 0;
    // минимально допустимая корреляция
    double minNormXCorr = 1;
    // максимально допустимое снижение средней спектральной амплитуды [дБ]
    double maxAverageAmpDrop = 0;
    // максимально допустимая разница в спектральной амплитуде [дБ]
    double maxAmpDifference = 0;
};

class QFilter
{
public:

    static constexpr double kMinQValue{10.0};
    static constexpr double kMaxQValue{1000.0};

    static constexpr auto kMinTWTravelTime{0.01};
    static constexpr auto kMaxTWTravelTime{50.0};

    void setQValue (double value)
    {
        q_value = std::max(std::min(value, kMaxQValue), kMinQValue);
    }

    double getQValue () const
    {
        return q_value;
    }

    void setTWTravelTime (double value)
    {
        two_way_travel_time = std::max(std::min(value, kMaxTWTravelTime), kMinTWTravelTime);
    }

    double getTWTravelTime () const
    {
        return two_way_travel_time;
    }

private:

    double q_value{kMinQValue};
    double two_way_travel_time{kMinTWTravelTime};
};

class WienerFilter
{
public:

    static constexpr size_t kMinLen{5};
    static constexpr size_t kMaxLen{200};

    static constexpr size_t kMinGap{1};
    static constexpr size_t kMaxGap{200};

    static constexpr double kMinWhiteLightPercentage{0};
    static constexpr double kMaxWhiteLightPercentage{10};

    void setLen (size_t value)
    {
        wiener_len = std::max(std::min(value, kMaxLen), kMinLen);
    }

    size_t getLen () const
    {
        return wiener_len;
    }

    void setGap (size_t value)
    {
        wiener_gap = std::max(std::min(value, kMaxGap), kMinGap);
    }

    size_t getGap () const
    {
        return wiener_gap;
    }

    void setWhiteLightPercentage (double value)
    {
        wiener_wl = std::max(std::min(value, kMaxWhiteLightPercentage), kMinWhiteLightPercentage);
    }

    double getWhiteLightPercentage () const
    {
        return wiener_wl;
    }

    bool spike_mode{false};
    bool scale_mode{false};

private:

    size_t wiener_len{kMinLen};
    size_t wiener_gap{kMinGap};
    double wiener_wl{kMinWhiteLightPercentage};
    // bool energy_scale{false};
};

class NewBandPassFilter
{
public:

    static constexpr double kMinLowFreqCut{0};
    static constexpr double kMaxLowFreqCut{10000};

    static constexpr double kMinHighFreqCut{0};
    static constexpr double kMaxHighFreqCut{25000};

    static constexpr double kMinLowSlope{6};
    static constexpr double kMaxLowSlope{72};

    static constexpr double kMinHighSlope{0};
    static constexpr double kMaxHighSlope{350};

    virtual void setLowFreqCut (double value)
    {
        low_freq_cut = std::max(std::min(value, kMaxLowFreqCut), kMinLowFreqCut);
    }

    virtual double getLowFreqCut () const
    {
        return low_freq_cut;
    }

    virtual void setHighFreqCut (double value)
    {
        high_freq_cut = std::max(std::min(value, kMaxHighFreqCut), kMinHighFreqCut);
    }

    virtual double getHighFreqCut () const
    {
        return high_freq_cut;
    }

    virtual void setLowSlope (double value)
    {
        low_slope = std::max(std::min(value, kMaxLowSlope), kMinLowSlope);
    }

    virtual double getLowSlope () const
    {
        return low_slope;
    }

    virtual void setHighSlope (double value)
    {
        high_slope = std::max(std::min(value, kMaxHighSlope), kMinHighSlope);
    }

    virtual double getHighSlope () const
    {
        return high_slope;
    }

    bool minphase{false};

protected:

    double low_slope{kMinLowSlope};
    double low_freq_cut{kMinLowFreqCut};

    double high_slope{kMinHighSlope};
    double high_freq_cut{kMinHighFreqCut};
};

class AntiAliasFilter : public NewBandPassFilter
{
public:

    static constexpr double kMinNyqPrcntHighFreqCut{50};
    static constexpr double kMaxNyqPrcntHighFreqCut{100};

    static constexpr double kMinHighSlope{24};
    static constexpr double kMaxHighSlope{500};

    AntiAliasFilter ()
    {
        NewBandPassFilter::minphase = true;
        NewBandPassFilter::low_slope = 0;
        NewBandPassFilter::low_freq_cut = 0;
    }

    void setHighNyqPrcntFreqCut (double value)
    {
        nyq_prcnt = std::max(std::min(value, kMaxNyqPrcntHighFreqCut), kMinNyqPrcntHighFreqCut);
        NewBandPassFilter::high_freq_cut = std::max(std::min(0.5 * 0.01 * nyq_prcnt / out_tsamp, kMaxHighFreqCut), kMinHighFreqCut);
    }

    double getHighNyqPrcntFreqCut () const
    {
        return nyq_prcnt;
    }

    void setHighSlope (double value) override
    {
        high_slope = std::max(std::min(value, kMaxHighSlope), kMinHighSlope);
    }

    double getLowSlope () const override
    {
        return 0;
    }

    double getLowFreqCut () const override
    {
        return 0;
    }

    void setHighFreqCut (double /* value */) override
    {
    }

    void setOutTSamp (double value)
    {
        out_tsamp = value;
        NewBandPassFilter::high_freq_cut = std::max(std::min(0.5 * 0.01 * nyq_prcnt / out_tsamp, kMaxHighFreqCut), kMinHighFreqCut);
    }

    double getOutTSamp () const
    {
        return out_tsamp;
    }

private:

    double out_tsamp{0.0005};
    double nyq_prcnt{kMinNyqPrcntHighFreqCut};
};
} // namespace gund_structs
