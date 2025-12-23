#pragma once

#include <filesystem>
#include <list>
#include <memory>
#include <optional>
#include <variant>
#include <vector>

#include <gund_structs.h>

namespace spectrum_solver_structs
{
// постоянные параметры солвера
struct SpectrumSolverParams
{
    // набор пушек, объединенных общей точкой наблюдения
    std::vector<gund_structs::Gun> gunArray;
    // физические параметры среды для генерации сигнатур пушек
    gund_structs::PhysicalParameters physParams;
};
// опциональные параметры солвера
struct SpectrumSolverOptions
{
    // параметры точки наблюдения
    gund_structs::ObservationPoint obsPoint;
    // наличие отражения от поверхности
    gund_structs::Reflection reflection;
    // характеристика фильтра
    gund_structs::Filter filter;
    // инструментальный фильтр
    gund_structs::Filter pre_filter;

    std::optional<gund_structs::QFilter> q_filter;
    std::optional<gund_structs::WienerFilter> wiener_filter;
    std::variant<gund_structs::BandpassFilter, std::filesystem::path> bandpass_filter;
    std::optional<gund_structs::AntiAliasFilter> antialias_filter;

    // шаг по времени/частоте из интерфейса, под который подстраиваются локальные значения для пушек
    gund_structs::SignatureParameters sigParams;
    gund_structs::SignatureParameters output_sigParams;
    // TODO: вынести passband и время начала поиска пузыря отдельно в параметры отчета
    // спектральная полоса пропускания, необходимая для вычисления некоторых инженерных параметров
    gund_structs::Passband passband;
    // время начала поиска пузыря - нужно для расчета метрик
    double bubble_search_start_time{0.04}; // сек
    // параметры отсева
    gund_structs::DropOutParameters dropParams;
    // параметры направленности
    gund_structs::DirectivityParameters dirParams;
    // подключение дифф модели
    bool diffModel = false;
    // сдвиг выходной сигнатуры по времени
    double sampleTimeShift{std::numeric_limits<double>::quiet_NaN()};

    // Вариация контроллера источника (с)
    double sourceControllerVariation{0};
};

struct IndividualGunSignals
{
    double timeShift;
    std::unique_ptr<std::vector<double>> gunSignal;

    IndividualGunSignals() = default;
    IndividualGunSignals (double timeShift_, std::unique_ptr<std::vector<double>> gunSignal_)
        : timeShift(timeShift_), gunSignal(std::move(gunSignal_))
    {
    }
};

struct SpectrumSolverTemporaryData
{
    struct SumValues
    {
        double sumVolume = 0.;
        // суммарное накопленное значение sum_i(p-p_i * x_i), sum_i(p-p_i * y_i), sum_i(p-p_i * z_i)
        // для получения координат нужно разделить на p-p_sum
        gund_structs::Coordinate pressureCoordSum = {0., 0., 0.};
        double peakToPeakSum = 0.;
        // квадрат накопленного значения ошибок err_i^2 + err_{i+1}^2...
        // для получения реальногоо значения ошибок нужно извлечь корень
        double sumSqPeakToPeakError = 0.;
        double sumSqZeroToPeakError = 0.;
        double sumSqPrToBubError = 0.;
        double sumSqBubPeriodError = 0.;
        // peak-to-peak каждой активной пушки
        std::vector<double> individualPeakToPeak;
    };
    // суммарные значения обновляемых параметров
    SumValues sumValues;

    // сигналы каждой пушки (для быстрого вычисления drop-out)
    // шаг по времени и количество элементов совпадают (после интерполяции)
    // учтено влияние отражения и т.наблюдения (но нет смещения по времени)
    // смещение по времени хранится отдельно во избежание обрезания сигнала
    std::vector<IndividualGunSignals> individualSignals;
    // отсортированный набор ближайших соседей для каждой пушки
    std::vector<std::vector<size_t>> gunsNeighbors; // std::vector<std::unordered_set<size_t>> gunsNeighbors;
    // отсортированный набор пневмоисточников, которые необходимо игнорировать из-за близкого расположения
    std::vector<size_t> dublicateGuns; // std::unordered_set<size_t> dublicateGuns;

    double minTimeShift;
    double maxDepth = 0;
};

struct SolverErrors
{
    std::string criticalErrors;
    std::list<std::string> messages;
};

// To Do: использовать expected для хранения ошибок солвера
struct GunSignalData
{
    GunSignalData() = default;
    GunSignalData (std::vector<double> gunSignal_, gund_structs::SignatureParameters gunSigParams_)
        : gunSignal(gunSignal_), gunSigParams(gunSigParams_)
    {
    }
    std::vector<double> gunSignal;
    gund_structs::SignatureParameters gunSigParams;
};

struct DropOut
{
    // пневмоисточники, участвующие в процедуре отсева
    const gund_structs::Gun* firstGunDrop = nullptr;
    const gund_structs::Gun* secondGunDrop = nullptr;
    const gund_structs::Gun* thirdGunDrop = nullptr;

    double peakToPeak = 0;
    double peakToPeakPercentDrop = 0;
    double primaryToBubble = 0;
    double primaryToBubbleDrop = 0;
    double normXCorrDrop = 0;
    double averageAmpDrop = 0;
    double maxAmpDrop = 0;
    double freqMaxAmpDrop = 0;

    bool succeedDropOut = true;

    bool peakToPeakSucceed = true;
    bool primaryToBubbleSucceed = true;
    bool normXCorrSucceed = true;
    bool averageAmpSucceed = true;
    bool maxAmpSucceed = true;
};

struct DropOutResult
{
    // процент успешных анализов отсева
    double singleSuccessfulPercent{0.0};
    double doubleSuccessfulPercent{0.0};
    double tripleSuccessfulPercent{0.0};

    std::vector<DropOut> singleDropOut;
    std::vector<DropOut> doubleDropOut;
    std::vector<DropOut> tripleDropOut;
};

struct ResultParameters
{
    // значение peak-to-peak в Бар на метр и децибелах
    double peakToPeakBar;
    double peakToPeakdB;
    // погрешность в Бар-м
    double peakToPeakError{0.0};

    double zeroToPeakBar;
    double zeroToPeakdB;

    double zeroToPeakError{0.0};
    // отношение величин - безразмерная величина
    double primaryToBubble;

    double primaryToBubbleError{0.0};
    // время до первого пузыря (с)
    double bubblePeriodToFirstPeak;

    double bubblePeriodToFirstPeakError{0.0};

    double rmsPressureBar;
    double rmsPressuredB;

    double maximumSpectralValue;
    double averageSpectralValue;
    double maximumSpectralRipple;

    double totalAcousticEnergy = 0.;
    double totalPotentialEnergy = 0.;
    double totalAcousticEfficiency = 0.;

    // другие относительные значение

    double sumVolume;
    gund_structs::Coordinate geometricCenter;
    gund_structs::Coordinate pressureCenter;
    gund_structs::Coordinate energyCenter{0.0, 0.0, 0.0};
};
// вклад единичного источника
struct SingleSourceContribution
{
    std::string gunName;
    // вклад peak-to-peak единичной пушки в результирующий сигнал (%)
    double singlePeakToPeakPct;
    // вклад в отдельного ПИ в акустическую энергию
    double singleEnergy;
};

struct SpectrumResult
{
    ResultParameters rParams;
    std::vector<SingleSourceContribution> singleContrib;

    std::vector<double> signal;
    std::vector<double> filteredSignal;
    std::vector<double> ampSpec;
    std::vector<double> phaseSpec;

    std::vector<double> bandpassfilterAmpSpec;
    std::vector<double> bandpassfilterPhaseSpec;

    gund_structs::SignatureParameters sigParams;
    gund_structs::SignatureParameters specParams;
};
} // namespace spectrum_solver_structs
