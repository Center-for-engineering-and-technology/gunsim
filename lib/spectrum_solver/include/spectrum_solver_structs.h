#pragma once

#include <vector>

#include <gund_structs.h>

namespace spectrum_solver_structs {
    // постоянные параметры солвера
    struct SpectrumSolverParams {
        // набор пушек, объединенных общей точкой наблюдения
        std::vector<gund_structs::Gun> gunArray;
        // физические параметры среды для генерации сигнатур пушек
        gund_structs::PhysicalParameters physParams;
    };

    struct SpectrumSolverOptions {
        // параметры точки наблюдения
        gund_structs::ObservationPoint obsPoint;
        // наличие отражения от поверхности
        gund_structs::Reflection reflection;
        // характеристика фильтра
        gund_structs::Filter filter;
        // To Do: тут будет структура, содержащая характеристику фильтра, рассчитанную или выгруженную

        // шаг по времени/частоте из интерфейса, под который подстраиваются локальные значения для пушек
        gund_structs::SignatureParameters sigParams;

    };

    // To Do: использовать expected для хранения ошибок солвера
    struct GunSignalData {
        GunSignalData() = default;
        GunSignalData(std::vector<double> gunSignal_, gund_structs::SignatureParameters gunSigParams_) : gunSignal(gunSignal_), gunSigParams(gunSigParams_) {}
        std::vector<double> gunSignal;
        gund_structs::SignatureParameters gunSigParams;
    };

    struct ResultParameters {
        // значение peak-to-peak в Бар на метр и децибелах
        double peakToPeakBar;
        double peakToPeakdB;
        // погрешность в Бар-м
        double peakToPeakError{ 0.0 };

        double zeroToPeakBar;
        double zeroToPeakdB;

        double zeroToPeakError{ 0.0 };
        // и тд

    };

    struct SpectrumResult {
        ResultParameters rParams;

        std::vector<double> signal;
        std::vector<double> filteredSignal;
        std::vector<double> ampSpec;
        std::vector<double> phaseSpec;

        gund_structs::SignatureParameters sigParams;
        gund_structs::SignatureParameters specParams;
    };
} // namespace spectrum_solver_structs

