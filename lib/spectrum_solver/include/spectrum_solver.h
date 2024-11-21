#pragma once

#include <vector>
#include <variant>

#include <gund_format_parser.h>
#include <gund_json_parser.h>
#include <gund_structs.h>
#include <fourier_solver.h>
#include <spectrum_solver_structs.h>

namespace spectrum_solver {
    class SpectrumSolver {
    public:
        SpectrumSolver(spectrum_solver_structs::SpectrumSolverParams params, spectrum_solver_structs::SpectrumSolverOptions options) :
            specSolverParams(params), specSolverOptions(options) {}

        // запись адреса папки с данными относительно сборки (для GUI)
        void setDataPath();
        // запись адреса папки с данными относительно исходников
        void setDataPath(const std::string& path) { dataPath = path; }

        // расчет спектра от ПИ (изменить тип вывода на класс выводных параметров, например)
        std::variant<spectrum_solver_structs::SpectrumResult, std::string> solve();
    protected:

        // промежуточные функции для вычислений итогового результата


        // добавление фильтра - домножение рассчитанного ранее спектра на частотную характеристику фильтра в частотной области
        void addFilter(fourier_solver::ConcreteFourierSolver& solver, spectrum_solver_structs::SpectrumResult& result);
        /// <summary>
        /// Добавление отражения от поверхности
        /// </summary>
        /// <returns> Строка с текстом ошибки (при наличии) </returns>
        std::string addReflection(const gund_structs::Gun& gun, std::vector<double>& gunSignal);
        /// <summary>
        /// Добавление влияния расстояния до точки наблюдения на исходный сигнал от ПИ при заданной точке наблюдения.
        /// </summary>
        /// <returns> Строка с текстом ошибки (при наличии) </returns>
        std::string addObservationPointInfluence(const gund_structs::Gun& gun, std::vector<double>& gunSignal);
        /// <summary>
        /// Добавление сигнала единичной пушки к общему сигналу, пришедшему на детектор
        /// </summary>
        /// <param name ="minTimeShift"> Минимальная задержка сигнала между источниками и детектором, принимается за нулевое значение по времени </param>
        void addNewGunSignal(double& minTimeShift, std::vector<double>& totalSignal, const gund_structs::Gun& gun, const std::vector<double>& gunSignal);

    private:
        spectrum_solver_structs::SpectrumSolverParams specSolverParams;
        spectrum_solver_structs::SpectrumSolverOptions specSolverOptions;
        // адрес файлов с характеристиками сигнала
        std::string dataPath;

        // To Do: сохранение полного спектра сигнала для дальнейшего ускорения процедуры подбора состояния пушек (отнимать спектр одной пушки быстрее чем считать заново)

        std::variant<std::unique_ptr<spectrum_solver_structs::GunSignalData>, std::string> getSignalFromFile(const gund_structs::Gun& gun);

        /// <summary>
        /// Проверка типа точки наблюдения.
        /// </summary>
        /// <returns> Возвращает true при точке наблюдения на бесконечности </returns>
        bool checkObservationPoint();
        // проверка отражения
        bool checkReflection();
        // проверка наличия фильтра
        bool checkFilter();

        // считает расстояние между пушкой и точкой наблюдения
        double computeDistanceFromGunToObsPoint(const gund_structs::Gun& gun);
        // считает задержку сигнала между источником и детектором
        double computeTimeShiftForSum(const gund_structs::Gun& gun);
        // считает задержку отраженного сигнала относительно прямого
        double computeRefTimeShift(const gund_structs::Gun& gun);
        // считает peak-to-peak и zero-to-peak параметры
        void computeToPeakValues(spectrum_solver_structs::SpectrumResult& result);
    };


} // namespace spectrum_solver