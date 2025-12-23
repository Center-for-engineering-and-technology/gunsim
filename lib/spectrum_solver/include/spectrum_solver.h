#pragma once

#include <variant>
#include <vector>

#include <direct_diag.h>
#include <energy_comp.h>
#include <fourier_solver.h>
#include <gun_model.h>
#include <gun_parameters_table.h>
#include <gund_format_parser.h>
#include <gund_json_parser.h>
#include <gund_structs.h>
#include <spectrum_solver_structs.h>

namespace spectrum_solver
{
using DropOutAndSpectrumResults = std::pair<spectrum_solver_structs::DropOutResult, spectrum_solver_structs::SpectrumResult>;

class SpectrumSolver
{
public:

    SpectrumSolver (const gun_parameters_table::GunMap& table, const spectrum_solver_structs::SpectrumSolverParams& params, const spectrum_solver_structs::SpectrumSolverOptions& options)
        : gunMap(table), specSolverParams(params), specSolverOptions(options)
    {
    }

    // запись адреса папки с данными относительно сборки (для GUI)
    void setDataPath();
    // запись адреса папки с данными относительно исходников
    void setDataPath (const std::string& path)
    {
        dataPath = path;
    }

    /// <summary>
    /// Расчет спектра от группы пневмоисточников
    /// </summary>
    /// <returns> Структура, содержащая сигнал, спектр и расчитанные инженерные параметры.
    /// В случае критической ошибки во время расчета - строка с текстом ошибки </returns>
    std::variant<spectrum_solver_structs::SpectrumResult, std::string> solve();
    /// <summary>
    /// Анализ отсева для группы пневмоисточников. Включает отсев одного, двух и(или) трех ПИ.
    /// Включает в себя расчет спектра от всей группы пневмоисточников - дополнительного применение solve() не требует.
    /// </summary>
    /// <returns> Пара структур: 1. DropOutResult - структура, содержащая информацию о всех проведенных отсевах,
    /// 2. SpectrumResult - структура, содержащая сигнал, спектр и расчитанные инженерные параметры.
    /// В случае критической ошибки во время расчета - строка с текстом ошибки </returns>
    std::variant<DropOutAndSpectrumResults, std::string> dropOut();

    const spectrum_solver_structs::SolverErrors& getErrors ()
    {
        return errors;
    }
    const std::list<std::string>& getMessages ()
    {
        return errors.messages;
    }

    /// <summary>
    /// Расчет диаграммы направленности
    /// </summary>
    /// <param name ="outDir"> Директория для выгрузки параметров </param>
    /// <param name ="angles"> Вектор углов (в градусах) для построения диаграммы направленности </param>
    /// <returns> Строка с текстом ошибки (при наличии) </returns>
    std::string computeDiagrams(const std::string& outDir, const std::unique_ptr<gun_model::GunModelResult>& gmr);
    std::string computeAzimutalDiagrams(const std::string& outDir, const std::unique_ptr<gun_model::GunModelResult>& gmr);

    /// <summary>
    /// Расчет сигнала одиночного пневмоисточника
    /// </summary>
    /// <param name ="gun"> Выбранный для моделирования пневмоисточник </param>
    /// <returns> Вектор с сигналом пневмоисточника (с тем же шагом по времени, что и указан в модели specSolverOptions) или строка с текстом ошибки (при наличии) </returns>
    std::variant<std::unique_ptr<std::vector<double>>, std::string> modelIndividualGunSignal(const gund_structs::Gun& gun);
    std::variant<std::unique_ptr<std::vector<double>>, std::string> modelIndividualGunSignal_diffModel(const gund_structs::Gun& gun, size_t i, const gun_model::GunModelResult& gun_model_result);
    // выдает расчитанные по дфииеренциальной модели сигнатуры, расчитанные в методе solve
    gun_model::GunModelResult getGunModelResult () const
    {
        return *gun_model_result;
    }

protected:

    // промежуточные функции для вычислений итогового результата

    // добавление фильтра - домножение рассчитанного ранее спектра на частотную характеристику фильтра в частотной области
    std::string addFilter(fourier_solver::ConcreteFourierSolver& solver, spectrum_solver_structs::SpectrumResult& result);
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

    // ссылка на таблицу параметров пушек
    const gun_parameters_table::GunMap& gunMap;
    const spectrum_solver_structs::SpectrumSolverParams& specSolverParams;
    spectrum_solver_structs::SpectrumSolverOptions specSolverOptions;
    spectrum_solver_structs::SpectrumSolverTemporaryData specSolverTmp;
    // расчитанные по дифференциальной модели сигнатуры
    std::unique_ptr<gun_model::GunModelResult> gun_model_result;
    // адрес файлов с характеристиками сигнала
    std::string dataPath;
    // хэш-таблица с относительными ошибками параметров (из файла)
    std::unordered_map<gund_structs::GunType, gund_structs::ErrorData> hashRelErrors;

    spectrum_solver_structs::SolverErrors errors;

    std::variant<std::unique_ptr<spectrum_solver_structs::GunSignalData>, std::string> getSignalFromFile(const gund_structs::Gun& gun);

    void fillRelErrors();
    // проверка того, что пушки (по номеру в массиве) не находятся в одной точке
    bool checkGunCoordinate(size_t vecPosition);
    // проверка несовпадения координаты точки наблюдения с координатами всех ПИ
    bool checkObservationCoordinate();

    /// <summary>
    /// Проверка типа точки наблюдения.
    /// </summary>
    /// <returns> Возвращает true при точке наблюдения на бесконечности </returns>
    bool checkObservationPoint();
    // проверка отражения
    bool checkReflection();
    // проверка наличия фильтра
    bool checkFilter();

    // добавляет к сигналу отражения от кабеля по определенному сценарию
    void addReflectionLineFromCable(std::vector<double>& gunCableSignal, const std::vector<double>& gunSignal, const double signalPath, const int power);

    // сохраняет индивидуальные сигналы пушек
    void fillIndividualSignals(const double gunTimeShift, const std::vector<double>& gunSignal);
    // обновляет значения во время прохода по циклу
    void updateSumValues(const gund_structs::Gun& gun, std::vector<double>& gunSignal);

    // считает расстояние между двумя пушками
    double computeDistanceBetweenGuns(const gund_structs::Gun& gunA, const gund_structs::Gun& gunB);
    // считает расстояние между пушкой и точкой наблюдения
    double computeDistanceFromGunToObsPoint(const gund_structs::Gun& gun);
    // считает задержку сигнала между источником и детектором
    double computeTimeShiftForSum(const gund_structs::Gun& gun);
    // считает задержку отраженного сигнала относительно прямого
    double computeRefTimeShift(const gund_structs::Gun& gun);
    // считает peak-to-peak и zero-to-peak параметры
    void computeToPeakValues(spectrum_solver_structs::SpectrumResult& result);
    //
    void computeContribAndErrorValues(spectrum_solver_structs::SpectrumResult& result);
    // считает центры: геометрический центр и центр давления
    void computeCenterValues(spectrum_solver_structs::SpectrumResult& result);
    // считает энергию и параметры, связанные с энергией
    void computeEnergyParams(spectrum_solver_structs::SpectrumResult& result);

    // заполнение соседних ПИ для каждой пушки
    void fillNeighborGuns();
    // заполнение максимальной глубины ПИ, выполняется для т. наблюдения на бесконечности
    void fillMaxDepth();

    // удаление сигнала единичной пушки из результирующего сигнала
    std::string removeGunSignal(std::vector<double>& totalSignal, const size_t gunIndex);
    // моделирование процедуры отсева по индексам отключенных пневмоисточников
    std::variant<spectrum_solver_structs::SpectrumResult, std::string> computeDropSignal(const std::vector<size_t>& dropGunIndexes, const std::vector<double>& fullSignal);
    //
    spectrum_solver_structs::DropOut compareDropAndOrigin(const std::vector<size_t>& dropGunIndexes, const spectrum_solver_structs::SpectrumResult& dropResult, const spectrum_solver_structs::SpectrumResult& allGunsResult);
};

} // namespace spectrum_solver