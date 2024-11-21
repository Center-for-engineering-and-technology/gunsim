#pragma once

#include <vector>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <string>

#include <gund_structs.h>
#include <gund_utility.h>

// для загрузки сигнатуры источника из файла с типом данных, предоставленных программой Gundalf
namespace gund_format_parser {

    enum FileType {
        sig, // разрешение файла, характеризующее сигнал из ПИ в точке наблюдения
        sgo, // разрешение файла, характеризующее сигнал из ПИ в точке наблюдения после фильтрации
        amp, // разрешение файла, характеризующее амплитуду спектра ПИ в точке наблюдения после фильтрации
        phs, // разрешение файла, характеризующее фазу спектра ПИ в точке наблюдения после фильтрации
        flt // разрешение файла, предназначенное для характеристики фильтра
    };

    class GundalfOutputParser {
    public:
        GundalfOutputParser(std::string const& inputFilename);

        const std::vector<double>& getData() { return data; }
        const gund_structs::SignatureParameters& getSigParams() { return sigParams; }
    private:
        // проверяем тип файла
        void checkType(std::string const& inputFilename);
        // сохраняем данные из "шапки" файла
        void parseFileHead();
        // читаем весь файл
        void parseFile();

        // данные для "шапки" файла, содержащей информацию о количестве элементов в массиве, шаге по частоте/времени
        gund_structs::SignatureParameters sigParams;
        // данные, полученные из файла: сигнал, фильтрованный сигнал, амплитуда, фаза или характеристика фильтра
        std::vector<double> data;
        // тип полученных данных
        FileType fileType = FileType::sig;
        // файл, заполняется в конструкторе
        std::ifstream file;
    };
} // namespace gund_format_parser

