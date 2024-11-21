#include <gund_format_parser.h>

using gund_utility::toUpperCase;

namespace gund_format_parser {

    GundalfOutputParser::GundalfOutputParser(std::string const& inputFilename) {
        file.open(inputFilename);
        if (!file.is_open()) {
            throw std::runtime_error("File could not be opened");
        }
        // проверка типа данных
        checkType(inputFilename);
        // читаем файл
        parseFile();

        file.close();
    }

    void GundalfOutputParser::parseFile() {
        // чтение "шапки" файла
        parseFileHead();
        // чтение массива
        double keyWord;
        double i = 0;
        while (!file.eof() && i < sigParams.sampleNum) {
            file >> keyWord;
            data.push_back(keyWord);
            i++;
        }
        // считая ноль
        if (data.size() != sigParams.sampleNum) {
            throw std::runtime_error("Wrong number of points");
        }
    }

    void GundalfOutputParser::checkType(std::string const& inputFilename) {
        // самый минимальный размер строки с путем для файла = 4, так как .sig (.sgo; .amp; .phs)
        size_t minStr = 4;
        if (inputFilename.size() < minStr) {
            throw std::runtime_error("Wrong filename");
        }

        size_t length = inputFilename.size();
        std::string type = inputFilename.substr(length - 3, length);
        if (type == "sig") { fileType = sig; }
        else if (type == "sgo") { fileType = sgo; }
        else if (type == "amp") { fileType = amp; }
        else if (type == "phs") { fileType = phs; }
        else if (type == "flt") { fileType = flt; }
        else {
            throw std::runtime_error("Wrong filename");
        }
    }

    void GundalfOutputParser::parseFileHead() {
        std::string indName;
        // для # и =
        std::string specialStr;
        double value;

        double lineNum;
        // для амплитуды и фазы спектра не указываются единицы измерения в файле
        if (fileType == sig || fileType == sgo) { lineNum = 4; }
        else { lineNum = 3; }
        for (int i = 0; i < lineNum; i++) {
            file >> specialStr;
            if (specialStr != "#")
                throw std::runtime_error("Wrong file type");
            file >> indName;
            toUpperCase(indName);
            if (indName == "DT" || indName == "DF") {
                file >> specialStr >> value;
                sigParams.sampleInterval = value;
            }
            else if (indName == "IZ") {
                file >> specialStr >> value;
                if (value != 0)
                    throw std::runtime_error("Wrong start point");
                // To Do: добавить возможность изменения стартовой точки по времени
            }
            else if (indName == "NS") {
                file >> specialStr >> value;
                sigParams.sampleNum = value;
            }
            else if (indName == "UN") {
                // игнорируем единицы измерения
                file >> specialStr >> specialStr;
            }
        }
    }
} // namespace gund_format_paser

