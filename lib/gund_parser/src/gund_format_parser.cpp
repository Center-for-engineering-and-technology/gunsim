#include "gund_format_parser.h"

#include <boost/scope_exit.hpp>
#include <string>

using gund_utility::toUpperCase;

namespace gund_format_parser
{

using Vec3 = Eigen::Vector3d;
using Cx = std::complex<double>;

GundalfOutputParser::GundalfOutputParser (std::string const& inputFilename)
{
    file.open(inputFilename);
    if (!file.is_open())
    {
        throw std::runtime_error("File " + inputFilename + " could not be opened");
    }
    // проверка типа данных
    checkType(inputFilename);
    // читаем файл
    parseFile();

    file.close();
}

void GundalfOutputParser::parseFile ()
{
    // чтение "шапки" файла
    parseFileHead();
    // чтение массива
    double keyWord;
    double i = 0;
    while (!file.eof() && i < sigParams.sampleNum)
    {
        file >> keyWord;
        data.push_back(keyWord);
        i++;
    }
    // считая ноль
    if (data.size() != sigParams.sampleNum)
    {
        throw std::runtime_error("Wrong number of points");
    }
}

void GundalfOutputParser::checkType (std::string const& inputFilename)
{
    // самый минимальный размер строки с путем для файла = 4, так как .sig (.sgo; .amp; .phs)
    size_t minStr = 4;
    if (inputFilename.size() < minStr)
    {
        throw std::runtime_error("Wrong filename");
    }

    size_t length = inputFilename.size();
    std::string type = inputFilename.substr(length - 3, length);
    if (type == "sig")
    {
        fileType = sig;
    }
    else if (type == "sgo")
    {
        fileType = sgo;
    }
    else if (type == "amp")
    {
        fileType = amp;
    }
    else if (type == "phs")
    {
        fileType = phs;
    }
    else if (type == "flt")
    {
        fileType = flt;
    }
    else
    {
        throw std::runtime_error("Wrong filename");
    }
}

void GundalfOutputParser::parseFileHead ()
{
    std::string indName;
    // для # и =
    std::string specialStr;
    double value;

    double lineNum;
    // для амплитуды и фазы спектра не указываются единицы измерения в файле
    if (fileType == sig || fileType == sgo)
    {
        lineNum = 4;
    }
    else
    {
        lineNum = 3;
    }
    for (int i = 0; i < lineNum; i++)
    {
        file >> specialStr;
        if (specialStr != "#")
            throw std::runtime_error("Wrong file type");
        file >> indName;
        toUpperCase(indName);
        if (indName == "DT" || indName == "DF")
        {
            file >> specialStr >> value;
            sigParams.sampleInterval = value;
        }
        else if (indName == "IZ")
        {
            file >> specialStr >> value;
            // if (value != 0)
            // throw std::runtime_error("Wrong start point");
            // To Do: добавить возможность изменения стартовой точки по времени
        }
        else if (indName == "NS")
        {
            file >> specialStr >> value;
            sigParams.sampleNum = value;
        }
        else if (indName == "UN")
        {
            // игнорируем единицы измерения
            file >> specialStr >> specialStr;
        }
    }
}

std::vector<gund_structs::RealSignature> parseNotionSourcesFile (
    const std::string& path
)
{
    std::ifstream file(path);
    if (!file.is_open())
    {
        throw std::runtime_error("File " + path + " could not be opened");
    }

    BOOST_SCOPE_EXIT(&file)
    {
        file.close();
    }
    BOOST_SCOPE_EXIT_END

    std::vector<gund_structs::RealSignature> result;

    // # dt =          0.0005
    // # ns = 2048
    // # nguns = 32
    //
    // # gun 1, vol=    45.0, x=   0.00, y= -10.40, z=   6.00

    std::string dummy;
    double dt;
    size_t ns, nguns;

    auto pass_values_n_times = [&] (int n_times)
    {
        for (auto i = 0; i < n_times; ++i)
        {
            file >> dummy;
            // printf("dummy: %s\n", dummy.c_str());
            // if (file.fail())
            //     printf("fail\n");
            // if (file.eof())
            //     printf("eof\n");
            // if (file.bad())
            //     printf("bad\n");
        }
    };

    pass_values_n_times(3);
    file >> dt;
    pass_values_n_times(3);
    file >> ns;
    pass_values_n_times(3);
    file >> nguns;

    result.reserve(nguns);

    std::string point_coord;
    for (size_t i = 0; i < nguns; ++i)
    {
        pass_values_n_times(2);
        file >> point_coord;
        point_coord.pop_back();
        int gun_num = stoi(point_coord);
        i = static_cast<size_t>(gun_num - 1);
        // printf("gun= %d ", gun_num);
        gund_structs::ObservationPoint point;
        point.observationType = gund_structs::ObservationPoint::ObservationPointType::SPECIFY;
        pass_values_n_times(3);
        file >> point_coord;
        // printf("x= %s ", point_coord.c_str());
        point_coord.pop_back();
        point.x = std::stod(point_coord);
        pass_values_n_times(1);
        file >> point_coord;
        // printf("y= %s ", point_coord.c_str());
        point_coord.pop_back();
        point.y = std::stod(point_coord);
        pass_values_n_times(1);
        file >> point_coord;
        // printf("z= %s\n", point_coord.c_str());
        point.z = std::stod(point_coord);

        gund_structs::NewSignatureParameters params(
            dt,
            std::move(point)
        );

        result.emplace_back(
            std::move(params),
            ns
        );

        auto& data = result.back().data;
        for (size_t j = 0; j < ns; ++j)
        {
            // 1 1       0.1325686
            pass_values_n_times(2);
            file >> data[j];
        }
    }
    return result;
}

void saveSignatureToSig (const std::string& path, gund_structs::RealSignature signature, const std::string& unit, double time_shift)
{
    std::ofstream file(path);
    if (!file.is_open())
    {
        throw std::runtime_error("File " + path + " could not be opened");
    }

    const auto& params = signature.params;
    const auto& data = signature.data;

    file << "# dt = " << params.sample_interval;
    file << "\n# iz = " << std::max(0, static_cast<int>(time_shift / params.sample_interval));
    file << "\n# ns = " << data.size();
    file << "\n# un = " << unit;
    for (auto val: data)
    {
        file << "\n"
             << val;
    }
}

} // namespace gund_format_parser
