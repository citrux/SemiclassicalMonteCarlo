#include <string>
#include "text_files.h"
#include "ini.h"
#include "logger.h"

using namespace std;

void array_from_file(const string & filename, Point * p_grid, double * Wer,
                     const Params & params) {
    logger(LOG_INFO, "load array from '" + filename + "'...");

    ifstream F_prob_in;
    F_prob_in.open(filename);
    for (int i = 0; i < (params.Nx + 1) * (params.Ny + 1); i++)
        F_prob_in >> p_grid[i].x >> p_grid[i].y >> Wer[i];
    F_prob_in.close();

    logger(LOG_OK, "\t[DONE]\n");
}

void array_to_file(const string & filename, Point * p_grid, double * Wer,
                   const Params & params) {
    logger(LOG_INFO, "write array to '" + filename + "'...");

    ofstream F_prob_out;
    F_prob_out.open(filename);
    for (int i = 0; i < (params.Nx + 1) * (params.Ny + 1); i++)
        F_prob_out << p_grid[i].x << "\t" << p_grid[i].y << "\t" << Wer[i]
                   << endl;
    F_prob_out.close();

    logger(LOG_OK, "\t[DONE]\n");
}

void write_logs(const string & filename_base, double var_value, double * px_log,
                double * py_log, int num_logs, const Params & params) {
    const string fn = filename_base + to_string(var_value) + ".txt";
    ostringstream tempStr;
    for (int i = 0; i < params.n_part; i++) {
        for (int j = 0; j < num_logs; j++)
            tempStr << px_log[i * num_logs + j] << "\t";
        tempStr << "\n";
        for (int j = 0; j < num_logs; j++)
            tempStr << py_log[i * num_logs + j] << "\t";
        tempStr << "\n";
    };
    ofstream file;
    file.open(fn);
    file << tempStr.str() << endl;
    file.close();
}

void results_to_file(const string & filename, int var_mas_count,
                     double * var_mas, double * result_value_mas_x,
                     double * result_value_mas_y, double * std_values_mas_x,
                     double * std_values_mas_y, double * result_av_time,
                     double * result_nOpt, double * result_nAc) {
    logger(LOG_INFO, "write results to '" + filename + "'...");

    ofstream f;
    f.open(filename);

    /*f << "Exc: " << params.Exc << "; Eyc: " << params.Eyc << "; wlo_max: " <<
    params.wlo_max << "; wla_max: " << params.wla_max << "; low: " << low << ";
    high: " << high << "; step: " << step << "; num_threads: " <<
    params.num_threads_openmp << endl;
    f << "Time info: " << total_time << endl;*/
    f << "#  E_y\tj_x\tj_y\tsigma_x\tsigma_y\ttau\tn_opt\tn_ac" << endl;
    for (int i = 0; i < var_mas_count; i++) {
        f << var_mas[i] << "\t" << result_value_mas_x[i] << "\t"
          << result_value_mas_y[i] << "\t" << std_values_mas_x[i] << "\t"
          << std_values_mas_y[i] << "\t" << result_av_time[i] << "\t"
          << result_nOpt[i] << "\t" << result_nAc[i] << endl;
    }
    f.close();

    logger(LOG_OK, "\t[DONE]\n");
}

double str_to_double(const string & str) {
    istringstream istream(str);
    double res;
    istream >> res;
    return res;
}

Point to_point(const string & str) {
    istringstream istream(str);
    Point p;
    istream >> p.x >> p.y;
    return p;
}

// на четных строках записаны комментарии, на нечетных строках записаны значения
// параметров
void load_config(const string & filename, Params & params, Graphic & graphic,
                 Probability & prob) {

    logger(LOG_INFO, "load config from '" + filename + "'...\t");

    INIReader reader(filename);

    if (reader.ParseError() < 0) {
        logger(LOG_ERROR, "[FAIL]\n");
        logger(LOG_ERROR, "Can't load '" + filename + "'\n");
        exit(1);
    }

    // разбор параметров
    params.beta = reader.GetReal("params", "beta", 0);
    params.wla_max =
        reader.GetReal("params", "wla_max", 0); // константа электрон-фононного
    // взаимодействия (акустические фононы)
    params.wlo_max =
        reader.GetReal("params", "wlo_max", 0); // константа электрон-фононного
    // взаимодействия (оптические фононы)
    params.Exc = reader.GetReal(
        "params", "Exc", 0); // постоянная составляющая электрического поля
                             // вдоль оси Ox
    params.Ex =
        reader.GetReal("params", "Ex", 0); // амплитуда переменной составляющей
    // электрического поля вдоль оси Ох
    params.Eyc = reader.GetReal(
        "params", "Eyc", 0); // постоянная составляющая электрического поля
                             // вдоль оси Oу
    params.Ey =
        reader.GetReal("params", "Ey", 0); // амплитуда переменной составляющей
    // электрического поля вдоль оси Оy
    params.H =
        reader.GetReal("params", "H", 0); // напряженность магнитного поля,
    // перпендикулярного к поверхности образца
    params.vx0 = reader.GetReal(
        "params", "vx0",
        0); // коэффициент в выражении для компоненты скорости (1)
    params.vy0 = reader.GetReal(
        "params", "vy0",
        0); // коэффициент в выражении для компоненты скорости (1)
    params.wx = reader.GetReal("params", "wx",
                               0); // частота электрического поля плоской волны,
    // поляризованной вдоль оси Ох
    params.wy = reader.GetReal("params", "wy",
                               0); // частота электрического поля плоской волны,
    // поляризованной вдоль оси Оу
    params.wy1 = reader.GetReal("params", "wy1", 1.0);
    params.wy2 = reader.GetReal("params", "wy2", 2.0);
    params.phi = reader.GetReal("params", "phi",
                                0); // относительная фаза волны, поляризованной
    // вдоль оси Оу, по сравнению с волной,
    // поляризованной вдоль оси Ох
    params.Ey1 = reader.GetReal("params", "Ey1", 0);
    params.Ey2 = reader.GetReal("params", "Ey2", 0);
    params.Nx =
        reader.GetInteger("params", "Nx", 1); // количество точек, на которые
    // разбивается первая зона Бриллюэна
    params.Ny = reader.GetInteger("params", "Ny", 1);
    params.px_max = reader.GetReal(
        "params", "px_max", 0); // максимальное значение компонент квазиимпульса
    params.py_max = reader.GetReal("params", "py_max", 0);
    params.A = to_point(
        reader.Get("params", "A", "0 0")); // координаты вершин параллелограмма,
    // ограничивающего первую зону Бриллюэна
    params.B = to_point(reader.Get("params", "B", "0 0"));
    params.D = to_point(reader.Get("params", "D", "0 0"));
    params.Newton_abs_error = reader.GetReal("params", "Newton_abs_error",
                                             0); // абсолютная ошибка вычисления
    // корня уравнения при
    // использовании метода Ньютона
    params.Newton_n_points = reader.GetInteger(
        "params", "Newton_n_points",
        0); // количество элементов, на которое делится отрезок, на
            // котором ищем решение уравнения методом Ньютона
    params.Simson_n = reader.GetInteger("params", "Simson_n",
                                        0); // количество элементов, на которое
    // разбивается область [0..2*Pi] при
    // интегрировании методом Симпсона
    params.dt = reader.GetReal(
        "params", "dt", 1e-4); // шаг по времени при моделировании (важен для
                               // метода Рунге-Кутты)
    params.all_time =
        reader.GetReal("params", "all_time", 0); // время всего моделирования
    params.num_threads_openmp = reader.GetInteger(
        "params", "num_threads_openmp",
        1); // количество ядер, на которое распараллеливается расчет
            // вероятности рассеяния и расчет среднего значения
            // плотности тока
    params.n_part =
        reader.GetInteger("params", "n_part", 1); // количество частиц n_part
    params.T = reader.GetReal("params", "T", -1); // температура
    params.Anorm = A_norm(params);
    params.max_prob = 1;

    graphic.low = reader.GetReal("plotting", "low", 0);
    graphic.high = reader.GetReal("plotting", "high", 0);
    graphic.step = reader.GetReal("plotting", "step", 0);
    graphic.num_var = reader.GetInteger("plotting", "num_var", -1);

    prob.flag_count_or_load =
        reader.GetInteger("mode", "flag_count_or_load", 0);
    prob.filename_opt = reader.Get("files", "optical", "UNKNOWN");
    prob.filename_ac = reader.Get("files", "acoustical", "UNKNOWN");
    prob.filename_res = reader.Get("files", "result", "UNKNOWN");

    logger(LOG_OK, "[DONE]\n");
}

int get_var_mas_count(const Graphic & graphic) {
    double low = graphic.low, high = graphic.high, step = graphic.step;
    int var_mas_count = (fabs(high - low) / step) + 1;
    return var_mas_count;
}

void set_var_mas(const Graphic & graphic, int var_mas_count, double * var_mas) {
    double low = graphic.low, step = graphic.step;
    for (int i = 0; i < var_mas_count; i++)
        var_mas[i] = low + step * i;
}

void var_value_graphic(int num_param, double var_value, Params & params) {
    switch (num_param) {
    case 0:
        params.beta = var_value;
        break;
    case 3:
        params.Exc = var_value;
        break;
    case 4:
        params.Ex = var_value;
        break;
    case 5:
        params.Eyc = var_value;
        break;
    case 6:
        params.Ey = var_value;
        break;
    case 7:
        params.H = var_value;
        break;
    case 14:
        params.phi = var_value;
        break;
    case 15:
        params.Ey1 = var_value;
        break;
    case 16:
        params.Ey2 = var_value;
        break;
    case 34:
        params.T = var_value;
        break;
    default:
        logger(LOG_ERROR, "Error!");
        exit(1);
    }
}