#include "text_files.h"
void array_from_file(char * filename, double * px_mas, double * py_mas,
                     double * Wer, Params params) {
    ifstream F_prob_in;
    F_prob_in.open(filename);
    for (int i = 0; i < (params.Nx + 1) * (params.Ny + 1); i++)
        F_prob_in >> px_mas[i] >> py_mas[i] >> Wer[i];
    F_prob_in.close();
}

void array_to_file(char * filename, double * px_mas, double * py_mas,
                   double * Wer, Params params) {
    ofstream F_prob_out;
    F_prob_out.open(filename);
    for (int i = 0; i < (params.Nx + 1) * (params.Ny + 1); i++)
        F_prob_out << px_mas[i] << "\t" << py_mas[i] << "\t" << Wer[i] << endl;
    F_prob_out.close();
}

void write_logs(char * filename_base, double var_value, double * px_log,
                double * py_log, int num_logs, Params params) {
    char fn[100] = "";
    sprintf(fn, "%s%.3f.txt", filename_base, var_value);
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

void results_to_file(char * filename_res, double * var_mas, int var_mas_count,
                     double * result_value_mas_x, double * result_value_mas_y,
                     double * std_values_mas_x, double * std_values_mas_y,
                     double * result_av_time, double * result_nOpt,
                     double * result_nAc) {
    ofstream f;
    f.open(filename_res);

    /*f << "Exc: " << params.Exc << "; Eyc: " << params.Eyc << "; wlo_max: " <<
    params.wlo_max << "; wla_max: " << params.wla_max << "; low: " << low << ";
    high: " << high << "; step: " << step << "; num_threads: " <<
    params.num_threads_openmp << endl;
    f << "Time info: " << total_time << endl;*/

    for (int i = 0; i < var_mas_count; i++) {
        f << var_mas[i] << "\t" << result_value_mas_x[i] << "\t"
          << result_value_mas_y[i] << "\t" << std_values_mas_x[i] << "\t"
          << std_values_mas_y[i] << "\t" << result_av_time[i] << "\t"
          << result_nOpt[i] << "\t" << result_nAc[i] << endl;
    }
    f.close();
}

int get_file_size(char * filename) {
    ifstream file_params_in;
    file_params_in.open(filename);
    file_params_in.seekg(0, ios::beg);
    int begin = file_params_in.tellg();
    file_params_in.seekg(0, ios::end);
    int end = file_params_in.tellg();
    file_params_in.close();
    return end - begin;
}

int get_file_num_lines(char * filename) {
    bool flag = true;
    char temp[1000];
    int curr;
    int num_lines = 0;
    int size = get_file_size(filename);
    ifstream file_params_in;
    file_params_in.open(filename);
    file_params_in.seekg(0, ios::beg);
    int begin = file_params_in.tellg();
    while (flag) {
        file_params_in.getline(temp, sizeof(temp));
        num_lines++;
        curr = (int) (file_params_in.tellg()) - begin;
        if (curr >= size)
            flag = false;
    };
    file_params_in.close();
    return num_lines;
}

double str_to_double(char * str) {
    std::istringstream istream(str);
    double res;
    istream >> res;
    return res;
}

// на четных строках записаны комментарии, на нечетных строках записаны значения
// параметров
void params_from_file(char * filename, Params * params, Graphic * graphic,
                      int * flag_count_or_load) {
    int size =
        get_file_size(filename); // определяем количество символов в файле
    int num_lines =
        get_file_num_lines(filename); // определяем количество строк в файле
    ifstream file_params_in; // создаем входной поток данных
    double * mas_param =
        new double[num_lines]; // массив, в который считываем параметры
    int num_param = 0; // номер параметра
    const int num_letters = 1000;
    char temp[num_letters]; // строка, в которую происходит считывание из файла
    char ** mas_comment = new char * [num_lines];
    for (int i = 0; i < num_lines; i++)
        mas_comment[i] = new char[num_letters];

    file_params_in.open(filename); // открываем файл
    file_params_in.seekg(0,
                         ios::beg); // устанавливаем указатель строк на начало
    int begin = file_params_in.tellg(); // позиция первого символа в файле
    bool flag = true; // флаг для определения конца файла
    bool even_odd = false; // флаг для определения - комментарий или параметр
    int curr = (int) (file_params_in.tellg()) -
               begin; // переменная для определения текущей позиции в файле

    while (flag) {
        file_params_in.getline(temp,
                               sizeof(temp)); // считываем одну строку из файла
        if (even_odd) // если строка нечетная (параметр)
        {
            mas_param[num_param] =
                str_to_double(temp); // считываем параметр в массив
            num_param++;
            even_odd = false; // следующая строка - комментарий
        } else // если строка четная (комментарий)
        {
            even_odd = true; // следующая строка - параметр
        }
        curr = (int) (file_params_in.tellg()) -
               begin; // вычисляем текущую позицию в файле
        if (curr >= size) // если дошли до конца файла, выходим
            flag = false;
    };
    file_params_in.close();

    // разбор параметров
    (*params).beta = mas_param[0];
    (*params).wla_max = mas_param[1]; // константа электрон-фононного
    // взаимодействия (акустические фононы)
    (*params).wlo_max = mas_param[2]; // константа электрон-фононного
    // взаимодействия (оптические фононы)
    (*params).Exc = mas_param[3]; // постоянная составляющая электрического поля
                                  // вдоль оси Ox
    (*params).Ex = mas_param[4]; // амплитуда переменной составляющей
    // электрического поля вдоль оси Ох
    (*params).Eyc = mas_param[5]; // постоянная составляющая электрического поля
                                  // вдоль оси Oу
    (*params).Ey = mas_param[6]; // амплитуда переменной составляющей
    // электрического поля вдоль оси Оy
    (*params).H = mas_param[7]; // напряженность магнитного поля,
    // перпендикулярного к поверхности образца
    (*params).vx0 =
        mas_param[8]; // коэффициент в выражении для компоненты скорости (1)
    (*params).vy0 =
        mas_param[9];             // коэффициент в выражении для компоненты скорости (1)
    (*params).wx = mas_param[10]; // частота электрического поля плоской волны,
    // поляризованной вдоль оси Ох
    (*params).wy = mas_param[11]; // частота электрического поля плоской волны,
    // поляризованной вдоль оси Оу
    (*params).wy1 = mas_param[12];
    (*params).wy2 = mas_param[13];
    (*params).phi = mas_param[14]; // относительная фаза волны, поляризованной
    // вдоль оси Оу, по сравнению с волной,
    // поляризованной вдоль оси Ох
    (*params).Ey1 = mas_param[15];
    (*params).Ey2 = mas_param[16];
    (*params).Nx = (int) (mas_param[17]); // количество точек, на которые
    // разбивается первая зона Бриллюэна
    (*params).Ny = (int) (mas_param[18]);
    (*params).px_max =
        mas_param[19]; // максимальное значение компонент квазиимпульса
    (*params).py_max = mas_param[20];
    (*params).Ax = mas_param[21]; // координаты вершин параллелограмма,
    // ограничивающего первую зону Бриллюэна
    (*params).Ay = mas_param[22];
    (*params).Bx = mas_param[23];
    (*params).By = mas_param[24];
    (*params).Dx = mas_param[25];
    (*params).Dy = mas_param[26];
    (*params).Newton_abs_error = mas_param[27]; // абсолютная ошибка вычисления
    // корня уравнения при
    // использовании метода Ньютона
    (*params).Newton_n_points =
        mas_param[28]; // количество элементов, на которое делится отрезок, на
                       // котором ищем решение уравнения методом Ньютона
    (*params).Simson_n = mas_param[29]; // количество элементов, на которое
    // разбивается область [0..2*Pi] при
    // интегрировании методом Симпсона
    (*params).dt = mas_param[30]; // шаг по времени при моделировании (важен для
                                  // метода Рунге-Кутты)
    (*params).all_time = mas_param[31]; // время всего моделирования
    (*params).num_threads_openmp =
        mas_param[32];                // количество ядер, на которое распараллеливается расчет
                                      // вероятности рассеяния и расчет среднего значения
                                      // плотности тока
    (*params).n_part = mas_param[33]; // количество частиц n_part
    (*params).T = mas_param[34];      // температура
    (*params).Anorm = A_norm(params);
    (*graphic).low = mas_param[35];
    (*graphic).high = mas_param[36];
    (*graphic).step = mas_param[37];
    (*graphic).num_var = (int) mas_param[38];
    *flag_count_or_load = (int) mas_param[39];
    delete[] mas_param;
}

int get_var_mas_count(Graphic graphic) {
    double low = graphic.low, high = graphic.high, step = graphic.step;
    int var_mas_count = (fabs(high - low) / step) + 1;
    return var_mas_count;
}

void set_var_mas(Graphic graphic, int var_mas_count, double * var_mas) {
    double low = graphic.low, high = graphic.high, step = graphic.step;
    for (int i = 0; i < var_mas_count; i++)
        var_mas[i] = low + step * i;
}

void var_value_graphic(int num_param, double var_value, Params * params) {
    switch (num_param) {
    case 0:
        (*params).beta = var_value;
        break;
    case 3:
        (*params).Exc = var_value;
        break;
    case 4:
        (*params).Ex = var_value;
        break;
    case 5:
        (*params).Eyc = var_value;
        break;
    case 6:
        (*params).Ey = var_value;
        break;
    case 7:
        (*params).H = var_value;
        break;
    case 14:
        (*params).phi = var_value;
        break;
    case 15:
        (*params).Ey1 = var_value;
        break;
    case 16:
        (*params).Ey2 = var_value;
        break;
    case 34:
        (*params).T = var_value;
        break;
    default:
        cout << "Error!";
    }
}

// на четных строках записаны комментарии, на нечетных строках записаны значения
// параметров
void comments_from_file(char * filename_source, char ** mas_comment,
                        int num_letters, int num_lines, int size) {
    ifstream file_params_in; // создаем входной поток данных
    double temp_double;
    char * temp = new char[num_letters]; // строка, в которую происходит
                                         // считывание из файла
    int num_comment = 0;
    file_params_in.open(filename_source); // открываем файл
    file_params_in.seekg(0,
                         ios::beg); // устанавливаем указатель строк на начало
    int begin = file_params_in.tellg(); // позиция первого символа в файле
    bool flag = true; // флаг для определения конца файла
    bool even_odd = false; // флаг для определения - комментарий или параметр
    int curr = (int) (file_params_in.tellg()) -
               begin; // переменная для определения текущей позиции в файле

    while (flag) {
        file_params_in.getline(mas_comment[num_comment],
                               num_letters); // считываем одну строку из файла
        num_comment++;
        curr = (int) (file_params_in.tellg()) -
               begin; // вычисляем текущую позицию в файле
        if (curr >= size) // если дошли до конца файла, выходим
            flag = false;
    };
    file_params_in.close();
}

void params_to_file(char * filename_dest, char * filename_source,
                    Params * params, Graphic * graphic, int total_time,
                    int flag_count_or_load) {
    int num_lines =
        get_file_num_lines(filename_source); // определяем количество строк в
    // исходном файле с параметрами
    int size = get_file_size(filename_source); // определяем количество символов
    // в исходном файле с параметрами
    int num_letters = 1000; // задаем количество символов в строке
    char ** mas_comment = new char * [num_lines]; // массив строк, каждая строка
    // в котором - соответствующий
    // комментарий
    for (int i = 0; i < num_lines; i++)
        mas_comment[i] = new char[num_letters];
    comments_from_file(filename_source, mas_comment, num_letters, num_lines,
                       size); // выбираем комментарии из исходного файла,
    // помещаем в массив строк
    // записываем комментарии и текущие значения параметров в файл
    ofstream file_params;
    file_params.open(filename_dest);
    file_params << (*params).beta << "\t" << mas_comment[0] << endl;
    file_params << (*params).wla_max << "\t" << mas_comment[2]
                << endl; // константа электрон-фононного взаимодействия
                         // (акустические фононы)
    file_params << (*params).wlo_max << "\t" << mas_comment[4]
                << endl; // константа электрон-фононного взаимодействия
                         // (оптические фононы)
    file_params
        << (*params).Exc << "\t" << mas_comment[6]
        << endl; // постоянная составляющая электрического поля вдоль оси Ox
    file_params << (*params).Ex << "\t" << mas_comment[8]
                << endl; // амплитуда переменной составляющей электрического
                         // поля вдоль оси Ох
    file_params
        << (*params).Eyc << "\t" << mas_comment[10]
        << endl; // постоянная составляющая электрического поля вдоль оси Oу
    file_params << (*params).Ey << "\t" << mas_comment[12]
                << endl; // амплитуда переменной составляющей электрического
                         // поля вдоль оси Оy
    file_params << (*params).H << "\t" << mas_comment[14]
                << endl; // напряженность магнитного поля, перпендикулярного к
                         // поверхности образца
    file_params << (*params).vx0 << "\t" << mas_comment[16]
                << endl; // коэффициент в выражении для компоненты скорости (1)
    file_params << (*params).vy0 << "\t" << mas_comment[18]
                << endl; // коэффициент в выражении для компоненты скорости (1)
    file_params << (*params).wx << "\t" << mas_comment[20]
                << endl; // частота электрического поля плоской волны,
                         // поляризованной вдоль оси Ох
    file_params << (*params).wy << "\t" << mas_comment[22]
                << endl; // частота электрического поля плоской волны,
                         // поляризованной вдоль оси Оу
    file_params << (*params).wy1 << "\t" << mas_comment[24] << endl;
    file_params << (*params).wy2 << "\t" << mas_comment[26] << endl;
    file_params << (*params).phi << "\t" << mas_comment[28]
                << endl; // относительная фаза волны, поляризованной вдоль оси
                         // Оу, по сравнению с волной, поляризованной вдоль оси
                         // Ох
    file_params << (*params).Ey1 << "\t" << mas_comment[30] << endl;
    file_params << (*params).Ey2 << "\t" << mas_comment[32] << endl;
    file_params << (*params).Nx << "\t" << mas_comment[34]
                << endl; // количество точек, на которые разбивается первая зона
                         // Бриллюэна
    file_params << (*params).Ny << "\t" << mas_comment[36] << endl;
    file_params << (*params).px_max << "\t" << mas_comment[38]
                << endl; // максимальное значение компонент квазиимпульса
    file_params << (*params).py_max << "\t" << mas_comment[40] << endl;
    file_params << (*params).Ax << "\t" << mas_comment[42]
                << endl; // координаты вершин параллелограмма, ограничивающего
                         // первую зону Бриллюэна
    file_params << (*params).Ay << "\t" << mas_comment[44] << endl;
    file_params << (*params).Bx << "\t" << mas_comment[46] << endl;
    file_params << (*params).By << "\t" << mas_comment[48] << endl;
    file_params << (*params).Dx << "\t" << mas_comment[50] << endl;
    file_params << (*params).Dy << "\t" << mas_comment[52] << endl;
    file_params << (*params).Newton_abs_error << "\t" << mas_comment[54]
                << endl; // абсолютная ошибка вычисления корня уравнения при
    // использовании метода Ньютона
    file_params << (*params).Newton_n_points << "\t" << mas_comment[56]
                << endl; // количество элементов, на которое делится отрезок, на
                         // котором ищем решение уравнения методом Ньютона
    file_params << (*params).Simson_n << "\t" << mas_comment[58]
                << endl; // количество элементов, на которое разбивается область
    // [0..2*Pi] при интегрировании методом Симпсона
    file_params << (*params).dt << "\t" << mas_comment[60]
                << endl; // шаг по времени при моделировании (важен для метода
                         // Рунге-Кутты)
    file_params << (*params).all_time << "\t" << mas_comment[62]
                << endl; // время всего моделирования
    file_params << (*params).num_threads_openmp << "\t" << mas_comment[64]
                << endl; // количество ядер, на которое распараллеливается
                         // расчет вероятности рассеяния и расчет среднего
                         // значения плотности тока
    file_params << (*params).n_part << "\t" << mas_comment[66]
                << endl; // количество частиц n_part
    file_params << (*params).T << "\t" << mas_comment[68]
                << endl; // температура
    file_params << (*graphic).low << "\t" << mas_comment[70] << endl; //
    file_params << (*graphic).high << "\t" << mas_comment[72] << endl;
    file_params << (*graphic).step << "\t" << mas_comment[74] << endl;
    file_params << (*graphic).num_var << "\t" << mas_comment[76] << endl;
    file_params << flag_count_or_load
                << "\t // если 1, вероятности загружали из файлов" << endl;
    file_params << total_time
                << "\t // общее время расчета для построения графика" << endl;
    file_params.close();
    for (int i = 0; i < num_lines; i++)
        delete[] mas_comment[i];
    delete[] mas_comment;
}

// выбираем названия файлов, в которых хранятся значения вероятности, а также
// помещаются результаты
void filename_from_file(char * filename_source, char * filename_params_after,
                        char * filename_opt, char * filename_ac,
                        char * filename_res) {
    ifstream file_with_filenames; // создаем входной поток данных
    file_with_filenames.open(filename_source); // открываем файл
    file_with_filenames.getline(filename_params_after,
                                1000); // считываем одну строку из файла
    cout << filename_params_after << endl;
    file_with_filenames.getline(filename_opt,
                                1000); // считываем одну строку из файла
    cout << filename_opt << endl;
    file_with_filenames.getline(filename_ac,
                                1000); // считываем одну строку из файла
    cout << filename_ac << endl;
    file_with_filenames.getline(filename_res,
                                1000); // считываем одну строку из файла
    cout << filename_res << endl;
    file_with_filenames.close();
}
