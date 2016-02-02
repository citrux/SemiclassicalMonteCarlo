#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
using std::string;
/*
    *
    * Структура, хранящая компоненты импульса
    *
*/
struct Point {
    double x;
    double y;
};

struct vec2 {
    double x;
    double y;
};

inline vec2 operator-(vec2 v) { return {-v.x, -v.y}; }

inline vec2 operator+(vec2 a, vec2 b) { return {a.x + b.x, a.y + b.y}; }

inline vec2 operator-(vec2 a, vec2 b) { return a + (-b); }

inline vec2 operator*(vec2 v, double scale) {
    return {scale * v.x, scale * v.y};
}

inline vec2 operator*(double scale, vec2 v) { return v * scale; }

inline vec2 operator/(vec2 v, double scale) { return v * (1. / scale); }

inline double dot(vec2 a, vec2 b) { return a.x * b.x + a.y * b.y; }

inline double cross(vec2 a, vec2 b) { return a.x * b.y - a.y * b.x; }

inline double len(vec2 v) { return sqrt(dot(v, v)); }

inline vec2 ort(vec2 v) { return v / len(v); }

inline vec2 operator-(Point end, Point start) {
    return {end.x - start.x, end.y - start.y};
}

inline Point operator+(Point start, vec2 shift) {
    return {start.x + shift.x, start.y + shift.y};
}

inline Point operator-(Point start, vec2 shift) { return start + (-shift); }

/*
    *
    * Структура, хранящая параметры задачи
    *
*/
struct Params {
    double beta;    // отношение энергии оптического фонона к ширине минизоны
    double wla_max; // константа электрон-фононного взаимодействия (акустические
                    // фононы)
    double wlo_max; // константа электрон-фононного взаимодействия (оптические
                    // фононы)

    double Exc; // постоянная составляющая электрического поля вдоль оси Ox
    double Ex;  // амплитуда переменной составляющей электрического поля вдоль
                // оси Ох
    double Eyc; // постоянная составляющая электрического поля вдоль оси Oу
    double Ey;  // амплитуда переменной составляющей электрического поля вдоль
                // оси Оy
    double H;   // напряженность магнитного поля, перпендикулярного к поверхности
                // образца
    double vx0; // коэффициент в выражении для компоненты скорости (1)
    double vy0; // коэффициент в выражении для компоненты скорости (1)
    double wx;  // частота электрического поля плоской волны, поляризованной
                // вдоль оси Ох
    double wy;  // частота электрического поля плоской волны, поляризованной
                // вдоль оси Оу
    double wy1;
    double wy2;
    double phi; // относительная фаза волны, поляризованной вдоль оси Оу, по
                // сравнению с
    // волной, поляризованной вдоль оси Ох

    double Ey1;
    double Ey2;

    int Nx; // количество точек, на которые разбивается первая зона Бриллюэна
    int Ny;
    double px_max; // максимальное значение компонент квазиимпульса
    double py_max;
    Point A; // координаты вершин параллелограмма, ограничивающего
    Point B; // первую зону Бриллюэна
    Point D;

    double Newton_abs_error; // абсолютная ошибка вычисления корня уравнения при
    // использовании метода Ньютона
    int Newton_n_points; // количество элементов, на которое делится отрезок, на
                         // котором ищем
    // решение уравнения методом Ньютона

    int Simson_n; // количество элементов, на которое разбивается область
                  // [0..2*Pi]
    // при интегрировании методом Симпсона

    double
        dt;          // шаг по времени при моделировании (важен для метода Рунге-Кутты)
    double all_time; // время всего моделирования

    int num_threads_openmp; // количество ядер, на которое распараллеливается
    // расчет вероятности рассеяния
    // и расчет среднего значения плотности тока
    int n_part; // количество частиц n_part

    double T; //Температура
    double Anorm;

    double max_prob;
};

/*
    *
    * Структура, определяющая, загружаем или вычисляем вероятности рассеяния
    *
*/
struct Probability {
    int flag_count_or_load; // если 1, загружаем значения вероятностей из
    // файлов, иначе - вычисляем
    string filename_opt;
    string filename_ac;
    string filename_res;
};

/*
    *
    * Структура, хранящая параметры построения графика
    *
*/
struct Graphic {
    double low;  // Минимальное значение переменной, от которой строим график
    double high; // Максимальное значение переменной, от которой строим график
    double step; // Шаг изменения переменной, в зависимости от которой строим
                 // график
    int num_var; // Номер переменной, в зависимости от которой строим график
};

void set_init_params(Params & params);

/*
    *
    * Выражение для энергетического спектра (в декартовых координатах)
    *
*/
double energy(Point p, const Params & params);
/*
    *
    * Выражение для энергетического спектра (в полярных координатах)
    *
*/
double energy_psi(double p, double psi, const Params & params);
/*
    *
    * Производная энергии по модулю импульса (в полярных координатах)
    *
*/
double d_energy_psi(double p, double psi, const Params & params);

/*
    *
    * Компоненты скорости
    *
*/
double vx_fun(Point p, const Params & params);
double vy_fun(Point p, const Params & params);

/*
    *
    * Правые части уравнений движения
    *
*/
double right_x(Point p, double t, const Params & params);
double right_y(Point p, double t, const Params & params);

/*
    *
    * Границы первой зоны Бриллюэна
    *
*/
double pmax(double psi, const Params & params);

/*
    *
    * Функция, приводящая квазиимпульс к первой зоне Бриллюэна
    *
*/
Point to_first_bz(Point p, const Params & params);

/*
    *
    * Возвращает точку с индексами (k,m)
    *
*/
Point point_k_m(int k, int m, const Params & params);

/*
    *
    * Массив координат точек в импульсном пространстве
    *
*/
void make_grid(Point * p_grid, const Params & params);