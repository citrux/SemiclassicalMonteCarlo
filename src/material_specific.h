#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "linalg.h"
using std::string;

/*
    *
    * Структура, определяющая, загружаем или вычисляем вероятности рассеяния
    *
*/
struct Files {
    bool load; // если 1, загружаем значения вероятностей из
    // файлов, иначе - вычисляем
    string probability, result;
};

struct Phonons {
    double beta, wla_max, wlo_max, T;
};

struct Fields {
    vec2 E0, E1, E2;

    double H, omega1, omega2, phi1, phi2, phi;
};

struct Bzone {
    Point A, B, C, D;
};

struct Probability {
    double p_error, *probability, *energy;

    int p_points, n_integral, e_points;
};

struct Model {
    double dt, all_time;

    int threads, particles;
};
/*
    *
    * Структура, хранящая параметры построения графика
    *
*/
struct Plot {
    double low;  // Минимальное значение переменной, от которой строим график
    double high; // Максимальное значение переменной, от которой строим график
    double step; // Шаг изменения переменной, в зависимости от которой строим
                 // график
    int num_var; // Номер переменной, в зависимости от которой строим график
};

/*
    *
    * Структура, хранящая параметры задачи
    *
*/
struct Params {
    Files files;

    bool load; // загружать из файла или считать вероятности

    Phonons phonons;
    Probability probability;
    Fields fields;
    Bzone bzone;
    Plot plot;
    Model model;
};

void set_init_params(Params & params);

/*
    *
    * Выражение для энергетического спектра (в декартовых координатах)
    *
*/
double energy(Point p);
/*
    *
    * Выражение для энергетического спектра (в полярных координатах)
    *
*/
double energy_psi(double p, double psi);
/*
    *
    * Производная энергии по модулю импульса (в полярных координатах)
    *
*/
// double d_energy_psi(double p, double psi);

/*
    *
    * Компоненты скорости
    *
*/
vec2 velocity(Point p);

/*
    *
    * Правые части уравнений движения
    *
*/
vec2 forces(Point p, double t, const Params & params);

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