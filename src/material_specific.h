#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "linalg.h"

using std::string;

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
    double momentum_error, probability_error, *probability, *energy;

    int momentum_samples, energy_samples;
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
    * Структура, определяющая, загружаем или вычисляем вероятности рассеяния
    *
*/

namespace config {
extern Files files;
extern Phonons phonons;
extern Probability probability;
extern Fields fields;
extern Bzone bzone;
extern Plot plot;
extern Model model;
};

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
    * Компоненты скорости
    *
*/
vec2 velocity(Point p);

/*
    *
    * Правые части уравнений движения
    *
*/
vec2 forces(Point p, double t);

/*
    *
    * Границы первой зоны Бриллюэна
    *
*/
double pmax(double psi);

/*
    *
    * Функция, приводящая квазиимпульс к первой зоне Бриллюэна
    *
*/
Point to_first_bz(Point p);

/*
    *
    * Возвращает точку с индексами (k,m)
    *
*/
Point point_k_m(int k, int m);