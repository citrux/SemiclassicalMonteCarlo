#pragma once

#include <string>
#include "linalg.h"

struct Files {
    bool load; // если 1, загружаем значения вероятностей из
    // файлов, иначе - вычисляем
    std::string probability, result;
};

struct Phonons {
    double beta, wla_max, wlo_max, T;
};

struct Fields {
    vec2 E0, E1, E2;

    double B0, B1, B2, omega1, omega2, phi1, phi2, phi;
};

struct Bzone {
    Point A, B, C, D;
    vec2 * basis, *dual_basis;
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
    double low;      // Минимальное значение переменной, от которой строим график
    double high;     // Максимальное значение переменной, от которой строим график
    double step;     // Шаг изменения переменной, в зависимости от которой строим
                     // график
    std::string var; // Переменная, в зависимости от которой строим график
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

void load_config(const std::string & filename);
void set_var(const std::string & varname, double value);

std::ostream & operator<<(std::ostream & os, const Fields & f);