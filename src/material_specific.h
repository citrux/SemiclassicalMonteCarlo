#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "linalg.h"

using std::string;

/*
    *
    * Выражение для энергетического спектра (в декартовых координатах)
    *
*/
double energy(const Point & p);
/*
    *
    * Выражение для энергетического спектра (в полярных координатах)
    *
*/
double energy_theta(double p, double theta);

vec2 energy_gradient(const Point & p);
/*
    *
    * Компоненты скорости
    *
*/
vec2 velocity(const Point & p);

/*
    *
    * Правые части уравнений движения
    *
*/
vec2 forces(const Point & p, double t);

/*
    *
    * Границы первой зоны Бриллюэна
    *
*/
double pmax(double theta);

/*
    *
    * Функция, приводящая квазиимпульс к первой зоне Бриллюэна
    *
*/
Point to_first_bz(const Point & p);