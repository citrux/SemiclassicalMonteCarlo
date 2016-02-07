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
double energy(Point p);
/*
    *
    * Выражение для энергетического спектра (в полярных координатах)
    *
*/
double energy_theta(double p, double theta);

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
double pmax(double theta);

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