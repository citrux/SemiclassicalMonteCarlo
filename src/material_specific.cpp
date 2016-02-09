#include <string>
#include "material_specific.h"
#include "config.h"
#include "logger.h"
using namespace std;
double vf = 1e6, hbar = 6.5e-16, d = 2e-8, eps0 = 0.059, eps1 = 0.029,
       g = pow(eps1 / eps0, 2);

/*
    *
    * Выражение для энергетического спектра (в декартовых координатах)
    *
*/
double energy(Point p) {
    double a = vf * hbar / eps0 / d;
    return eps0 * (sqrt(1 + a * a * p.x * p.x) +
                   g * (1 - cos(p.y)) / sqrt(1 + a * a * p.x * p.x));
}
/*
    *
    * Выражение для энергетического спектра (в полярных координатах)
    *
*/
double energy_theta(double p, double theta) {
    return energy({p * cos(theta), p * sin(theta)});
}

/*
    *
    * Компоненты скорости
    *
*/
vec2 velocity(Point p) {
    double a = vf * hbar / eps0 / d;
    return {vf * a * p.x / sqrt(1 + a * a * p.x * p.x) *
                (1 - g * (1 - cos(p.y)) / (1 + a * a * p.y * p.y)),
            g * d * eps0 / hbar / sqrt(1 + a * a * p.x * p.x) * sin(p.y)};
}

/*
    *
    * Правые части уравнений движения
    *
*/
vec2 forces(Point p, double t) {
    double phi = config::fields.phi, phi1 = config::fields.phi1,
           phi2 = config::fields.phi2;
    double omega1 = config::fields.omega1, omega2 = config::fields.omega2;
    vec2 E0 = config::fields.E0;
    vec2 E1 = {config::fields.E1.x * cos(omega1 * t),
               config::fields.E1.y * cos(omega1 * t + phi1)};
    vec2 E2 = {config::fields.E2.x * cos(omega2 * t + phi),
               config::fields.E2.y * cos(omega2 * t + phi + phi2)};
    vec2 ov = {velocity(p).y, -velocity(p).x};
    return E0 + E1 + E2 + config::fields.H * ov;
}

/*
    *
    * Границы первой зоны Бриллюэна
    *
*/
double pmax(double theta) {
    // требуется найти пересечение луча с границей четырёхугольника
    Point O = {0, 0};
    // Считаем расстояние от начала координат до точки пересечения луча с
    // отрезками
    vec2 OA = config::bzone.A - O;
    vec2 OB = config::bzone.B - O;
    vec2 OC = config::bzone.C - O;
    vec2 OD = config::bzone.D - O;
    vec2 l = {cos(theta), sin(theta)};

    double res = -1;

    vec2 vs[5] = {OA, OB, OC, OD, OA};
    for (int i = 0; i < 4; ++i) {
        double p = cross(vs[i], vs[i + 1]) / cross(l, vs[i + 1] - vs[i]);

        if (p < 0)
            continue;

        if (res < 0 || p < res)
            res = p;
    }

    return res;
}

/*
    *
    * Функция, приводящая квазиимпульс к первой зоне Бриллюэна
    *
*/
Point to_first_bz(Point p) {

    vec2 pv = p - config::bzone.A;

    // находим разложение по базису, используя взаимный базис
    int n1 = floor(dot(pv, config::bzone.dual_basis[0]));
    int n2 = floor(dot(pv, config::bzone.dual_basis[1]));

    return p - n1 * config::bzone.basis[0] - n2 * config::bzone.basis[1];
}
