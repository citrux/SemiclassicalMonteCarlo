#include <string>
#include "material_specific.h"
#include "config.h"
#include "logger.h"
using namespace std;
double vf = 1e6, hbar = 6.5e-16, d = 2e-8, eps0 = 0.059, eps1 = 0.029,
       g = pow(eps1 / eps0, 2), c = 3e8;

/*
    *
    * Выражение для энергетического спектра (в декартовых координатах)
    *
*/
double energy(Point p) {
    double a = vf * hbar / eps0 / d;
    double root = sqrt(1 + a * a * p.x * p.x);
    return eps0 * (root + g * (1 - cos(p.y)) / root);
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
    * Градиент энергии в импульсном пространстве
    *
*/
vec2 energy_gradient(Point p) {
    double a = vf * hbar / eps0 / d;
    double b = 1 + a * a * p.x * p.x;
    double root = sqrt(b);
    double s, c;
    sincos(p.y, &s, &c);
    return {eps0 * a * a * p.x / root * (1 - g * (1 - c) / b),
            g * eps0 / root * s};
}
/*
    *
    * Скорость
    *
*/
vec2 velocity(Point p) {
    return (d / hbar / c) * energy_gradient(p);
}

/*
    *
    * Правые части уравнений движения
    *
*/
vec2 forces(Point p, double t) {
    using namespace config;
    vec2 force = fields.E0;
    force += vec2(fields.E1.x == 0 ? 0 : fields.E1.x * cos(fields.omega1 * t),
                  fields.E1.y == 0 ? 0 : fields.E1.y * cos(fields.omega1 * t +
                                                           fields.phi1));
    force += vec2(
        fields.E2.x == 0 ? 0
                         : fields.E2.x * cos(fields.omega2 * t + fields.phi),
        fields.E2.y == 0 ? 0 : fields.E2.y * cos(fields.omega2 * t +
                                                 fields.phi + fields.phi2));
    if (fields.H != 0) {
        vec2 v = velocity(p);
        v = vec2(v.y, -v.x);
        force += fields.H * v;
    }
    return force;
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
