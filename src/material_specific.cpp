#include <string>
#include "material_specific.h"
#include "config.h"
#include "logger.h"
using namespace std;
/*
    *
    * Выражение для энергетического спектра (в декартовых координатах)
    *
*/
double energy(Point p) { return (1 - cos(p.x) * cos(p.y)); }
/*
    *
    * Выражение для энергетического спектра (в полярных координатах)
    *
*/
double energy_theta(double p, double theta) {
    return (1 - cos(p * cos(theta)) * cos(p * sin(theta)));
}

/*
    *
    * Компоненты скорости
    *
*/
vec2 velocity(Point p) { return {sin(p.x) * cos(p.y), cos(p.x) * sin(p.y)}; }

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
    // базис обратной решётки
    vec2 b = config::bzone.B - config::bzone.A;
    vec2 d = config::bzone.D - config::bzone.A;

    // строим взаимный базис
    vec2 bc = ort(b - dot(b, ort(d)) * ort(d));
    bc = bc / dot(bc, b);

    vec2 dc = ort(d - dot(d, ort(b)) * ort(b));
    dc = dc / dot(dc, d);

    vec2 pr = p - config::bzone.A;

    // находим разложение по базису, используя взаимный базис
    int nb = floor(dot(pr, bc));
    int nd = floor(dot(pr, dc));

    return p - nb * b - nd * d;
}
