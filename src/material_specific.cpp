#include "material_specific.h"
#include "logger.h"
/*
    *
    * Выражение для энергетического спектра (в декартовых координатах)
    *
*/
double energy(Point p) {
    return (1 - cos(p.x) * cos(p.y));
}
/*
    *
    * Выражение для энергетического спектра (в полярных координатах)
    *
*/
double energy_psi(double p, double psi) {
    return (1 - cos(p * cos(psi)) * cos(p * sin(psi)));
}

/*
    *
    * Производная энергии по модулю импульса (в полярных координатах)
    *
*/
double d_energy_psi(double p, double psi) {
    return (sin(p * cos(psi)) * cos(psi) * cos(p * sin(psi)) +
            cos(p * cos(psi)) * sin(p * sin(psi)) * sin(psi));
}

/*
    *
    * Компоненты скорости
    *
*/
vec2 velocity(Point p) {
    return {sin(p.x) * cos(p.y), cos(p.x) * sin(p.y)};
}

/*
    *
    * Правые части уравнений движения
    *
*/
vec2 forces(Point p, double t, const Params & params) {
    double phi = params.fields.phi, phi1 = params.fields.phi1, phi2 = params.fields.phi2;
    double omega1 = params.fields.omega1, omega2 = params.fields.omega2;
    vec2 E0 = params.fields.E0;
    vec2 E1 = {params.fields.E1.x * cos(omega1 * t),
               params.fields.E1.y * cos(omega1 * t + phi1)};
    vec2 E2 = {params.fields.E2.x * cos(omega2 * t + phi),
               params.fields.E2.y * cos(omega2 * t + phi + phi2)};
    vec2 ov = {velocity(p).y, -velocity(p).x};
    return E0 + E1 + E2 + params.fields.H * ov;
}

/*
    *
    * Границы первой зоны Бриллюэна
    *
*/
double pmax(double psi, const Params & params) {
    // требуется найти пересечение луча с границей четырёхугольника
    Point O = {0, 0};
    // Считаем расстояние от начала координат до точки пересечения луча с
    // отрезками
    vec2 OA = params.bzone.A - O;
    vec2 OB = params.bzone.B - O;
    vec2 OC = params.bzone.C - O;
    vec2 OD = params.bzone.D - O;
    vec2 l = {cos(psi), sin(psi)};

    double res = -1;

    vec2 vs[5] = {OA, OB, OC, OD, OA};
    for (int i = 0; i < 5; ++i) {
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
Point to_first_bz(Point p, const Params & params) {
    // базис обратной решётки
    vec2 b = params.bzone.B - params.bzone.A;
    vec2 d = params.bzone.D - params.bzone.A;

    // строим взаимный базис
    vec2 bc = ort(b - dot(b, ort(d)) * ort(d));
    bc = bc / dot(bc, b);

    vec2 dc = ort(d - dot(d, ort(b)) * ort(b));
    dc = dc / dot(dc, d);

    vec2 pr = p - params.bzone.A;

    // находим разложение по базису, используя взаимный базис
    int nb = floor(dot(pr, bc));
    int nd = floor(dot(pr, dc));

    return p - nb * b - nd * d;
}