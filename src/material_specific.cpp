#include "material_specific.h"
#include "logger.h"
/*
    *
    * Выражение для энергетического спектра (в декартовых координатах)
    *
*/
double energy(Point p, const Params & params) {
    return (1 - cos(p.x) * cos(p.y));
}
/*
    *
    * Выражение для энергетического спектра (в полярных координатах)
    *
*/
double energy_psi(double p, double psi, const Params & params) {
    return (1 - cos(p * cos(psi)) * cos(p * sin(psi)));
}

/*
    *
    * Производная энергии по модулю импульса (в полярных координатах)
    *
*/
double d_energy_psi(double p, double psi, const Params & params) {
    return (sin(p * cos(psi)) * cos(psi) * cos(p * sin(psi)) +
            cos(p * cos(psi)) * sin(p * sin(psi)) * sin(psi));
}

/*
    *
    * Компоненты скорости
    *
*/
double vx_fun(Point p, const Params & params) {
    return params.vx0 * sin(p.x) * cos(p.y);
}

double vy_fun(Point p, const Params & params) {
    return params.vy0 * cos(p.x) * sin(p.y);
}

/*
    *
    * Правые части уравнений движения
    *
*/

double right_x(Point p, double t, const Params & params) {
    return params.Exc;
}

double right_y(Point p, double t, const Params & params) {
    return params.Ey1 * cos(t * params.wy1) +
           params.Ey2 * cos(t * params.wy2 + params.phi);
}

/*
    *
    * Границы первой зоны Бриллюэна
    *
*/
double pmax(double psi, const Params & params) {
    // требуется найти пересечение луча с границей четырёхугольника
    Point C = params.B + (params.D - params.A);
    Point O = {0, 0};
    // Считаем расстояние от начала координат до точки пересечения луча с отрезками
    vec2 OA = params.A - O;
    vec2 OB = params.B - O;
    vec2 OC = C - O;
    vec2 OD = params.D - O;
    vec2 l = {cos(psi), sin(psi)};

    double res = -1;

    vec2 vs[5] = {OA, OB, OC, OD, OA};
    for (int i = 0; i < 5; ++i)
    {
        double p = cross(vs[i], vs[i+1]) / cross(l, vs[i+1] - vs[i]);

        if (p < 0) continue;

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
    vec2 b = params.B - params.A;
    vec2 d = params.D - params.A;

    // строим взаимный базис
    vec2 bc = ort(b - dot(b, ort(d)) * ort(d));
    bc = bc / dot(bc, b);

    vec2 dc = ort(d - dot(d, ort(b)) * ort(b));
    dc = dc / dot(dc, d);

    vec2 pr = p - params.A;

    // находим разложение по базису, используя взаимный базис
    int nb = floor(dot(pr, bc));
    int nd = floor(dot(pr, dc));

    return p - nb * b - nd * d;
}

/*
    *
    * Возвращает точку с индексами (k,m)
    *
*/
Point point_k_m(int k, int m, const Params & params) {
    double Ax = params.A.x;
    double Ay = params.A.y;
    double Bx = params.B.x;
    double By = params.B.y;
    double Dx = params.D.x;
    double Dy = params.D.y;

    double Kx = Ax + (Dx - Ax) * k / params.Nx;
    double Ky = Ay + (Dy - Ay) * k / params.Nx;
    double Mx = Ax + (Bx - Ax) * m / params.Ny;
    double My = Ay + (By - Ay) * m / params.Ny;
    double Ox =
        ((By - Ay) / (Bx - Ax) * Kx - (Dy - Ay) / (Dx - Ax) * Mx + My - Ky) /
        ((By - Ay) / (Bx - Ax) - (Dy - Ay) / (Dx - Ax));
    double Oy = (By - Ay) / (Bx - Ax) * (Ox - Kx) + Ky;

    return {Ox, Oy};
}

/*
    *
    * Массив координат точек в импульсном пространстве
    *
*/
void make_grid(Point * p_grid, const Params & params) {
    logger(LOG_INFO, "creating grid...");
    int Nx = params.Nx;
    int Ny = params.Ny;
    for (int i = 0; i < Nx + 1; i++)
        for (int j = 0; j < Ny + 1; j++) {
            p_grid[j + i * (Ny + 1)] = point_k_m(i, j, params);
        };
    logger(LOG_OK, "\t[DONE]\n");
}