#include "material_specific.h"
#include "logger.h"
/*
    *
    * Выражение для энергетического спектра (в декартовых координатах)
    *
*/
double energy(double px, double py, const Params & params) {
    return (1 - cos(px) * cos(py));
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
double vx_fun(double px, double py, const Params & params) {
    return params.vx0 * sin(px) * cos(py);
}

double vy_fun(double px, double py, const Params & params) {
    return params.vy0 * cos(px) * sin(py);
}

/*
    *
    * Правые части уравнений движения
    *
*/

double right_x(double px, double py, double t, const Params & params) {
    // return (params.Exc +
    // params.Ex*cos(t*params.wx)+params.H*vy_fun(px,py,params));
    return params.Exc;
}

double right_y(double px, double py, double t, const Params & params) {
    // return (params.Eyc + params.Ey*cos(t*params.wy+params.phi) -
    // params.H*vx_fun(px,py,params));
    return params.Ey1 * cos(t * params.wy1) +
           params.Ey2 * cos(t * params.wy2 + params.phi);
}

/*
    *
    * Границы первой зоны Бриллюэна
    *
*/
double pmax(double psi, const Params & params) {
    double res;
    if ((psi >= 0) && (psi < M_PI_2))
        res = M_PI / (1 + tan(psi)) / cos(psi);
    if (((psi >= M_PI_2) && (psi < M_PI)))
        res = -M_PI / (1 - tan(psi)) / cos(psi);
    if (((psi >= M_PI) && (psi < M_PI * 3 / 2)))
        res = -M_PI / (1 + tan(psi)) / cos(psi);
    if ((psi >= 3 * M_PI / 2) && (psi < 2 * M_PI))
        res = M_PI / (1 - tan(psi)) / cos(psi);
    return res;
}

/*
    *
    * Функция, приводящая квазиимпульс к первой зоне Бриллюэна
    *
*/
Point ToFirstBand(Point p0, const Params & params) {
    double px = p0.x;
    double py = p0.y;
    double p = sqrt(px * px + py * py);
    double phi = atan2(py, px);
    double OK = M_PI / sqrt(2.0); // орт e1 (модуль этого орта)
    double OM = M_PI / sqrt(2.0); // орт e2 (модуль этого орта)
    double OK1 =
        p * cos(phi - M_PI / 4); // проецируем исходный радиус-вектор на орт е1
    double OM1 =
        p * sin(phi - M_PI / 4); // проецируем исходный радиус-вектор на орт е2
    double pxi = floor(
        OK1 / OK); // целая часть от деления проекции исходного радиус-вектора
    // на орт е1 на модуль этого вектора
    double pyi = floor(
        OM1 / OM); // целая часть от деления проекции исходного радиус-вектора
    // на орт е2 на модуль этого вектора
    double pxf =
        OK1 -
        pxi * OK; // дробная часть от деления проекции исходного радиус-вектора
                  // на орт е1 на модуль этого вектора
    double pyf =
        OM1 -
        pyi * OM; // дробная часть от деления проекции исходного радиус-вектора
                  // на орт е1 на модуль этого вектора

    double pxN, pyN;
    if (fabs(pxi - floor(pxi / 2) * 2) < 0.000001) // если pxi - четная
        pxN = pxf;
    else
        pxN = -(OK - pxf);

    if (fabs(pyi - floor(pyi / 2) * 2) < 0.000001) // если pyi - четная
        pyN = pyf;
    else
        pyN = -(OM - pyf);
    double pxN1 = pxN / sqrt(2.0) - pyN / sqrt(2.0);
    double pyN1 = pxN / sqrt(2.0) + pyN / sqrt(2.0);
    Point p_res;
    p_res.x = pxN1;
    p_res.y = pyN1;
    return p_res;
}

/*
    *
    * Возвращает точку с индексами (k,m)
    *
*/
Point point_k_m(int k, int m, const Params & params) {
    double Ax = params.Ax;
    double Ay = params.Ay;
    double Bx = params.Bx;
    double By = params.By;
    double Dx = params.Dx;
    double Dy = params.Dy;

    double Kx = Ax + (Dx - Ax) * k / params.Nx;
    double Ky = Ay + (Dy - Ay) * k / params.Nx;
    double Mx = Ax + (Bx - Ax) * m / params.Ny;
    double My = Ay + (By - Ay) * m / params.Ny;
    double Ox =
        ((By - Ay) / (Bx - Ax) * Kx - (Dy - Ay) / (Dx - Ax) * Mx + My - Ky) /
        ((By - Ay) / (Bx - Ax) - (Dy - Ay) / (Dx - Ax));
    double Oy = (By - Ay) / (Bx - Ax) * (Ox - Kx) + Ky;
    Point res;
    res.x = Ox;
    res.y = Oy;
    return res;
}

/*
    *
    * Массив координат точек в импульсном пространстве
    *
*/
void make_grid(double * px_mas, double * py_mas, const Params & params) {
    logger(LOG_INFO, "creating grid...");
    int Nx = params.Nx;
    int Ny = params.Ny;
    Point temp_p;
    for (int i = 0; i < Nx + 1; i++)
        for (int j = 0; j < Ny + 1; j++) {
            temp_p = point_k_m(i, j, params);
            px_mas[j + i * (Ny + 1)] = temp_p.x;
            py_mas[j + i * (Ny + 1)] = temp_p.y;
        };
    logger(LOG_OK, "\t[DONE]\n");
}