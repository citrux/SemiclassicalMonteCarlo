#include "find_probability.h"
/*
    *
    * процедура проверки, меняет ли функция знак (сравниваем leftBound и
  *rightBound - значения
    * функции на границах рассматриваемого промежутка, если эти значения имеют
  *разные знаки, то на промежутке есть корень)
    *
*/
bool sign2(double leftBound, double rightBound) {
    if ((leftBound < 0 && rightBound > 0) ||
        (leftBound > 0 && rightBound < 0)) {
        return true;
    } else {
        return false;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Функции, используемые для нахождения вероятности рассеяния в полярной системе
// координат ////
///////////////////////////////////////////////////////////////////////////////////////////////

/*
    *
    * Функция для нахождения p1 из уравнения energy(p1,psi) = energy_value на
  *промежутке (p1_low, p1_high),
    * на котором гарантированно есть корень
    * Используется для определения границ первой зоны Бриллюэна (при построении
  *изоэнергетических поверхностей)
    *
*/
double Newton_psi_energy(double p1_low, double p1_high, double psi,
                         double energy_value, const Params & params) {
    double pn, pn1;
    bool usl1 = true;

    // выбираем в качестве начального приближения середину отрезка
    if (p1_low < p1_high) {
        pn = (p1_high - p1_low) / 2 + p1_low;
    } else if (p1_low > p1_high) {
        pn = (p1_low - p1_high) / 2 + p1_high;
    } else {
        return p1_low; // длина интервала равна нулю, поэтому выводим нижнюю
                       // границу
    }
    while (usl1) {
        pn1 = pn -
              (energy_psi(pn, psi, params) - energy_value) /
                  d_energy_psi(pn, psi, params);
        if (fabs(pn1 - pn) > params.Newton_abs_error) {
            pn = pn1;
        } else {
            usl1 = false;
        }
    }
    return pn1;
}

/*
    *
    * Функция, находящая первый корень уравнения energy(p1, psi)=energy_value на
  *промежутке (0, p1max).
    * Используется для построения изоэнергетических поверхностей
    *
*/
double apply_Newton_psi_energy(double psi, bool & flag, double energy_value,
                               double pmax, const Params & params) {
    double tempres = 0;
    // объявляем и заполняем массив значений p1, на котором будем искать решение
    // p1 берем в промежутке 0..p1max
    double * mas_p1 = new double[params.Newton_n_points + 1];
    double step = pmax / params.Newton_n_points;
    double lb, rb;

    for (int i = 0; i < params.Newton_n_points + 1; i++) {
        mas_p1[i] = i * step;
    }
    for (int i = 0; i < params.Newton_n_points; i++) {
        // если на интервале от mas_py1[i] до mas_py1[i+1] функция
        // eps(p1)-energy_value меняет знак,
        // то ищем на этом промежутке решение
        lb = energy_psi(mas_p1[i], psi, params) - energy_value;
        rb = energy_psi(mas_p1[i + 1], psi, params) - energy_value;
        if (sign2(lb, rb)) {
            tempres = Newton_psi_energy(mas_p1[i], mas_p1[i + 1], psi,
                                        energy_value, params);
            flag = true;
            goto L1;
        };
    }
L1:
    delete[] mas_p1;
    return tempres;
}

/*
    * Вычисление вероятности рассеяния электрона с импульсом p при помощи
    * билинейной интерполяции на основе нахождения пересечения вертикальной
 * прямой и плоскости,
    * построенной по трем точкам
    * p - импульс, для которого ищем вероятность. Должен быть приведен к первой
 * зоне Бриллюэна.
    * px_max, py_max - границы первой зоны Бриллюэна (предполагается
 * прямоугольная зона)
    * Nx, Ny - количество точек по соответствующим осям в сетке
    * angle - угол поворота системы координат (отсчитанный против часовой
 * стрелки) относительно
    *         оси сверхрешетки
    * points - массив значений импульса - точек сетки
    * Wer - массив значений вероятности, посчитанной на узлах выбранной сетки
*/
double getWer(double px, double py, double * px_mas, double * py_mas,
              double * Wer, const Params & params) {
    double Ax = params.Ax;
    double Ay = params.Ay;
    double Bx = params.Bx;
    double By = params.By;
    double Dx = params.Dx;
    double Dy = params.Dy;
    int Nx = params.Nx;
    int Ny = params.Ny;
    // высчитываем номер клеточки по x и y
    double K1x =
        ((By - Ay) / (Bx - Ax) * px - (Dy - Ay) / (Dx - Ax) * Ax - py + Ay) /
        ((By - Ay) / (Bx - Ax) - (Dy - Ay) / (Dx - Ax));
    double K1y = (By - Ay) / (Bx - Ax) * (K1x - px) + py;
    double M1x =
        ((Dy - Ay) / (Dx - Ax) * px - (By - Ay) / (Bx - Ax) * Ax - py + Ay) /
        ((Dy - Ay) / (Dx - Ax) - (By - Ay) / (Bx - Ax));
    double M1y = (Dy - Ay) / (Dx - Ax) * (M1x - px) + py;
    double stepAD = (Dx - Ax) / Nx;
    int px_i_temp = (int) floor((K1x - Ax) / stepAD);

    // !!! костыль !!!
    if (px_i_temp < 0) px_i_temp = 0;

    double stepAB = (Bx - Ax) / Ny;
    int py_i_temp = (int) floor((M1x - Ax) / stepAB);
    // определяем вершины клеточки, внутрь которой попала наша точка
    Point p1;
    p1.x = px_mas[py_i_temp + px_i_temp * (Ny + 1)];
    p1.y = py_mas[py_i_temp +
                  px_i_temp * (Ny + 1)]; // cout<<p1.x<<"\t"<<p1.y<<endl;
    Point p2;
    p2.x = px_mas[py_i_temp + px_i_temp * (Ny + 1) + 1];
    p2.y = py_mas[py_i_temp + px_i_temp * (Ny + 1) +
                  1]; // cout<<p2.x<<"\t"<<p2.y<<endl;
    Point p3;
    p3.x = px_mas[py_i_temp + px_i_temp * (Ny + 1) + Ny + 2];
    p3.y = py_mas[py_i_temp + px_i_temp * (Ny + 1) + Ny +
                  2]; // cout<<p3.x<<"\t"<<p3.y<<endl;
    Point p4;
    p4.x = px_mas[py_i_temp + px_i_temp * (Ny + 1) + Ny + 1];
    p4.y = py_mas[py_i_temp + px_i_temp * (Ny + 1) + Ny +
                  1]; // cout<<p4.x<<"\t"<<p4.y<<endl;

    /*
    определяем, где находится точка по отношению к диагонали клетки,
    проходящей через точки

   p1

           p4

   Уравнение этой диагонали

   (py-p1.y)/(p4.y-p1.y)=(px-p1.x)/(p4.x-p1.x),

   py=p1.y+(p4.y-p1.y)/(p4.x-p1.x)*(p.x-p1.x)
   */
    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    if (py > p1.y + (p4.y - p1.y) / (p4.x - p1.x) * (px - p1.x)) {
        // выше диагонали
        x1 = p1.x;
        y1 = p1.y;
        z1 = Wer[py_i_temp + px_i_temp * (Ny + 1)];
        x2 = p2.x;
        y2 = p2.y;
        z2 = Wer[py_i_temp + px_i_temp * (Ny + 1) + 1];
        x3 = p4.x;
        y3 = p4.y;
        z3 = Wer[py_i_temp + px_i_temp * (Ny + 1) + Ny + 2];
        /*
        1 p1	2 p2

                3 p4
        */
    } else {
        // ниже диагонали
        x1 = p1.x;
        y1 = p1.y;
        z1 = Wer[py_i_temp + px_i_temp * (Ny + 1)];
        x2 = p3.x;
        y2 = p3.y;
        z2 = Wer[py_i_temp + px_i_temp * (Ny + 1) + Ny + 1];
        x3 = p4.x;
        y3 = p4.y;
        z3 = Wer[py_i_temp + px_i_temp * (Ny + 1) + Ny + 2];
        /*
        1 p1

        2 p3    3 p4
        */
    }
    /* решаем систему линейных уравнений
        A*x1+B*y1+C*z1=1
        A*x2+B*y2+C*z2=1
        A*x3+B*y3+C*z3=1
        Находим A, B, C.
        Тогда z_int = (1-A*p.x-B*p.x)/C
    */
    double res;
    if (fabs(x3 * (y1 - y2) + x1 * (y2 - y3) + x2 * (-y1 + y3)) <
        0.0001) // условие, необходимое для обработки скачков
        res = Wer[py_i_temp + px_i_temp * (Ny + 1)];
    else
        res =
            ((-x3) * y2 * z1 + x2 * y3 * z1 + x3 * y1 * z2 - x1 * y3 * z2 -
             x2 * y1 * z3 + x1 * y2 * z3 +
             ((-y3) * z1 - y1 * z2 + y3 * z2 + y2 * (z1 - z3) + y1 * z3) * px +
             ((-x2) * z1 + x3 * z1 + x1 * z2 - x3 * z2 - x1 * z3 + x2 * z3) *
                 py) /
            (x3 * (y1 - y2) + x1 * (y2 - y3) + x2 * (-y1 + y3));
    if (res > 1000)
        res = 10;
    if (res < 0.00001)
        res = 0;
    return res;
    /*ofstream temp;
    temp.open("temp.txt");
    temp<<p.x<<" "<<p.y<<endl;
    temp<<p1.x<<" "<<p1.y<<endl;
    temp<<p2.x<<" "<<p2.y<<endl;
    temp<<p3.x<<" "<<p3.y<<endl;
    temp<<p4.x<<" "<<p4.y<<endl;
    temp.close();*/
}

/*
    *
    * Выражение, интеграл от которого берется методом Симпсона
    *
*/
double simpson_function(double psi, double px, double py, const Params & params) {
    bool flag = false;
    double result = 0.0;

    double p1max = pmax(psi, params);
    double energy_value = energy(px, py, params) - params.beta;

    double p = apply_Newton_psi_energy(psi, flag, energy_value, p1max, params);

    if (!flag)
        result = 0;
    else
        result = p / (fabs(d_energy_psi(p, psi, params)));

    return result;
}

/*
    *
    * Реализация метода Симпсона
    *
*/
double simpson(double px, double py, const Params & params) {
    double result = 0.0, a = 0, b = 2 * M_PI;
    double h = (b - a) / ((double) params.Simson_n);

    for (int i = 1; i < params.Simson_n - 1; i += 2) {
        double psi_1 = a + h * (i - 1);
        double psi_2 = a + h * i;
        double psi_3 = a + h * (i + 1);

        result += simpson_function(psi_1, px, py, params) +
                  4 * simpson_function(psi_2, px, py, params) +
                  simpson_function(psi_3, px, py, params);
    }

    return result * h / 3;
}

/*
    *
    * Собственно вычисление вероятности рассеяния
    *
*/
void full_probability_psi(double * px_mas, double * py_mas, double * res_mas,
                          const Params & params) {
    omp_set_num_threads(params.num_threads_openmp);
#pragma omp parallel for
    for (int i = 0; i < (params.Nx + 1) * (params.Ny + 1); i++) {
        res_mas[i] = simpson(px_mas[i], py_mas[i], params);
        // printf("%d \t %f \t %f \t %f\n", i, px_mas[i], py_mas[i],
        // res_mas[i]);
    };
}

double distrib_function(double p, double psi, const Params & params) {
    return exp((-1.0 / params.T) * energy_psi(p, psi, params)) * p;
}

double ditrib_function_int_p(double psi, const Params & params) {
    double result = 0.0, a = 0, b = pmax(psi, params);
    double h = (b - a) / ((double) ((params.Simson_n) * 5.0));

    for (int i = 1; i < (params.Simson_n) * 5.0 - 1; i += 2) {
        double p_1 = a + h * (i - 1);

        double p_2 = a + h * i;

        double p_3 = a + h * (i + 1);

        result += distrib_function(p_1, psi, params) +
                  4 * distrib_function(p_2, psi, params) +
                  distrib_function(p_3, psi, params);
    }
    return result * h / 3;
}

double A_norm(const Params & params) {
    double result = 0.0, a = 0, b = 2 * M_PI;
    double h = (b - a) / ((double) ((params.Simson_n) * 5.0));

    for (int i = 1; i < params.Simson_n * 5.0 - 1; i += 2) {
        double psi_1 = a + h * (i - 1);
        double psi_2 = a + h * i;
        double psi_3 = a + h * (i + 1);

        result += ditrib_function_int_p(psi_1, params) +
                  4 * ditrib_function_int_p(psi_2, params) +
                  ditrib_function_int_p(psi_3, params);
    }

    return result * h / 3;
}