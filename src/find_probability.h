#include "material_specific.h"
#include <omp.h>
/*
    *
    * процедура проверки, меняет ли функция знак (сравниваем leftBound и
  *rightBound - значения
    * функции на границах рассматриваемого промежутка, если эти значения имеют
  *разные знаки, то на промежутке есть корень)
    *
*/
bool sign2(double leftBound, double rightBound);

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
                         double energy_value, const Params & params);

/*
    *
    * Функция, находящая первый корень уравнения energy(p1, psi)=energy_value на
  *промежутке (0, p1max).
    * Используется для построения изоэнергетических поверхностей
    *
*/
double apply_Newton_psi_energy(double psi, bool & flag, double energy_value,
                               double pmax, const Params & params);

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
              double * Wer, const Params & params);

/*
    *
    * Выражение, интеграл от которого берется методом Симпсона
    *
*/
double simpson_function(double psi, double px, double py, const Params & params);

/*
    *
    * Реализация метода Симпсона
    *
*/
double simpson(double px, double py, const Params & params);

/*
    *
    * Собственно вычисление вероятности рассеяния
    *
*/
void full_probability_psi(double * px_mas, double * py_mas, double * res_mas,
                          const Params & params);

double distrib_function(double p, double psi, const Params & params);

double A_norm(const Params & params);

double ditrib_function_int_p(double phi, const Params & params);
