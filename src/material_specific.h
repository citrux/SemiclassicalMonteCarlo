#define _USE_MATH_DEFINES
#include <math.h>
/*
	*
	* Структура, хранящая компоненты импульса
	*
*/
struct __attribute__((target(mic:0))) Point
{
	double x;
	double y;
};
/*
	*
	* Структура, хранящая параметры задачи
	*
*/
struct __attribute__((target(mic:0))) Params
{
	double beta;    // отношение энергии оптического фонона к ширине минизоны
	double wla_max; // константа электрон-фононного взаимодействия (акустические фононы)
	double wlo_max; // константа электрон-фононного взаимодействия (оптические фононы)

	double Exc;     // постоянная составляющая электрического поля вдоль оси Ox
	double Ex;      // амплитуда переменной составляющей электрического поля вдоль оси Ох
	double Eyc;     // постоянная составляющая электрического поля вдоль оси Oу
	double Ey;      // амплитуда переменной составляющей электрического поля вдоль оси Оy
	double H;       // напряженность магнитного поля, перпендикулярного к поверхности образца 
	double vx0;     // коэффициент в выражении для компоненты скорости (1)
	double vy0;     // коэффициент в выражении для компоненты скорости (1)
	double wx;      // частота электрического поля плоской волны, поляризованной вдоль оси Ох 
	double wy;		// частота электрического поля плоской волны, поляризованной вдоль оси Оу
	double wy1;
	double wy2;
	double phi;		// относительная фаза волны, поляризованной вдоль оси Оу, по сравнению с
					// волной, поляризованной вдоль оси Ох

	double Ey1;
	double Ey2;

	int Nx;         // количество точек, на которые разбивается первая зона Бриллюэна 
	int Ny;	
	double px_max;  // максимальное значение компонент квазиимпульса
	double py_max;
	double Ax;      // координаты вершин параллелограмма, ограничивающего
	double Ay;		// первую зону Бриллюэна
	double Bx;
	double By;
	double Dx;
	double Dy;

	double Newton_abs_error; // абсолютная ошибка вычисления корня уравнения при использовании метода Ньютона
	int Newton_n_points;	 // количество элементов, на которое делится отрезок, на котором ищем 
							 // решение уравнения методом Ньютона

	int Simson_n;    // количество элементов, на которое разбивается область [0..2*Pi]
					 // при интегрировании методом Симпсона

	double dt;       // шаг по времени при моделировании (важен для метода Рунге-Кутты)
	double all_time; // время всего моделирования

	int num_threads_openmp; // количество ядер, на которое распараллеливается расчет вероятности рассеяния
							// и расчет среднего значения плотности тока 
	int n_part; // количество частиц n_part

	double T;		//Температура
	double Anorm;
	
	double max_prob;
};

/*
	*
	* Структура, хранящая параметры построения графика
	*
*/
struct Probability
{
	int flag_count_or_load; // если 1, загружаем значения вероятностей из файлов, иначе - вычисляем
	char* filename_opt;
	char* filename_ac;
	char* filename_res;
};

/*
	*
	* Структура, определяющая, загружаем или вычисляем вероятности рассеяния
	*
*/
struct Graphic
{
	double low; // Минимальное значение переменной, от которой строим график
	double high; // Максимальное значение переменной, от которой строим график
	double step; // Шаг изменения переменной, в зависимости от которой строим график
	int num_var; // Номер переменной, в зависимости от которой строим график
};



void set_init_params(Params *params);


/*
	*
	* Выражение для энергетического спектра (в декартовых координатах)
	*
*/
double __attribute__((target(mic:0))) energy(double px, double py, Params* params);
/*
	*
	* Выражение для энергетического спектра (в полярных координатах)
	*
*/
double __attribute__((target(mic:0))) energy_psi(double p, double psi, Params* params);
/*
	*
	* Производная энергии по модулю импульса (в полярных координатах)
	*
*/
double __attribute__((target(mic:0))) d_energy_psi(double p, double psi, Params* params);

/*
	*
	* Компоненты скорости
	*
*/
double __attribute__((target(mic:0))) vx_fun(double px, double py, Params* params);
double __attribute__((target(mic:0))) vy_fun(double px, double py, Params* params);

/*
	*
	* Правые части уравнений движения
	*
*/
double __attribute__((target(mic:0))) right_x(double px, double py, double t, Params* params);
double __attribute__((target(mic:0))) right_y(double px, double py, double t, Params* params);

/*
	*
	* Границы первой зоны Бриллюэна
	*
*/
double __attribute__((target(mic:0))) pmax(double psi, Params* params);

/*
	*
	* Функция, приводящая квазиимпульс к первой зоне Бриллюэна
	*
*/
Point __attribute__((target(mic:0))) ToFirstBand(Point p0, Params* params);

/*
	*
	* Возвращает точку с индексами (k,m)
	*
*/
Point __attribute__((target(mic:0))) point_k_m(int k, int m, Params* params);

/*
	*
	* Массив координат точек в импульсном пространстве
	*
*/
void __attribute__((target(mic:0))) points_mas(double* px_mas, double* py_mas, Params* params);