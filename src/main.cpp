#pragma offload_attribute(push, target(mic))
#include "modelling_tools.h"
#pragma offload_attribute(pop)
using namespace std;

int main()
{
	// задаем параметры задачи
	Params params;
	Graphic graphic;
	int flag_count_or_load; // если 1, загружаем значения вероятностей из файлов, иначе - вычисляем	
	set_init_params(&params);
	char filename_params[] = "params.txt";
	char fn[] = "filenames.txt";

	char filename_params_after[] = "params_after_calculations.txt";
	char filename_opt[] = "wopt.txt";
	char filename_ac[] = "wac.txt";
	char filename_res[] = "result.txt";
	
	filename_from_file(fn, filename_params_after, filename_opt, filename_ac, filename_res);
	params_from_file(filename_params, &params, &graphic, &flag_count_or_load);
	double beta = params.beta;		
	Result_one_point result_one_point;

	// задаем массивы, хранящие значения вероятности
	double* px_mas = new double[(params.Nx+1)*(params.Ny+1)];
	double* py_mas = new double[(params.Nx+1)*(params.Ny+1)];
	double* WerOpt = new double[(params.Nx+1)*(params.Ny+1)]; // массив значений вероятности рассеяния на оптических фононах
	double* WerAc = new double[(params.Nx+1)*(params.Ny+1)]; // массив значений вероятности рассеяния на акустических фононах
	
	// строим сетку
	points_mas(px_mas, py_mas, &params);

       time_t ttt;
	if(flag_count_or_load == 1) // если 1, загружаем вероятности из файла
	{
		cout << "Start loading... ";
		array_from_file(filename_opt, px_mas, py_mas, WerOpt, params);
		array_from_file(filename_ac, px_mas, py_mas, WerAc, params);
		cout << "done" << endl;
	}
	else // иначе вычисляем вероятности
	{	
		cout << "Start calculations of scatt. prob. on optical phonons...";
		ttt = time(NULL);
		params.beta = beta;
		#pragma offload target(mic:0) inout(WerOpt:length((params.Nx+1)*(params.Ny+1))) in(px_mas, py_mas: length((params.Nx+1)*(params.Ny+1))) in(params)
		{
			full_probability_psi(px_mas, py_mas, WerOpt, &params);
		};
		cout<<time(NULL)-ttt<<endl;	
		cout << "Writing the scatt. prob. on optical phonons in text file ... ";
		array_to_file(filename_opt, px_mas, py_mas, WerOpt, params);
		cout<<"done."<<endl<<" File name: " << filename_opt <<endl;
		
		ttt = time(NULL);
		cout << "Start calculations of scatt. prob. acoustical phonons...";
		params.beta = 0;
		#pragma offload target(mic:0) inout(WerAc:length((params.Nx+1)*(params.Ny+1))) in(px_mas, py_mas: length((params.Nx+1)*(params.Ny+1))) in(params)
		{
			full_probability_psi(px_mas, py_mas, WerAc, &params);
		}
		cout<<time(NULL)-ttt<<endl;
		cout << "Writing the scatt. prob. on acoustical phonons in text file ... ";
		array_to_file(filename_ac, px_mas, py_mas, WerAc, params);
		cout<<"done."<<endl<<" File name: " << filename_ac <<endl;
	};
	
/* массивы, необходимые для построения графика зависимости постоянной  составляющей тока от одной из переменных */

	// массив значений переменной, в зависимости от которой строим график плотности тока
	int var_mas_count = get_var_mas_count(graphic);
	double *var_mas = new double[var_mas_count];
	set_var_mas(graphic, var_mas_count, var_mas);
	
	double * result_value_mas_x = new double[var_mas_count]; // плотность тока jx, усредненная по времени и по ансамблю
	double * result_value_mas_y = new double[var_mas_count]; // плотность тока jy, усредненная по времени и по ансамблю
	double * std_values_mas_x = new double[var_mas_count];   // стандартная ошибка тока jx (по ансамблю)
	double * std_values_mas_y = new double[var_mas_count];   // стандартная ошибка тока jy (по ансамблю)
	double * result_av_time = new double[var_mas_count];     // среднее время свободного пробега (по ансамблю)
	double * result_nOpt = new double[var_mas_count];        // среднее количество рассеяний на оптических фононах (по ансамблю)
	double * result_nAc = new double[var_mas_count];         // среднее количество рассеяний на акустических фононах (по ансамблю)
	for(int i = 0; i < var_mas_count; i++)
	{
		result_value_mas_x[i] = 0.0;
		result_value_mas_y[i] = 0.0;
		std_values_mas_x[i] = 0.0;
		std_values_mas_y[i] = 0.0;
		result_av_time[i] = 0.0;
		result_nOpt[i] = 0.0;
		result_nAc[i] = 0.0;
	}
/* ----------------------------------------------------------------------------------------------- */
	time_t total_time = time(NULL), time_load  = time(NULL);
	for(int i = 0; i < var_mas_count; i++) // для всех значений переменной, в зависимости от которой строим график
	{		
		var_value_graphic(graphic.num_var, var_mas[i], &params);
		result_one_point = one_graphic_point(&params, beta, px_mas, py_mas, WerOpt, WerAc, var_mas[i], filename_params_after);
		cout << time(NULL) - time_load << endl;
		time_load = time(NULL);
		result_value_mas_x[i] = result_one_point.result_value_mas_x;
		result_value_mas_y[i] = result_one_point.result_value_mas_y;
		std_values_mas_x[i] = result_one_point.std_values_mas_x;
		std_values_mas_y[i] = result_one_point.std_values_mas_y;
		result_av_time[i] = result_one_point.result_av_time;
		result_nOpt[i] = result_one_point.result_nOpt;
		result_nAc[i] = result_one_point.result_nAc;
	};
	total_time = time(NULL) - total_time;
	results_to_file(filename_res, 
					var_mas, var_mas_count,
					result_value_mas_x, result_value_mas_y, 
					std_values_mas_x, std_values_mas_y,
					result_av_time,
					result_nOpt, result_nAc
					);
	params_to_file(filename_params_after, filename_params, &params, &graphic, (int)total_time, flag_count_or_load);
}
