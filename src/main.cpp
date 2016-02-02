#include "modelling_tools.h"
#include "logger.h"
using namespace std;

int main() {
    // задаем параметры задачи
    Params params;
    string config = "config.ini";

    logger(LOG_INFO, "> Loading config\n");
    load_config(config, params);

    set_probabilities(params);

    // Result_one_point result_one_point;

    /* массивы, необходимые для построения графика зависимости постоянной
     * составляющей тока от одной из переменных */
    /*
        // массив значений переменной, в зависимости от которой строим график
        // плотности тока
        int var_mas_count = get_var_mas_count(params);
        double * var_mas = new double[var_mas_count];
        set_var_mas(plot, var_mas_count, var_mas);

        double * result_value_mas_x =
            new double[var_mas_count]; // плотность тока jx, усредненная по
       времени
        // и по ансамблю
        double * result_value_mas_y =
            new double[var_mas_count]; // плотность тока jy, усредненная по
       времени
        // и по ансамблю
        double * std_values_mas_x =
            new double[var_mas_count]; // стандартная ошибка тока jx (по
       ансамблю)
        double * std_values_mas_y =
            new double[var_mas_count]; // стандартная ошибка тока jy (по
       ансамблю)
        double * result_av_time = new double[var_mas_count]; // среднее время
        // свободного пробега
        // (по ансамблю)
        double * result_nOpt = new double[var_mas_count]; // среднее количество
        // рассеяний на оптических
        // фононах (по ансамблю)
        double * result_nAc =
            new double[var_mas_count]; // среднее количество рассеяний на
        // акустических фононах (по ансамблю)
    */
    /* -----------------------------------------------------------------------------------------------
     */
    /*
logger(LOG_INFO, "> Start calculations\n");
time_t total_time = time(NULL);
for (int i = 0; i < var_mas_count; i++) // для всех значений переменной, в
// зависимости от которой строим
// график
{
    var_value_plot(plot.num_var, var_mas[i], params);
    result_one_point =
        one_plot_point(params, beta, p_grid, WerOpt, WerAc, var_mas[i], "");
    result_value_mas_x[i] = result_one_point.result_value_mas_x;
    result_value_mas_y[i] = result_one_point.result_value_mas_y;
    std_values_mas_x[i] = result_one_point.std_values_mas_x;
    std_values_mas_y[i] = result_one_point.std_values_mas_y;
    result_av_time[i] = result_one_point.result_av_time;
    result_nOpt[i] = result_one_point.result_nOpt;
    result_nAc[i] = result_one_point.result_nAc;
};
total_time = time(NULL) - total_time;
results_to_file(prob.filename_res, var_mas_count, var_mas,
                result_value_mas_x, result_value_mas_y, std_values_mas_x,
                std_values_mas_y, result_av_time, result_nOpt, result_nAc);
                */
}
