#include <ctime>
#include "modelling_tools.h"
#include "logger.h"

double random_uniform(unsigned int & x1, unsigned int & y1, unsigned int & z1,
                      unsigned int & w1) {
    unsigned int t1 = (x1 ^ (x1 << 11));
    x1 = y1;
    y1 = z1;
    z1 = w1;
    w1 = (w1 ^ (w1 >> 19)) ^ (t1 ^ (t1 >> 8));
    return ((double) w1) / ((double) UINT_MAX);
}

/*void RNGtest(double * mas, unsigned int * seed, int idx)
{
    unsigned int x1, y1, z1, w1; x1 = seed[idx]; y1 = 362436069; z1 = 521288629;
w1 = 88675123;
    mas[idx] = random_uniform(x1, y1, z1, w1);
}*/

void runge(Point & p, double t, const Params & params) {
    vec2 k1, k2, k3, k4;
    double dt = params.model.dt;

    k1 = forces(p, t, params);
    k2 = forces(p + dt / 2 * k1, t + dt / 2, params);
    k3 = forces(p + dt / 2 * k2, t + dt / 2, params);
    k4 = forces(p + dt * k3, t + dt, params);

    p += 1.0 / 6.0 * dt * (k1 + 2 * k2 + 2 * k3 + k4);
}

Point init_dist(unsigned int & x1, unsigned int & y1, unsigned int & z1,
                unsigned int & w1, const Params & params) {
    double psi;
    double p;
    double z;
    double f;
    Point point;
    do {
        psi = 2 * M_PI * random_uniform(x1, y1, z1, w1);
        p = pmax(psi, params) * random_uniform(x1, y1, z1, w1);
        z = random_uniform(x1, y1, z1, w1);
        f = distrib_function(p, psi, params);
    } while (z >= f);
    point.x = p * cos(psi);
    point.y = p * sin(psi);
    return point;
}

double Mean(double * arr, int count) {
    double result = 0.0;

    for (int i = 0; i < count; i++) {
        result += arr[i];
    }

    return result / count;
}

double Mean(unsigned int * arr, int count) {
    int result = 0;

    for (int i = 0; i < count; i++) {
        result += arr[i];
    }

    return ((double) result) / ((double) count);
}

double Std(double * arr, int count) {
    double sum = 0.0;
    double mean = Mean(arr, count);

    for (int i = 0; i < count; i++) {
        sum += pow((arr[i] - mean), 2);
    }

    return sqrt(sum / (count - 1));
}

void jobKernel(double * dev_average_value_x, double * dev_average_value_y,
               double * dev_average_time_array, unsigned int * dev_nAc,
               unsigned int * dev_nOpt,

               const Params & params, double beta, unsigned int rand_init_value,
               int idx,

               Point * p_grid, double * res_ac, double * res_opt,
               double * px_log, double * py_log, int num_logs) {
    // logger(LOG_INFO, "Start jobKernel\n");
    double * px_log_local = new double[num_logs];
    double * py_log_local = new double[num_logs];

    unsigned int x1, y1, z1, w1;
    x1 = rand_init_value;
    y1 = 362436069;
    z1 = 521288629;
    w1 = 88675123; // RNGInit

    bool flag;
    int iCount;
    double energy_value;
    double p_max;

    vec2 v;
    double value_x = 0.0;
    double value_y = 0.0;

    Point p = init_dist(x1, y1, z1, w1, params);
    double p1 = sqrt(p.x * p.x + p.y * p.y);
    double psi = atan2(p.y, p.x);

    double t = 0.0;
    double wsum = 0.0;
    double r1 = 0.0;
    double r2 = 0.0;
    double wla = 0.0;
    double wlo = 0.0;
    double r = -log(random_uniform(x1, y1, z1, w1));

    int n0 = 0;
    int nAc = 0;
    int nOpt = 0;

    double beta_opt = params.phonons.beta, beta_ac = 0.0;

    int t_num = 0;

    while (t < params.model.all_time) {
        v = velocity(p);

        value_x += v.x * params.model.dt;
        value_y += v.y * params.model.dt;

        runge(p, t, params); // решаем уравнения движения

        /* приводим импульс к зоне */
        p = to_first_bz(p, params);
        px_log_local[t_num] = p.x;
        py_log_local[t_num] = p.y;
        t += params.model.dt;

        r = -log(random_uniform(x1, y1, z1, w1));
        wsum += get_probability(p, params) * params.model.dt;

        if (wsum > r) {
            n0++; // наращиваем счетчик общего числа рассеяний
            wsum = 0.0;
            r = random_uniform(x1, y1, z1, w1);
            if (wla > r * wsum) {
                nAc++; // наращиваем счетчик рассеяний на акустических фононах
                energy_value = energy(p);
                flag = false;
                iCount = 0;
                while (!flag && iCount < 15) {
                    psi = 2 * M_PI *
                          random_uniform(x1, y1, z1, w1); // случайным образом
                    // разыгрываем фазу
                    // квазиимпульса
                    p_max = pmax(psi, params); // максимальное значение
                    // модуля квазиимпульса
                    // в направлении угла psi
                    p = momentums_with_energy_in_direction(psi, energy_value,
                                                           params)[0];
                    // если p1 существует, то flag = true, и мы правильно
                    // подобрали угол рассеяния, поэтому выходим из цикла

                    // если за 15 попыток не нашли решение, выходим из цикла
                    iCount++;
                }
            } else {
                if ((wlo > 0.0001) && (energy(p) >= beta_opt)) {
                    nOpt++; // наращиваем счетчик рассеяний на оптических
                            // фононах
                    energy_value = energy(p) - beta_opt;
                    flag = false;
                    iCount = 0;

                    while (!flag && iCount < 15) {
                        psi = 2 * M_PI * random_uniform(x1, y1, z1,
                                                        w1); // случайным
                                                             // образом
                        // разыгрываем фазу
                        // квазиимпульса
                        p_max = pmax(psi, params); // максимальное значение
                        // модуля квазиимпульса
                        // в направлении угла psi
                        p = momentums_with_energy_in_direction(
                            psi, energy_value, params)[0];
                        // если p1 существует, то flag = true, и мы правильно
                        // подобрали угол рассеяния, поэтому выходим из цикла

                        // если за 15 попыток не нашли решение, выходим из цикла
                        iCount++;
                    }
                }
            }
        }
        t_num++;
    }
    dev_average_value_x[idx] = value_x / t;
    dev_average_value_y[idx] = value_y / t;
    dev_average_time_array[idx] = t / n0;
    dev_nAc[idx] = nAc;
    dev_nOpt[idx] = nOpt;
    for (int i = 0; i < num_logs; i++) {
        px_log[i + idx * num_logs] = px_log_local[i];
        py_log[i + idx * num_logs] = py_log_local[i];
    };
    delete[] px_log_local;
    delete[] py_log_local;
    // logger(LOG_INFO, "End jobKernel\n");
}

Result_one_point one_plot_point(const Params & params, double beta,
                                Point * p_grid, double * WerOpt, double * WerAc,
                                double var_value,
                                const string & filename_base) {
    logger(LOG_INFO,
           "Calculate current density for " + to_string(var_value) + "...");

    time_t time_load = time(NULL);

    double * values_x =
        new double[params.model.particles]; // плотность постоянного
    // тока вдоль Ох (до
    // усреднения по ансамблю)
    double * values_y =
        new double[params.model.particles]; // плотность постоянного
    // тока вдоль Оу (до
    // усреднения по ансамблю)
    double * av_time =
        new double[params.model.particles]; // среднее время релаксации
    // (до усреднения по
    // ансамблю)
    unsigned int * nOpt =
        new unsigned int[params.model.particles]; // количество рассеяний на
    // оптических фононах
    unsigned int * nAc =
        new unsigned int[params.model.particles]; // количество рассеяний на
    // акустических фононах

    //Массивы-логи
    int num_logs = (int) (params.model.all_time / (params.model.dt)) + 1;
    double * px_log = new double[num_logs * (params.model.particles)];
    double * py_log = new double[num_logs * (params.model.particles)];

    // Инициализируем генератор случайных чисел
    srand(time(NULL));
    unsigned int * seed = new unsigned int[params.model.particles];
    for (int j = 0; j < params.model.particles; j++) {
        seed[j] = ((unsigned int) rand()) % 100000000 + 100000000;
    };

    omp_set_num_threads(params.model.threads);
#pragma omp parallel for
    for (int j = 0; j < params.model.particles; j++) {
        jobKernel(values_x, values_y, av_time, nAc, nOpt, params, beta, seed[j],
                  j, p_grid, WerAc, WerOpt, px_log, py_log, num_logs);
    };

    // Записываем значения компонент квазиимпульса в каждый момент времени для
    // каждой частицы
    // cout << "Writing logs for " << var_value << " .. ";
    // write_logs(filename_base, var_value, px_log, py_log, num_logs, *params);
    // cout << "done " << endl;

    // Помещаем результат расчета для одной точки графика в структуру
    Result_one_point result;
    result.result_value_mas_x = Mean(values_x, params.model.particles);
    result.result_value_mas_y = Mean(values_y, params.model.particles);
    result.std_values_mas_x = Std(values_x, params.model.particles);
    result.std_values_mas_y = Std(values_y, params.model.particles);
    result.result_av_time = Mean(av_time, params.model.particles);
    result.result_nOpt = Mean(nOpt, params.model.particles);
    result.result_nAc = Mean(nAc, params.model.particles);

    // Удаляем массивы - логи
    delete[] px_log;
    delete[] py_log;

    // Удаляем массивы, хранящие значения скорости, среднего времени релаксации
    // и т.д. до усреднения по ансамблю
    delete[] values_x;
    delete[] values_y;
    delete[] av_time;
    delete[] nOpt;
    delete[] nAc;

    logger(LOG_OK,
           "\t[DONE in " + to_string(time(NULL) - time_load) + " sec]\n");
    return result;
}