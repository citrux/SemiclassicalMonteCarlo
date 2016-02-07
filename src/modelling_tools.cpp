#include <ctime>
#include "modelling_tools.h"
#include "logger.h"
#include "config.h"
#include "material_specific.h"
#include "find_probability.h"

double distrib_function(double p, double psi) {
    return exp(-energy_psi(p, psi) / config::phonons.T) * p;
}

double random_uniform(unsigned int & x1, unsigned int & y1, unsigned int & z1,
                      unsigned int & w1) {
    unsigned int t1 = (x1 ^ (x1 << 11));
    x1 = y1;
    y1 = z1;
    z1 = w1;
    w1 = (w1 ^ (w1 >> 19)) ^ (t1 ^ (t1 >> 8));
    return ((double) w1) / ((double) UINT_MAX);
}

double * angle_int_distrib;

double inverse_function_angle(double prob) {
    int n = 1000;
    int m = 500;
    if (angle_int_distrib == nullptr) {
        angle_int_distrib = new double[n+1];
        angle_int_distrib[0] = 0;
        for (int i = 1; i <= n; ++i)
        {
            double psi = 2 * M_PI * i / n;
            double momentum_int = 0;
            double step = pmax(psi) / m;
            for (int j = 0; j < m; ++j)
                momentum_int += distrib_function((j + 0.5) * step, psi) * step;

            angle_int_distrib[i] = angle_int_distrib[i-1] + momentum_int;
        }
    }

    // нормировка
    for (int i = 1; i <= n; ++i) {
        angle_int_distrib[i] /= angle_int_distrib[n];
    }

    // бинпоиск
    int i = 0, j = n;
    while (j - i > 1) {
        int k = (i + j) / 2;
        if (angle_int_distrib[k] > prob)
            j = k;
        else
            i = k;
    }

    double w = (prob - angle_int_distrib[i]) / (angle_int_distrib[j] - angle_int_distrib[i]);
    return (i + w) * 2 * M_PI / n;
}

double inverse_function_momentum(double psi, double prob) {
    int n = 1000;
    double * int_distrib = new double[n + 1];
    double p = pmax(psi), step = p / n;

    int_distrib[0] = 0;
    // табулируем интегральную функцию
    for (int i = 1; i <= n; ++i) {
        int_distrib[i] =
            int_distrib[i - 1] + distrib_function((i - 0.5) * step, psi) * step;
    }

    // нормировка
    for (int i = 1; i <= n; ++i) {
        int_distrib[i] /= int_distrib[n];
    }

    // бинпоиск
    int i = 0, j = n;
    while (j - i > 1) {
        int k = (i + j) / 2;
        if (int_distrib[k] > prob)
            j = k;
        else
            i = k;
    }
    double w = (prob - int_distrib[i]) / (int_distrib[j] - int_distrib[i]);
    return (i + w) * step;
}

void runge(Point & p, double t) {
    vec2 k1, k2, k3, k4;
    double dt = config::model.dt;

    k1 = forces(p, t);
    k2 = forces(p + dt / 2 * k1, t + dt / 2);
    k3 = forces(p + dt / 2 * k2, t + dt / 2);
    k4 = forces(p + dt * k3, t + dt);

    p += 1.0 / 6.0 * dt * (k1 + 2 * k2 + 2 * k3 + k4);
}

Point init_dist(unsigned int & x1, unsigned int & y1, unsigned int & z1,
                unsigned int & w1) {
    double r = random_uniform(x1, y1, z1, w1);
    double psi = inverse_function_angle(r);
    double z = random_uniform(x1, y1, z1, w1);
    double p = inverse_function_momentum(psi, z);
    return {p * cos(psi), p * sin(psi)};
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

               double beta, unsigned int rand_init_value, int idx,

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

    Point p = init_dist(x1, y1, z1, w1);
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

    double beta_opt = config::phonons.beta, beta_ac = 0.0;

    int t_num = 0;

    while (t < config::model.all_time) {
        v = velocity(p);

        value_x += v.x * config::model.dt;
        value_y += v.y * config::model.dt;

        runge(p, t); // решаем уравнения движения

        /* приводим импульс к зоне */
        p = to_first_bz(p);
        px_log_local[t_num] = p.x;
        py_log_local[t_num] = p.y;
        t += config::model.dt;

        r = -log(random_uniform(x1, y1, z1, w1));
        double e = energy(p);
        wsum += (config::phonons.wla_max * get_probability(e - beta_ac) +
                 config::phonons.wlo_max * get_probability(e - beta_opt)) *
                config::model.dt;

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
                    p_max = pmax(psi); // максимальное значение
                    // модуля квазиимпульса
                    // в направлении угла psi
                    p = momentums_with_energy_in_direction(psi,
                                                           energy_value)[0];
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
                        p_max = pmax(psi); // максимальное значение
                        // модуля квазиимпульса
                        // в направлении угла psi
                        p = momentums_with_energy_in_direction(psi,
                                                               energy_value)[0];
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

Result_one_point one_plot_point(double beta, Point * p_grid, double * WerOpt,
                                double * WerAc, double var_value,
                                const string & filename_base) {
    logger(LOG_INFO,
           "Calculate current density for " + to_string(var_value) + "...");

    time_t time_load = time(NULL);

    double * values_x =
        new double[config::model.particles]; // плотность постоянного
    // тока вдоль Ох (до
    // усреднения по ансамблю)
    double * values_y =
        new double[config::model.particles]; // плотность постоянного
    // тока вдоль Оу (до
    // усреднения по ансамблю)
    double * av_time =
        new double[config::model.particles]; // среднее время релаксации
    // (до усреднения по
    // ансамблю)
    unsigned int * nOpt =
        new unsigned int[config::model.particles]; // количество рассеяний на
    // оптических фононах
    unsigned int * nAc =
        new unsigned int[config::model.particles]; // количество рассеяний на
    // акустических фононах

    //Массивы-логи
    int num_logs = (int) (config::model.all_time / (config::model.dt)) + 1;
    double * px_log = new double[num_logs * (config::model.particles)];
    double * py_log = new double[num_logs * (config::model.particles)];

    // Инициализируем генератор случайных чисел
    srand(time(NULL));
    unsigned int * seed = new unsigned int[config::model.particles];
    for (int j = 0; j < config::model.particles; j++) {
        seed[j] = ((unsigned int) rand()) % 100000000 + 100000000;
    };

    omp_set_num_threads(config::model.threads);
#pragma omp parallel for
    for (int j = 0; j < config::model.particles; j++) {
        jobKernel(values_x, values_y, av_time, nAc, nOpt, beta, seed[j], j,
                  p_grid, WerAc, WerOpt, px_log, py_log, num_logs);
    };

    // Записываем значения компонент квазиимпульса в каждый момент времени для
    // каждой частицы
    // cout << "Writing logs for " << var_value << " .. ";
    // write_logs(filename_base, var_value, px_log, py_log, num_logs, *params);
    // cout << "done " << endl;

    // Помещаем результат расчета для одной точки графика в структуру
    Result_one_point result;
    result.result_value_mas_x = Mean(values_x, config::model.particles);
    result.result_value_mas_y = Mean(values_y, config::model.particles);
    result.std_values_mas_x = Std(values_x, config::model.particles);
    result.std_values_mas_y = Std(values_y, config::model.particles);
    result.result_av_time = Mean(av_time, config::model.particles);
    result.result_nOpt = Mean(nOpt, config::model.particles);
    result.result_nAc = Mean(nAc, config::model.particles);

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