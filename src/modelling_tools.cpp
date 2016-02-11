#include <ctime>
#include <fstream>
#include "modelling_tools.h"
#include "logger.h"
#include "material_specific.h"
#include "find_probability.h"

double distrib_function(double p, double theta) {
    return exp(-energy_theta(p, theta) / config::phonons.T) * p;
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
        angle_int_distrib = new double[n + 1];
        angle_int_distrib[0] = 0;
        for (int i = 1; i <= n; ++i) {
            double theta = 2 * M_PI * i / n;
            double momentum_int = 0;
            double step = pmax(theta) / m;
            for (int j = 0; j < m; ++j)
                momentum_int +=
                    distrib_function((j + 0.5) * step, theta) * step;

            angle_int_distrib[i] = angle_int_distrib[i - 1] + momentum_int;
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

    double w = (prob - angle_int_distrib[i]) /
               (angle_int_distrib[j] - angle_int_distrib[i]);
    return (i + w) * 2 * M_PI / n;
}

double inverse_function_momentum(double theta, double prob) {
    int n = 1000;
    double * int_distrib = new double[n + 1];
    double p = pmax(theta), step = p / n;

    int_distrib[0] = 0;
    // табулируем интегральную функцию
    for (int i = 1; i <= n; ++i) {
        int_distrib[i] = int_distrib[i - 1] +
                         distrib_function((i - 0.5) * step, theta) * step;
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
    delete[] int_distrib;
    return (i + w) * step;
}

Point init_dist(unsigned int & x1, unsigned int & y1, unsigned int & z1,
                unsigned int & w1) {
    double r = random_uniform(x1, y1, z1, w1);
    double theta = inverse_function_angle(r);

    r = random_uniform(x1, y1, z1, w1);
    double p = inverse_function_momentum(theta, r);

    return {p * cos(theta), p * sin(theta)};
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

double mean(double * arr, int count) {
    double result = 0;

    for (int i = 0; i < count; i++) {
        result += arr[i];
    }

    return result / (double) count;
}

double mean(unsigned int * arr, int count) {
    double result = 0;

    for (int i = 0; i < count; i++) {
        result += arr[i];
    }

    return result / (double) count;
}

vec2 mean(vec2 * arr, int count) {
    vec2 result;

    for (int i = 0; i < count; i++) {
        result += arr[i];
    }

    return result / (double) count;
}

double sd(double * arr, int count) {
    double sum = 0.0;
    double m = mean(arr, count);

    for (int i = 0; i < count; i++) {
        sum += pow((arr[i] - m), 2);
    }

    return sqrt(sum / (count - 1));
}

vec2 sd(vec2 * arr, int count) {
    double sum1 = 0.0, sum2 = 0.0;
    vec2 m = mean(arr, count);

    for (int i = 0; i < count; i++) {
        sum1 += pow((arr[i].x - m.x), 2);
        sum2 += pow((arr[i].y - m.y), 2);
    }

    return {sqrt(sum1 / (count - 1)), sqrt(sum2 / (count - 1))};
}

void job_kernel(const Point & init_condition, unsigned int seed, vec2 & current,
               double & tau, unsigned int & n_ac, unsigned int & n_opt) {
    // logger(LOG_INFO, "Start job_kernel\n");

    int count;

    vec2 v;

    unsigned int x1 = seed, y1 = 362436069, z1 = 521288629,
                 w1 = 88675123; // RNGInit

    Point p = init_condition;

    double t = 0.0;
    double wla = 0.0;
    double wlo = 0.0;

    int n0 = 0;
    n_ac = 0, n_opt = 0;

    double beta_opt = config::phonons.beta;

    int t_num = 0;
    while (t < config::model.all_time) {
        v = velocity(p);

        current += v * config::model.dt;

        runge(p, t); // решаем уравнения движения

        /* приводим импульс к зоне */
        p = to_first_bz(p);

        t += config::model.dt;


        double e = energy(p);
        wlo += config::phonons.wlo_max * get_probability(e - beta_opt) *
               config::model.dt;
        wla += config::phonons.wla_max * get_probability(e) * config::model.dt;

        if (wla > random_uniform(x1, y1, z1, w1)) {
            ++n_ac; // наращиваем счетчик рассеяний на акустических фононах
            wla = 0.0;
            e = energy(p);
            count = 15;
            while (count) {
                double theta =
                    2 * M_PI *
                    random_uniform(x1, y1, z1, w1); // случайным образом
                // разыгрываем направление квазиимпульса
                auto ps = momentums_with_energy_in_direction(theta, e);
                if (ps.size()) {
                    p = ps[0];
                    break;
                }
                // если p существует, то мы правильно
                // подобрали угол рассеяния, поэтому выходим из цикла
                // если за 15 попыток не нашли решение, выходим из цикла
                --count;
            }
        } else if ((wlo > random_uniform(x1, y1, z1, w1)) &&
                   (energy(p) >= beta_opt)) {
            ++n_opt; // наращиваем счетчик рассеяний на оптических
                     // фононах
            wlo = 0;
            e = energy(p) - beta_opt;
            count = 15;

            while (count) {
                double theta = 2 * M_PI * random_uniform(x1, y1, z1,
                                                         w1); // случайным
                                                              // образом
                // разыгрываем направление
                auto ps = momentums_with_energy_in_direction(theta, e);
                if (ps.size()) {
                    p = ps[0];
                    break;
                }
                --count;
            }
        }
        t_num++;
    }
    n0 = n_ac + n_opt;
    current = current / t;
    tau = t / (n0 + 1);
    // logger(LOG_INFO, "End job_kernel\n");
}

Result result() {

    time_t time_load = time(nullptr);

    vec2 * current =
        new vec2[config::model.particles]; // плотность постоянного тока
    // (до усреднения по ансамблю)
    double * tau = new double[config::model.particles];
    // среднее время релаксации (до усреднения по ансамблю)
    unsigned int * n_opt =
        new unsigned int[config::model.particles]; // количество рассеяний на
    // оптических фононах
    unsigned int * n_ac =
        new unsigned int[config::model.particles]; // количество рассеяний на
    // акустических фононах

    // Инициализируем генератор случайных чисел
    srand(time(nullptr));
    Point * init_condition = new Point[config::model.particles];
    unsigned int * seed = new unsigned int[config::model.particles];

    unsigned int x1 = rand();
    unsigned int y1 = 362436069;
    unsigned int z1 = 521288629;
    unsigned int w1 = 88675123; // RNGInit

    for (int j = 0; j < config::model.particles; j++) {
        init_condition[j] = init_dist(x1, y1, z1, w1);
        seed[j] = rand();
    };

    omp_set_num_threads(config::model.threads);
#pragma omp parallel for
    for (int j = 0; j < config::model.particles; j++) {
        job_kernel(init_condition[j], seed[j], current[j], tau[j], n_ac[j],
                  n_opt[j]);
    };

    // Помещаем результат расчета для одной точки графика в структуру
    Result result;
    result.current_mean = mean(current, config::model.particles);
    result.current_std = sd(current, config::model.particles);
    result.tau = mean(tau, config::model.particles);
    result.n_opt = mean(n_opt, config::model.particles);
    result.n_ac = mean(n_ac, config::model.particles);

    // Удаляем массивы, хранящие значения скорости, среднего времени релаксации
    // и т.д. до усреднения по ансамблю
    delete[] current;
    delete[] tau;
    delete[] n_opt;
    delete[] n_ac;
    delete[] seed;
    delete[] init_condition;

    logger(LOG_OK,
           "\t[DONE in " + to_string(time(nullptr) - time_load) + " sec]\n");
    return result;
}