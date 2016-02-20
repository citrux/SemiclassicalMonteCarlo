#include <fstream>
#include "find_probability.h"
#include "material_specific.h"
#include "logger.h"
#include "config.h"
#include "text_files.h"

using namespace std;

vector<Point> momentums_with_energy_in_direction(double theta,
                                                 double energy_value) {
    int n = config::probability.momentum_samples;
    vec2 dir = {cos(theta), sin(theta)};
    vec2 step = pmax(theta) * dir / n;
    vector<Point> ps;
    Point O = {0, 0};

    for (int i = 0; i < n; ++i) {
        Point left = O + i * step, right = left + step;
        if ((energy(left) - energy_value) * (energy(right) - energy_value) <
            0) {

            while (len(right - left) > config::probability.momentum_error) {
                Point middle = left + (right - left) / 2;
                if ((energy(left) - energy_value) *
                        (energy(middle) - energy_value) <
                    0)
                    right = middle;
                else
                    left = middle;
            }

            ps.push_back(left + (right - left) / 2);
        }
    }
    return ps;
}

void get_energy_limits(double & emin, double & emax) {
    using namespace config;
    emin = 100;
    emax = -100;
    int n = 1000;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double e = energy(bzone.A + bzone.basis[0] * i / n +
                              bzone.basis[1] * j / n);
            if (e < emin)
                emin = e;
            if (e > emax)
                emax = e;
        }
    }
}

Point point_on_contour(Point a, double nrg) {
    vec2 dir = ort(energy_gradient(a));
    Point b = a - dir * (energy(a) - nrg) / dot(energy_gradient(a), dir);
    while (len(a - b) > config::probability.momentum_error) {
        a = b;
        b = a - dir * (energy(a) - nrg) / dot(energy_gradient(a), dir);
    }

    return b;
}

double integrate(Point a, Point b, double nrg) {
    Point c = point_on_contour(a + (b - a) / 2, nrg);
    double l1 = len(c - a), l2 = len(b - c), l = l1 + l2;
    double wa = (1. / 3 - l2 / l1 / 6) * l, wc = l * l * l / 6 / l1 / l2,
           wb = (1. / 3 - l1 / l2 / 6) * l;
    return wa / len(energy_gradient(a)) + wb / len(energy_gradient(b)) +
           wc / len(energy_gradient(c));
}

void calculate_probability(double * output) {
    double e_min, e_max;
    get_energy_limits(e_min, e_max);
    double step = (e_max - e_min) / config::probability.energy_samples;

    omp_set_num_threads(config::model.threads);
#pragma omp parallel for
    for (int i = 0; i < config::probability.energy_samples; ++i) {
        double e = e_min + i * step;
        config::probability.energy[i] = e;
        vector<Point> prev, curr;
        double output_old = 0;
        output[i] = 1;
        int n = 500;
        int count = 12;
        while (fabs(output[i] - output_old) / output_old > 1e-5 && count) {
            double dtheta = 2 * M_PI / n;
            output_old = output[i];
            output[i] = 0;
            // logger(LOG_INFO, to_string(e) + " " +to_string(n) + " " +
            // to_string(output_old) + "\n");
            for (double theta = 0; theta < 2 * M_PI; theta += dtheta) {
                // logger(LOG_INFO, to_string(1.0 * j /
                // config::output.n_integral) + "\n");
                curr = momentums_with_energy_in_direction(theta, e);

                if (theta > 0) {
                    for (size_t k = curr.size(); k < prev.size(); ++k) {
                        vec2 dir = {cos(theta - dtheta), sin(theta - dtheta)};
                        Point end = pmax(theta - dtheta) * dir;
                        output[i] += integrate(prev[k], end, e);
                    }

                    for (size_t k = prev.size(); k < curr.size(); ++k) {
                        vec2 dir = {cos(theta), sin(theta)};
                        Point start = pmax(theta) * dir;
                        output[i] += integrate(start, curr[k], e);
                    }

                    if (prev.size() == curr.size()) {
                        for (size_t k = 0; k < curr.size(); ++k)
                            output[i] += integrate(prev[k], curr[k], e);
                    }
                }

                swap(prev, curr);
            }
            n *= 2;
            --count;
        }
        // logger(LOG_INFO, "epoint" +to_string(i) +"complete\n");
    }
    // logger(LOG_INFO, "calc finished\n");
}

void set_probabilities() {
    if (config::files.load) {
        load_arrays(
            config::files.probability,
            {config::probability.energy, config::probability.probability},
            config::probability.energy_samples);
    } else {
        calculate_probability(config::probability.probability);
        save_arrays(
            config::files.probability,
            {config::probability.energy, config::probability.probability},
            config::probability.energy_samples);
    }
}

double get_probability(double energy) {
    using config::probability;
    int i = 0, j = probability.energy_samples - 1;

    if (energy < probability.energy[i] || energy > probability.energy[j])
        return 0;

    // бисекция на произвольной сетке
    // while (j - i > 1) {
    //     int k = (i + j) / 2;
    //     if (probability.energy[k] < energy)
    //         i = k;
    //     else
    //         j = k;
    // }
    // double w = (energy - probability.energy[i]) /
    //            (probability.energy[j] - probability.energy[i]);

    // оптимизация для равномерной сетки
    double step = probability.energy[1] - probability.energy[0];
    double c = (energy - probability.energy[0]) / step;
    i = floor(c);
    j = i + 1;
    double w = c - i;

    return (1 - w) * config::probability.probability[i] +
           w * config::probability.probability[j];
}