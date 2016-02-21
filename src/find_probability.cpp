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
        while (fabs(output[i] - output_old) / output_old >
                   config::probability.probability_error &&
               count) {
            double dtheta = 2 * M_PI / n;
            output_old = output[i];
            output[i] = 0;

            for (double theta = 0; theta < 2 * M_PI; theta += dtheta) {

                curr = momentums_with_energy_in_direction(theta, e);

                for (size_t k = 0; k < min(prev.size(), curr.size()); ++k)
                    output[i] +=
                        len(prev[k] - curr[k]) / len(energy_gradient(curr[k]));

                swap(prev, curr);
            }
            n *= 2;
            --count;
        }
    }
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