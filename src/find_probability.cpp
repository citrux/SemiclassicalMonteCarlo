#include <fstream>
#include "find_probability.h"
#include "logger.h"

using namespace std;

vector<Point> momentums_with_energy_in_direction(double psi,
                                                 double energy_value) {
    int n = config::probability.momentum_samples;
    vec2 dir = {cos(psi), sin(psi)};
    vec2 step = pmax(psi) * dir / n;
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
    emin = 0;
    emax = 2;
}

Point point_on_contour(Point a, double nrg) {
    vec2 dir = ort(velocity(a));
    Point b = a - dir * (energy(a) - nrg) / dot(velocity(a), dir);
    while (len(a - b) > config::probability.momentum_error) {
        a = b;
        b = a - dir * (energy(a) - nrg) / dot(velocity(a), dir);
    }

    return b;
}

double integrate(Point a, Point b, double nrg) {
    Point c = point_on_contour(a + (b - a) / 2, nrg);
    double l1 = len(c - a), l2 = len(b - c), l = l1 + l2;
    double wa = (1. / 3 - l2 / l1 / 6) * l, wc = l * l * l / 6 / l1 / l2,
           wb = (1. / 3 - l1 / l2 / 6) * l;
    return wa / len(velocity(a)) + wb / len(velocity(b)) +
           wc / len(velocity(c));
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
            double dpsi = 2 * M_PI / n;
            output_old = output[i];
            output[i] = 0;
            // logger(LOG_INFO, to_string(e) + " " +to_string(n) + " " +
            // to_string(output_old) + "\n");
            for (double psi = 0; psi < 2 * M_PI; psi += dpsi) {
                // logger(LOG_INFO, to_string(1.0 * j /
                // config::output.n_integral) + "\n");
                curr = momentums_with_energy_in_direction(psi, e);

                if (psi > 0) {
                    for (size_t k = curr.size(); k < prev.size(); ++k) {
                        vec2 dir = {cos(psi - dpsi), sin(psi - dpsi)};
                        Point end = pmax(psi - dpsi) * dir;
                        output[i] += integrate(prev[k], end, e);
                    }

                    for (size_t k = prev.size(); k < curr.size(); ++k) {
                        vec2 dir = {cos(psi), sin(psi)};
                        Point start = pmax(psi) * dir;
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

void load_arrays(const string & filename, vector<double *> arrays, int n) {
    ifstream f;
    f.open(filename);
    for (int i = 0; i < n; ++i)
        for (double * a : arrays)
            f >> a[i];
    f.close();
}

void save_arrays(const string & filename, vector<double *> arrays, int n) {
    ofstream f;
    f.open(filename);
    for (int i = 0; i < n; ++i) {
        for (double * a : arrays)
            f << a[i] << "\t";
        f << "\n";
    }
    f.close();
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
    while (j - i > 1) {
        int k = (i + j) / 2;
        if (probability.energy[k] < energy)
            i = k;
        else
            j = k;
    }
    double w = (energy - probability.energy[i]) /
               (probability.energy[j] - probability.energy[i]);

    return (1 - w) * config::probability.probability[i] +
           w * config::probability.probability[j];
}

double distrib_function(double p, double psi) {
    return exp(-energy_psi(p, psi) / config::phonons.T) * p;
}