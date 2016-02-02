#include <fstream>
#include "find_probability.h"
#include "logger.h"

using namespace std;

vector<Point> momentums_with_energy_in_direction(double psi,
                                                 double energy_value,
                                                 const Params & params) {
    int n = params.probability.p_points;
    vec2 dir = {cos(psi), sin(psi)};
    vec2 step = pmax(psi, params) * dir / n;
    // logger(LOG_INFO, "pmax: " + to_string(pmax(psi, params)) + " \n");
    // logger(LOG_INFO, "step: " + to_string(step.x) + ", " + to_string(step.y)
    // + " \n");
    vector<Point> ps;
    Point O = {0, 0};
    for (int i = 0; i < n; ++i) {
        Point left = O + i * step, right = left + step;
        // logger(LOG_INFO, "momentum: left = " + to_string(left.x) + ", right =
        // " + to_string(right.x) + " \n");
        // logger(LOG_INFO, "energy: left = " + to_string(energy(left)) + ",
        // right = " + to_string(energy(right)) + " \n");
        if ((energy(left) - energy_value) * (energy(right) - energy_value) <
            0) {
            // logger(LOG_INFO, "different\n");
            while (len(right - left) > params.probability.p_error) {
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

void calculate_probability(double energy_shift, double * probability,
                           Params & params) {
    double e_min, e_max;
    get_energy_limits(e_min, e_max);
    double step = (e_max - e_min) / params.probability.e_points;
    double dpsi = 2 * M_PI / params.probability.n_integral;

    // logger(LOG_INFO, "dpsi = " + to_string(dpsi) + "\n");

    for (int i = 0; i < params.probability.e_points; ++i) {
        double e = e_min + i * step;
        params.probability.energy[i] = e;
        vector<Point> prev, curr;
        for (int j = 0; j < params.probability.n_integral; ++j) {
            curr = momentums_with_energy_in_direction(j * dpsi,
                                                      e - energy_shift, params);
            if (prev.size())
                for (int k = 0; k < min(prev.size(), curr.size()); ++k)
                    probability[i] +=
                        len(curr[k] - prev[k]) /
                        len(velocity(curr[k] + (prev[k] - curr[k]) / 2));
            swap(prev, curr);
        }
    }
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

void set_probabilities(Params & params) {
    if (params.files.load) {
        load_arrays(params.files.filename_ac,
                    {params.probability.energy, params.probability.acoustical},
                    params.probability.e_points);
        load_arrays(params.files.filename_opt,
                    {params.probability.energy, params.probability.optical},
                    params.probability.e_points);
    } else {
        calculate_probability(0, params.probability.acoustical, params);
        calculate_probability(params.phonons.beta, params.probability.optical,
                              params);
        save_arrays(params.files.filename_ac,
                    {params.probability.energy, params.probability.acoustical},
                    params.probability.e_points);
        save_arrays(params.files.filename_opt,
                    {params.probability.energy, params.probability.optical},
                    params.probability.e_points);
    }
}

double get_probability(Point p, const Params & params) {
    double e_min, e_max;
    get_energy_limits(e_min, e_max);
    double step = (e_max - e_min) / params.probability.e_points;

    double e = energy(p);
    double pos = (e - e_min) / step;
    int i = floor(pos);
    double w = pos - i;

    return (1 - w) *
               (params.phonons.wla_max * params.probability.acoustical[i] +
                params.phonons.wlo_max * params.probability.optical[i]) +
           w * (params.phonons.wla_max * params.probability.acoustical[i + 1] +
                params.phonons.wlo_max * params.probability.optical[i + 1]);
}

double distrib_function(double p, double psi, const Params & params) {
    return exp((-1.0 / params.phonons.T) * energy_psi(p, psi)) * p;
}