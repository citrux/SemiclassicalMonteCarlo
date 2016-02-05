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

Point point_on_contour(Point a, double nrg, const Params & params) {
    vec2 dir = ort(velocity(a));
    Point b = a - dir * (energy(a) - nrg) / dot(velocity(a), dir);
    while (len(a-b) > params.probability.p_error) {
        a = b;
        b = a - dir * (energy(a) - nrg) / dot(velocity(a), dir);
    }

    return b;
}

double integrate(Point a, Point b, double nrg, const Params & params) {
    // double step = 0.00001;
    // double result = 0;

    // // logger(LOG_INFO, "points: \n");

    // Point c = a;
    // int i = 0;
    // while (len(c - b) > step && i < 20) {
    //     // logger(LOG_INFO, "> " + to_string(c.x) + ", " + to_string(c.y) + "; speed: " + to_string(len(velocity(c))) + " \n");
    //     vec2 dir = ort(b - c);
    //     Point d = point_on_contour(c + step * dir, nrg, params);
    //     result += len(d-c)/len(velocity(c));
    //     c = d;
    //     ++i;
    // }
    // result += len(c-b)/len(velocity(c));
    // // logger(LOG_INFO, "result: " + to_string(result) + "\n");
    // return result;
    Point c = point_on_contour(a + (b-a)/2, nrg, params);
    double l1 = len(c-a), l2 = len(b-c), l = l1 + l2;
    double wa = (1./3 - l2 / l1 / 6) * l,
           wc = l * l * l / 6 / l1 / l2,
           wb = (1./3 - l1 / l2 / 6) * l;
    return wa / len(velocity(a)) + wb / len(velocity(b)) + wc / len(velocity(c));
}

void calculate_probability(double * probability, Params & params) {
    double e_min, e_max;
    get_energy_limits(e_min, e_max);
    double step = (e_max - e_min) / params.probability.e_points;

    // logger(LOG_INFO, "dpsi = " + to_string(dpsi) + "\n");
    omp_set_num_threads(params.model.threads);
    #pragma omp parallel for
    for (int i = 0; i < params.probability.e_points; ++i) {
        double e = e_min + i * step;
        params.probability.energy[i] = e;
        vector<Point> prev, curr;
        double probability_old = 0;
        probability[i] = 1;
        int n = params.probability.n_integral;
        int count = 12;
        while (fabs(probability[i] - probability_old) / probability_old > 1e-5 && count) {
            double dpsi = 2 * M_PI / n;
            probability_old = probability[i];
            probability[i] = 0;
            // logger(LOG_INFO, to_string(e) + " " +to_string(n) + " " + to_string(probability_old) + "\n");
            for (double psi = 0; psi < 2 * M_PI; psi += dpsi) {
                // logger(LOG_INFO, to_string(1.0 * j / params.probability.n_integral) + "\n");
                curr = momentums_with_energy_in_direction(psi, e, params);

                if (psi > 0) {
                    for (int k = curr.size(); k < prev.size(); ++k) {
                        vec2 dir = {cos(psi-dpsi), sin(psi - dpsi)};
                        Point end = pmax(psi - dpsi, params) * dir;
                        probability[i] += integrate(prev[k], end, e, params);
                    }

                    for (int k = prev.size(); k < curr.size(); ++k) {
                        vec2 dir = {cos(psi), sin(psi)};
                        Point start = pmax(psi, params) * dir;
                        probability[i] += integrate(start, curr[k], e, params);
                    }

                    if (prev.size() == curr.size()) {
                        for (int k = 0; k < curr.size(); ++k)
                            probability[i] += integrate(prev[k], curr[k], e, params);
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

void set_probabilities(Params & params) {
    if (params.files.load) {
        load_arrays(params.files.probability,
                    {params.probability.energy, params.probability.probability},
                    params.probability.e_points);
    } else {
        calculate_probability(params.probability.probability, params);
        save_arrays(params.files.probability,
                    {params.probability.energy, params.probability.probability},
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

    return (1 - w) * params.probability.probability[i] +
           w * params.probability.probability[i + 1];
}

double distrib_function(double p, double psi, const Params & params) {
    return exp((-1.0 / params.phonons.T) * energy_psi(p, psi)) * p;
}