#include <vector>
#include <omp.h>
#include "find_probability.h"

using namespace std;

vector<Point> momentums_with_energy_in_direction(double psi, double energy_value,
                                          const Params & params) {
    int n = params.probability.p_points;
    vec2 dir = {cos(psi), sin(psi)};
    vec2 step = pmax(psi, params) * dir / n;
    vector<Point> ps;
    Point left = {0, 0}, right = left + step;
    for (int i = 0; i < n; ++i)
    {
        if ((energy(left) - energy_value) *
            (energy(right) - energy_value) < 0) {

            while(len(right - left) > params.probability.p_error) {
                Point middle = left + (right - left) / 2;
                if ((energy(left) - energy_value) * (energy(middle) - energy_value) < 0)
                    right = middle;
                else
                    left = middle;
            }

            ps.push_back(left + (right - left) / 2);
        }
        swap(left, right);
        right = left + step;
    }
    return ps;
}

void get_energy_limits(double & emin, double & emax) {
    emin = 0;
    emax = 1;
}

void calculate_probability(double energy_shift, double * probability, Params & params) {
    double e_min, e_max;
    get_energy_limits(e_min, e_max);
    double step = (e_max - e_min) / params.probability.e_points;
    double dpsi = 2 * M_PI / params.probability.n_integral;

    for (int i = 0; i < params.probability.e_points; ++i) {
        double e = e_min + i * step;
        vector<Point> prev, curr;
        for (int j = 0; j < params.probability.n_integral; ++j) {
            curr = momentums_with_energy_in_direction(j * dpsi, e - energy_shift, params);
            if (prev.size())
                for (int k = 0; k < min(prev.size(), curr.size()); ++k)
                    probability[i] += len(curr[k] - prev[k]) / len(velocity(curr[k] + (prev[k] - curr[k])/2));
            swap(prev, curr);
        }
    }
}

void set_probabilities(Params & params) {
    calculate_probability(0, params.probability.acoustical, params);
    calculate_probability(params.phonons.beta, params.probability.optical, params);
}


double get_probability(Point p, Params & params) {
    double e_min, e_max;
    get_energy_limits(e_min, e_max);
    double step = (e_max - e_min) / params.probability.e_points;

    double e = energy(p);
    double pos = (e - e_min) / step;
    int i = floor(pos);
    double w = pos - i;

    return (1 - w) * (params.probability.acoustical[i] + params.probability.optical[i]) +
               w * (params.probability.acoustical[i+1] + params.probability.optical[i+1]);
}
