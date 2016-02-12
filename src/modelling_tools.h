#include <climits>
#include <cstdio>
#include "text_files.h"
#include "linalg.h"
#include "config.h"

struct Result {
    vec2 current_mean;
    vec2 current_std;
    double tau;
    double n_opt;
    double n_ac;
};

double random_uniform(unsigned int & x1, unsigned int & y1, unsigned int & z1,
                      unsigned int & w1);
void runge(Point & p, double t);
Point init_dist(unsigned int & x1, unsigned int & y1, unsigned int & z1,
                unsigned int & w1);
double mean(double * arr, int count);
double mean(unsigned int * arr, int count);
double sd(double * arr, int count);

void job_kernel(const Point & init_condition, unsigned int seed, vec2 & current,
                double & tau, unsigned int & n_ac, unsigned int & n_opt);

Result result();