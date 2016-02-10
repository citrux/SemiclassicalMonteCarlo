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

void jobKernel(double * dev_average_value_x, double * dev_average_value_y,
               double * dev_average_time_array, unsigned int * dev_nAc,
               unsigned int * dev_nOpt,

               double beta, unsigned int rand_init_value, int idx,

               Point * p_grid, double * res_ac, double * res_opt,
               double * px_log, double * py_log, int num_logs);

Result result();