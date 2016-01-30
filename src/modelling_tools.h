#include "text_files.h"
#include <limits.h>
#include <stdio.h>

struct Result_one_point {
    double result_value_mas_x;
    double result_value_mas_y;
    double std_values_mas_x;
    double std_values_mas_y;
    double result_av_time;
    double result_nOpt;
    double result_nAc;
};

double random_uniform(unsigned int & x1, unsigned int & y1, unsigned int & z1,
                      unsigned int & w1);
void runge(double & px, double & py, double t, Params * params);
Point init_dist(unsigned int & x1, unsigned int & y1, unsigned int & z1,
                unsigned int & w1, Params * params);
double Mean(double * arr, int count);
double Mean(unsigned int * arr, int count);
double Std(double * arr, int count);

void jobKernel(double * dev_average_value_x, double * dev_average_value_y,
               double * dev_average_time_array, unsigned int * dev_nAc,
               unsigned int * dev_nOpt,

               Params * params, double beta, unsigned int rand_init_value,
               int idx,

               double * px_mas, double * py_mas, double * res_ac,
               double * res_opt, double * px_log, double * py_log, int num_logs
               /*int *mas_ac_opt, double *mas_t, int *mas_iCount, int n_step*/
               );

Result_one_point one_graphic_point(Params * params, double beta,
                                   double * px_mas, double * py_mas,
                                   double * WerOpt, double * WerAc,
                                   double var_value, char * filename_base);