#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <time.h>
#include "find_probability.h"
using namespace std;

void results_to_file(const string & filename_res, int var_mas_count,
                     double * var_mas, double * result_value_mas_x,
                     double * result_value_mas_y, double * std_values_mas_x,
                     double * std_values_mas_y, double * result_av_time,
                     double * result_nOpt, double * result_nAc);
double str_to_double(const string & str);
void load_config(const string & filename);
int get_var_mas_count();
void set_var_mas(int var_mas_count, double * var_mas);
void var_value_plot(int num_param, double var_value);