#pragma once
#include <string>
#include <vector>

void results_to_file(const std::string & filename_res, int var_mas_count,
                     double * var_mas, double * result_value_mas_x,
                     double * result_value_mas_y, double * std_values_mas_x,
                     double * std_values_mas_y, double * result_av_time,
                     double * result_nOpt, double * result_nAc);

void load_arrays(const std::string & filename, std::vector<double *> arrays, int n);
void save_arrays(const std::string & filename, std::vector<double *> arrays, int n);

double str_to_double(const std::string & str);
int get_var_mas_count();
void set_var_mas(int var_mas_count, double * var_mas);
void var_value_plot(int num_param, double var_value);