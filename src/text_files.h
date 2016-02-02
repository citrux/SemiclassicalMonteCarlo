#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <time.h>
#include "find_probability.h"
using namespace std;
void array_from_file(const string & filename, Point * p_grid, double * Wer,
                     const Params & params);
void array_to_file(const string & filename, Point * p_grid, double * Wer,
                   const Params & params);
void write_logs(const string & filename_base, double var_value, double * px_log,
                double * py_log, int num_logs, const Params & params);
void results_to_file(const string & filename_res, int var_mas_count,
                     double * var_mas, double * result_value_mas_x,
                     double * result_value_mas_y, double * std_values_mas_x,
                     double * std_values_mas_y, double * result_av_time,
                     double * result_nOpt, double * result_nAc);
double str_to_double(const string & str);
void load_config(const string & filename, Params & params, Graphic & graphic,
                 Probability & prob);
int get_var_mas_count(const Graphic & graphic);
void set_var_mas(const Graphic & graphic, int var_mas_count, double * var_mas);
void var_value_graphic(int num_param, double var_value, Params & params);