#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <time.h>
#include "find_probability.h"
using namespace std;
void array_from_file(const string & filename, double * px_mas, double * py_mas,
                     double * Wer, Params params);
void array_to_file(const string & filename, double * px_mas, double * py_mas,
                   double * Wer, Params params);
void write_logs(const string & filename_base, double var_value, double * px_log,
                double * py_log, int num_logs, Params params);
void results_to_file(const string & filename_res, double * var_mas, int var_mas_count,
                     double * result_value_mas_x, double * result_value_mas_y,
                     double * std_values_mas_x, double * std_values_mas_y,
                     double * result_av_time, double * result_nOpt,
                     double * result_nAc);
double str_to_double(string & str);
void params_from_file(const string & filename, Params * params, Graphic * graphic,
                      int * flag_count_or_load);
int get_var_mas_count(Graphic graphic);
void set_var_mas(Graphic graphic, int var_mas_count, double * var_mas);
void var_value_graphic(int num_param, double var_value, Params * params);
void comments_from_file(const string & filename_source, string * mas_comment,
                        int num_letters, int num_lines, int size);
void params_to_file(const string & filename_dest, string & filename_source,
                    Params * params, Graphic * graphic, int total_time,
                    int flag_count_or_load);
void filename_from_file(const string & filename_source, string & filename_params_after,
                        string & filename_opt, string & filename_ac,
                        string & filename_res);