#include <fstream>
#include <sstream>
#include "text_files.h"
#include "logger.h"
#include "config.h"

using namespace std;

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

void results_to_file(const string & filename, int var_mas_count,
                     double * var_mas, double * result_value_mas_x,
                     double * result_value_mas_y, double * std_values_mas_x,
                     double * std_values_mas_y, double * result_av_time,
                     double * result_nOpt, double * result_nAc) {
    logger(LOG_INFO, "write results to '" + filename + "'...");

    save_arrays(filename, {var_mas, result_value_mas_x, result_value_mas_y, std_values_mas_x, std_values_mas_y, result_av_time, result_nOpt, result_nAc}, var_mas_count);

    logger(LOG_OK, "\t[DONE]\n");
}

double str_to_double(const string & str) {
    istringstream istream(str);
    double res;
    istream >> res;
    return res;
}


int get_var_mas_count() {
    double low = config::plot.low, high = config::plot.high,
           step = config::plot.step;
    int var_mas_count = (fabs(high - low) / step) + 1;
    return var_mas_count;
}

void set_var_mas(int var_mas_count, double * var_mas) {
    double low = config::plot.low, step = config::plot.step;
    for (int i = 0; i < var_mas_count; i++)
        var_mas[i] = low + step * i;
}

void var_value_plot(int num_param, double var_value) {
    switch (num_param) {
    default:
        logger(LOG_ERROR, "Error!");
        exit(1);
    }
}