#include <string>
#include "text_files.h"
#include "ini.h"
#include "logger.h"

using namespace std;

void results_to_file(const string & filename, int var_mas_count,
                     double * var_mas, double * result_value_mas_x,
                     double * result_value_mas_y, double * std_values_mas_x,
                     double * std_values_mas_y, double * result_av_time,
                     double * result_nOpt, double * result_nAc) {
    logger(LOG_INFO, "write results to '" + filename + "'...");

    ofstream f;
    f.open(filename);

    /*f << "Exc: " << params.Exc << "; Eyc: " << params.Eyc << "; wlo_max: " <<
    params.wlo_max << "; wla_max: " << params.wla_max << "; low: " << low << ";
    high: " << high << "; step: " << step << "; num_threads: " <<
    params.thread << endl;
    f << "Time info: " << total_time << endl;*/
    f << "#  E_y\tj_x\tj_y\tsigma_x\tsigma_y\ttau\tn_opt\tn_ac" << endl;
    for (int i = 0; i < var_mas_count; i++) {
        f << var_mas[i] << "\t" << result_value_mas_x[i] << "\t"
          << result_value_mas_y[i] << "\t" << std_values_mas_x[i] << "\t"
          << std_values_mas_y[i] << "\t" << result_av_time[i] << "\t"
          << result_nOpt[i] << "\t" << result_nAc[i] << endl;
    }
    f.close();

    logger(LOG_OK, "\t[DONE]\n");
}

double str_to_double(const string & str) {
    istringstream istream(str);
    double res;
    istream >> res;
    return res;
}

Point to_point(const string & str) {
    istringstream istream(str);
    Point p;
    istream >> p.x >> p.y;
    return p;
}

vec2 to_vec2(const string & str) {
    istringstream istream(str);
    vec2 v;
    istream >> v.x >> v.y;
    return v;
}


void load_config(const string & filename, Params & params) {

    logger(LOG_INFO, "load config from '" + filename + "'...\t");

    INIReader reader(filename);

    if (reader.ParseError() < 0) {
        logger(LOG_ERROR, "[FAIL]\n");
        logger(LOG_ERROR, "Can't load '" + filename + "'\n");
        exit(1);
    }

    params.files.load = reader.GetInteger("files", "load", 0);
    params.files.filename_opt = reader.Get("files", "optical", "UNKNOWN");
    params.files.filename_ac = reader.Get("files", "acoustical", "UNKNOWN");
    params.files.filename_res = reader.Get("files", "result", "UNKNOWN");

    params.phonons.beta = reader.GetReal("phonons", "beta", 0);
    params.phonons.wla_max = reader.GetReal("phonons", "wla_max", 0);
    params.phonons.wlo_max = reader.GetReal("phonons", "wlo_max", 0);
    params.phonons.T = reader.GetReal("phonons", "T", -1);

    params.fields.E0 = to_vec2(reader.Get("fields", "E0", "0 0"));
    params.fields.E1 = to_vec2(reader.Get("fields", "E1", "0 0"));
    params.fields.E2 = to_vec2(reader.Get("fields", "E2", "0 0"));
    params.fields.H = reader.GetReal("fields", "H",0);
    params.fields.omega1 = reader.GetReal("fields", "omega1",0);
    params.fields.omega2 = reader.GetReal("fields", "omega2",0);
    params.fields.phi1 = reader.GetReal("fields", "phi1",0);
    params.fields.phi2 = reader.GetReal("fields", "phi2",0);
    params.fields.phi = reader.GetReal("fields", "phi",0);

    params.bzone.A = to_point(reader.Get("bzone", "A", "0 0"));
    params.bzone.B = to_point(reader.Get("bzone", "B", "0 0"));
    params.bzone.D = to_point(reader.Get("bzone", "D", "0 0"));
    params.bzone.C = params.bzone.B + (params.bzone.D - params.bzone.A);

    params.probability.p_error = reader.GetReal("probability", "p_error",0);
    params.probability.p_points = reader.GetInteger("probability", "p_points",0);
    params.probability.n_integral = reader.GetInteger("probability", "n_integral",0);
    params.probability.e_points = reader.GetInteger("probability", "e_points",0);
    params.probability.optical = new double[params.probability.e_points];
    params.probability.acoustical = new double[params.probability.e_points];

    params.model.dt = reader.GetReal("model", "dt", 1e-4);
    params.model.all_time = reader.GetReal("model", "all_time", 0);
    params.model.threads = reader.GetInteger("model", "threads", 1);
    params.model.particles = reader.GetInteger("model", "particles", 1);

    params.plot.low = reader.GetReal("plotting", "low", 0);
    params.plot.high = reader.GetReal("plotting", "high", 0);
    params.plot.step = reader.GetReal("plotting", "step", 0);
    params.plot.num_var = reader.GetInteger("plotting", "num_var", -1);


    logger(LOG_OK, "[DONE]\n");
}

int get_var_mas_count(const Params & params) {
    double low = params.plot.low, high = params.plot.high, step = params.plot.step;
    int var_mas_count = (fabs(high - low) / step) + 1;
    return var_mas_count;
}

void set_var_mas(const Params & params, int var_mas_count, double * var_mas) {
    double low = params.plot.low, step = params.plot.step;
    for (int i = 0; i < var_mas_count; i++)
        var_mas[i] = low + step * i;
}

void var_value_plot(int num_param, double var_value, Params & params) {
    switch (num_param) {
    default:
        logger(LOG_ERROR, "Error!");
        exit(1);
    }
}