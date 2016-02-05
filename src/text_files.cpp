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

    /*f << "Exc: " << config::Exc << "; Eyc: " << config::Eyc << "; wlo_max: "
    <<
    config::wlo_max << "; wla_max: " << config::wla_max << "; low: " << low <<
    ";
    high: " << high << "; step: " << step << "; num_threads: " <<
    config::model.threads << endl;
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

void load_config(const string & filename) {

    using namespace config;

    logger(LOG_INFO, "load config from '" + filename + "'...\t");

    INIReader reader(filename);

    if (reader.ParseError() < 0) {
        logger(LOG_ERROR, "[FAIL]\n");
        logger(LOG_ERROR, "Can't load '" + filename + "'\n");
        exit(1);
    }

    files.load = reader.GetInteger("files", "load", 0);
    files.probability = reader.Get("files", "probability", "UNKNOWN");
    files.result = reader.Get("files", "result", "UNKNOWN");

    phonons.beta = reader.GetReal("phonons", "beta", 0);
    phonons.wla_max = reader.GetReal("phonons", "wla_max", 0);
    phonons.wlo_max = reader.GetReal("phonons", "wlo_max", 0);
    phonons.T = reader.GetReal("phonons", "T", -1);

    fields.E0 = to_vec2(reader.Get("fields", "E0", "0 0"));
    fields.E1 = to_vec2(reader.Get("fields", "E1", "0 0"));
    fields.E2 = to_vec2(reader.Get("fields", "E2", "0 0"));
    fields.H = reader.GetReal("fields", "H", 0);
    fields.omega1 = reader.GetReal("fields", "omega1", 0);
    fields.omega2 = reader.GetReal("fields", "omega2", 0);
    fields.phi1 = reader.GetReal("fields", "phi1", 0);
    fields.phi2 = reader.GetReal("fields", "phi2", 0);
    fields.phi = reader.GetReal("fields", "phi", 0);

    bzone.A = to_point(reader.Get("bzone", "A", "0 0"));
    bzone.B = to_point(reader.Get("bzone", "B", "0 0"));
    bzone.D = to_point(reader.Get("bzone", "D", "0 0"));
    bzone.C = bzone.B + (bzone.D - bzone.A);

    probability.momentum_samples =
        reader.GetReal("probability", "momentum_samples", 0);
    probability.momentum_error =
        reader.GetReal("probability", "momentum_error", 0);
    probability.probability_error =
        reader.GetInteger("probability", "probability_error", 0);
    probability.energy_samples =
        reader.GetInteger("probability", "energy_samples", 0);
    probability.probability = new double[probability.energy_samples];
    probability.energy = new double[probability.energy_samples];

    model.dt = reader.GetReal("model", "dt", 1e-4);
    model.all_time = reader.GetReal("model", "all_time", 0);
    model.threads = reader.GetInteger("model", "threads", 1);
    model.particles = reader.GetInteger("model", "particles", 1);

    plot.low = reader.GetReal("plotting", "low", 0);
    plot.high = reader.GetReal("plotting", "high", 0);
    plot.step = reader.GetReal("plotting", "step", 0);
    plot.num_var = reader.GetInteger("plotting", "num_var", -1);

    logger(LOG_OK, "[DONE]\n");
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