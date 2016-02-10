#include <cstdlib>
#include <sstream>
#include "config.h"
#include "logger.h"
#include "ini.h"

using namespace std;

Files config::files;
Phonons config::phonons;
Probability config::probability;
Fields config::fields;
Bzone config::bzone;
Plot config::plot;
Model config::model;

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

void set_var(const string & varname, double value) {
    using namespace config;
    if (varname == "E0.x")
        fields.E0.x = value;
    if (varname == "E0.y")
        fields.E0.y = value;
    if (varname == "E1.x")
        fields.E1.x = value;
    if (varname == "E1.y")
        fields.E1.y = value;
    if (varname == "E2.x")
        fields.E2.x = value;
    if (varname == "E2.y")
        fields.E2.y = value;
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
    bzone.basis = new vec2[2];
    bzone.basis[0] = bzone.B - bzone.A;
    bzone.basis[1] = bzone.D - bzone.A;
    bzone.dual_basis = new vec2[2];
    // строим взаимный базис
    vec2 v =
        ort(bzone.basis[0] -
            dot(bzone.basis[0], ort(bzone.basis[1])) * ort(bzone.basis[1]));
    bzone.dual_basis[0] = v / dot(v, bzone.basis[0]);
    v = ort(bzone.basis[1] -
            dot(bzone.basis[1], ort(bzone.basis[0])) * ort(bzone.basis[0]));
    bzone.dual_basis[1] = v / dot(v, bzone.basis[1]);

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

    plot.low = reader.GetReal("plot", "low", 0);
    plot.high = reader.GetReal("plot", "high", 0);
    plot.step = reader.GetReal("plot", "step", 0);
    plot.var = reader.Get("plot", "var", "");

    logger(LOG_OK, "[DONE]\n");
}
