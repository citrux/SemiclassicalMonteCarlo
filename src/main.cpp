#include <fstream>
#include "modelling_tools.h"
#include "logger.h"
#include "config.h"
#include "find_probability.h"
using namespace std;

void plot(const string & fname) {
    ofstream f;
    f.open(fname);
    f << "# j(" << config::plot.var << "), " << config::plot.var
      << " = np.arange(" << config::plot.low << ", " << config::plot.high
      << ", " << config::plot.step << ")" << endl;
    f << "#" << setw(15) << "E0.x" << setw(15) << "E0.y" << setw(15) << "E1.x"
      << setw(15) << "E1.y" << setw(15) << "E2.x" << setw(15) << "E2.y"
      << setw(15) << "B0" << setw(15) << "B1" << setw(15) << "B2" << setw(15)
      << "omega1" << setw(15) << "omega2" << setw(15) << "phi" << setw(15)
      << "jx" << setw(15) << "djx" << setw(15) << "jy" << setw(15) << "djy"
      << setw(15) << "nopt" << setw(15) << "nac" << setw(15) << "tau" << endl;
    for (double value = config::plot.low; value <= config::plot.high;
         value += config::plot.step) {
        set_var(config::plot.var, value);
        Result res = result();
        f << " " << config::fields << res << endl;
    }

    f.close();
}

void test_zero(const string & fname, int n) {
    using namespace config;
    ofstream f;
    f.open(fname);

    for (int i = 0; i < n; ++i) {
        double dt_new = 1e-2, time_new = 1e-2;
        swap(dt_new, model.dt);
        swap(time_new, model.all_time);
        Result res = result();

        f << res.current.x << "\t" << res.current.y << "\t" << res.current_std.x
          << "\t" << res.current_std.y << "\t";

        swap(dt_new, model.dt);
        swap(time_new, model.all_time);
        res = result();

        f << res.current.x << "\t" << res.current.y << "\t" << res.current_std.x
          << "\t" << res.current_std.y << endl;
    }
    f.close();
}

int main(int argc, char const * argv[]) {
    // задаем параметры задачи
    string config = "config.ini";

    if (argc == 2)
        config = argv[1];

    logger(LOG_INFO, "> Loading config\n");

    load_config(config);

    set_probabilities();

    // plot(config::files.result);

    test_zero("data/test_zero.dat", 100);
}
