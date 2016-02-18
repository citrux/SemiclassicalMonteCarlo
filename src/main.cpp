#include <fstream>
#include <iomanip>
#include "modelling_tools.h"
#include "logger.h"
#include "config.h"
#include "find_probability.h"
using namespace std;

int main(int argc, char const * argv[]) {
    // задаем параметры задачи
    string config = "config.ini";

    if (argc == 2)
        config = argv[1];

    logger(LOG_INFO, "> Loading config\n");

    load_config(config);

    set_probabilities();

    ofstream f;
    f.open(config::files.result);
    f << "#" << setw(15) << config::plot.var << setw(15) << "jx" << setw(15)
      << "djx" << setw(15) << "jy" << setw(15) << "djy" << setw(15) << "nopt"
      << setw(15) << "nac" << setw(15) << "tau" << endl;
    for (double value = config::plot.low; value <= config::plot.high;
         value += config::plot.step) {
        set_var(config::plot.var, value);
        Result res = result();
        f << scientific << " " << setw(15) << value << setw(15)
          << res.current_mean.x << setw(15) << res.current_std.x << setw(15)
          << res.current_mean.y << setw(15) << res.current_std.y << fixed
          << setw(15) << res.n_opt << setw(15) << res.n_ac << scientific
          << setw(15) << res.tau << endl;
    }

    f.close();
}