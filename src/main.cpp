#include <fstream>
#include "modelling_tools.h"
#include "logger.h"
#include "config.h"
#include "find_probability.h"
using namespace std;

int main() {
    // задаем параметры задачи
    string config = "config.ini";

    logger(LOG_INFO, "> Loading config\n");

    load_config(config);

    set_probabilities();

    ofstream f;
    f.open(config::files.result);
    f << "# " << config::plot.var
      << "\tjx\tjy\tsigma(jx)\tsigma(jy)\tnopt\tnac\ttau\n";
    for (double value = config::plot.low; value <= config::plot.high;
         value += config::plot.step) {
        set_var(config::plot.var, value);
        Result res = result();
        f << value << "\t" << res.current_mean.x << "\t" << res.current_mean.y
          << "\t" << res.current_std.x << "\t" << res.current_std.y << "\t"
          << res.n_opt << "\t" << res.n_ac << "\t" << res.tau << endl;
    }

    f.close();
}