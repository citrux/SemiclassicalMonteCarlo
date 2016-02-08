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
    f.open("test.dat");

    for (int i = 0; i < 30; ++i) {
        double e = 2. * i / 30;
        Result res = result(VAR_E0X, e);
        f << e << "\t" << res.current_mean.x << "\t" << res.current_mean.y
          << "\t" << res.current_std.x << "\t" << res.current_std.y << "\t"
          << res.n_opt << "\t" << res.n_ac << "\t" << res.tau << endl;
    }

    f.close();
}