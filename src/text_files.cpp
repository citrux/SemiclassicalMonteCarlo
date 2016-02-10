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