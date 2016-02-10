#pragma once
#include <string>
#include <vector>

void load_arrays(const std::string & filename, std::vector<double *> arrays,
                 int n);
void save_arrays(const std::string & filename, std::vector<double *> arrays,
                 int n);
