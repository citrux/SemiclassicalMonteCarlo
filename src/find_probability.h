#pragma once

#include <vector>
#include <omp.h>
#include "linalg.h"

std::vector<Point> momentums_with_energy_in_direction(double psi,
                                                      double energy_value);

void set_probabilities();

double get_probability(double energy);