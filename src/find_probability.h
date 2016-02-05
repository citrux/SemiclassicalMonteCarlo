#include <vector>
#include <omp.h>
#include "material_specific.h"

std::vector<Point> momentums_with_energy_in_direction(double psi,
                                                      double energy_value);

void set_probabilities();

double get_probability(double energy);

double distrib_function(double p, double psi);