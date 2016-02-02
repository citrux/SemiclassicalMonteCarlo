#include <vector>
#include <omp.h>
#include "material_specific.h"

std::vector<Point> momentums_with_energy_in_direction(double psi,
                                                      double energy_value,
                                                      const Params & params);

void set_probabilities(Params & params);

double get_probability(Point p, const Params & params);

double distrib_function(double p, double psi, const Params & params);