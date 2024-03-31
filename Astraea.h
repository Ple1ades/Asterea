#include "lib/Quartic/quartic.h"
#include <vector>
#include "TaylorConstants.h"

std::vector<double> taylor(std::vector<double> taylorFactorsN, std::vector<double> taylorFactorsZ, double deltaX, double deltaY, double velocity);
double angleFromPosition(int deltaX, int deltaY, int velocity);