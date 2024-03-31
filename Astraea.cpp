#include <iostream>
#include <vector>
#include <complex>
#include <map>
#include "Astraea.h"
#define PI 3.14159265358979323846

std::vector<double> taylor(std::vector<double> taylorFactorsN, std::vector<double> taylorFactorsZ, double deltaX, double deltaY, double velocity)
{
    std::vector<double> a(5);
    int factorial = 1;
    for (int i = 0; i < 5; ++i)
    {
        factorial *= std::max(1, i);
        a[4 - i] = (taylorFactorsN[4 - i] * deltaX - taylorFactorsZ[4 - i] * deltaX * deltaX / (velocity * velocity)) / factorial;
    }
    a[4] -= deltaY;
    return a;
}
double validateSolutions(std::vector<std::complex<double>> solutions, double ang, double A)
{
    double minSol = 1000;
    for (int i = 0; i < 4; ++i)
    {
        // if (solutions[i].real() > 0)
        // {
        //     std::cout << "_" << solutions[i].real() << std::endl;
        // }
        if (solutions[i].real() > 0 && solutions[i].real() < minSol)
        {

            minSol = solutions[i].real() + ang;
        }
    }
    return minSol;
}
double validateAngle(double angle, double x, double v)
{
    return x * 16 * (angle * PI / 180) * (PI - angle * PI / 180) * (PI * PI + (angle * PI / 180) * (angle * PI / 180)) / ((5 * PI * PI - 4 * angle * PI / 180 * (PI - angle * PI / 180)) * (PI * PI - 4 * (angle * PI / 180) * (angle * PI / 180))) - 9.81 / 2 * (x * (PI * PI + (angle * PI / 180) * (angle * PI / 180)) * (x * (PI * PI + (angle * PI / 180) * (angle * PI / 180)) / (v * (PI * PI - 4 * (angle * (PI / 180) * angle * PI / 180)) * v * (PI * PI - 4 * (angle * (PI / 180) * angle * PI / 180)))));
}
double angleFromPosition(int deltaX, int deltaY, int velocity)
{
    std::map<int, std::vector<double>> diffs;
    std::vector<double> deg_15 = taylor(std::vector<double>{AN1, BN1, CN1, DN1, EN1}, std::vector<double>{AZ1, BZ1, CZ1, DZ1, EZ1}, deltaX, deltaY, velocity);
    std::complex<double> *deg_15sol = solve_quartic(deg_15[1] / deg_15[0], deg_15[2] / deg_15[0], deg_15[3] / deg_15[0], deg_15[4] / deg_15[0]);
    std::vector<std::complex<double>> deg_15sols = {deg_15sol[0], deg_15sol[1], deg_15sol[2], deg_15sol[3]};
    double deg_15ang = validateSolutions(deg_15sols, 15, deg_15[0]);
    diffs[15] = std::vector<double>{validateAngle(deg_15ang, deltaX, velocity), deg_15ang};
    delete[] deg_15sol;
    std::vector<double> deg_45 = taylor(std::vector<double>{AN2, BN2, CN2, DN2, EN2}, std::vector<double>{AZ2, BZ2, CZ2, DZ2, EZ2}, deltaX, deltaY, velocity);
    std::complex<double> *deg_45sol = solve_quartic(deg_45[1] / deg_45[0], deg_45[2] / deg_45[0], deg_45[3] / deg_45[0], deg_45[4] / deg_45[0]);
    std::vector<std::complex<double>> deg_45sols = {deg_45sol[0], deg_45sol[1], deg_45sol[2], deg_45sol[3]};
    double deg_45ang = validateSolutions(deg_45sols, 45, deg_45[0]);
    diffs[45] = std::vector<double>{validateAngle(deg_45ang, deltaX, velocity), deg_45ang};
    delete[] deg_45sol;
    std::vector<double> deg_60 = taylor(std::vector<double>{AN3, BN3, CN3, DN3, EN3}, std::vector<double>{AZ3, BZ3, CZ3, DZ3, EZ3}, deltaX, deltaY, velocity);
    std::complex<double> *deg_60sol = solve_quartic(deg_60[1] / deg_60[0], deg_60[2] / deg_60[0], deg_60[3] / deg_60[0], deg_60[4] / deg_60[0]);
    std::vector<std::complex<double>> deg_60sols = {deg_60sol[0], deg_60sol[1], deg_60sol[2], deg_60sol[3]};
    double deg_60ang = validateSolutions(deg_60sols, 60, deg_60[0]);
    diffs[60] = std::vector<double>{validateAngle(deg_60ang, deltaX, velocity), deg_60ang};
    delete[] deg_60sol;
    std::vector<double> deg_75 = taylor(std::vector<double>{AN4, BN4, CN4, DN4, EN4}, std::vector<double>{AZ4, BZ4, CZ4, DZ4, EZ4}, deltaX, deltaY, velocity);
    std::complex<double> *deg_75sol = solve_quartic(deg_75[1] / deg_75[0], deg_75[2] / deg_75[0], deg_75[3] / deg_75[0], deg_75[4] / deg_75[0]);
    std::vector<std::complex<double>> deg_75sols = {deg_75sol[0], deg_75sol[1], deg_75sol[2], deg_75sol[3]};
    double deg_75ang = validateSolutions(deg_75sols, 75, deg_75[0]);
    diffs[75] = std::vector<double>{validateAngle(deg_75ang, deltaX, velocity), deg_75ang};
    delete[] deg_75sol;
    double minAngle = 180;
    double minError = 100;
    // for (int i = 0; i < deg_15.size(); ++i)
    // {
    //     std::cout << deg_15[i] << std::endl;
    // }
    for (auto const &diff : diffs)
    {
        // std::cout << abs(diff.second[0] - deltaY) << std::endl;
        if (abs(diff.second[0] - deltaY) < minError)
        {
            minAngle = diff.second[1];
            minError = abs(diff.second[0] - deltaY);
        }
    }

    return minAngle;
}
