#ifndef WINDOWFUNCTION_H
#define WINDOWFUNCTION_H
#include "math.h"
#include "assert.h"
#include <vector>

class Windowfunction{
public:
    static void funcHammingTukey(std::vector<double> *winFunction);
    static void funcHamming(std::vector<double> *winFunction, int size);
    static void funcTukey(double *array, int size);
    static const int lenProbeSignal;
};

#endif // WINDOWFUNCTION_H
