#ifndef WINDOWFUNCTION_H
#define WINDOWFUNCTION_H
#include "math.h"
#include "assert.h"
#include <vector>

class Windowfunction{
public:
    static void funcHammingTukey(std::vector<float> *winFunction);
    static void funcHamming(std::vector<float> *winFunction, int size);
    static void funcTukey(float *array, int size);
};

#endif // WINDOWFUNCTION_H
