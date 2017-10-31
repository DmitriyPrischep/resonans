#ifndef WINDOWFUNCTION_H
#define WINDOWFUNCTION_H
#include "math.h"
#include "assert.h"
#include <vector>

class Windowfunction{
public:
    static void funcHammingTukey(std::vector<double> *winFunction);
    static void funcHamming(double *array, int size);
    static void funcTukey(double *array, int size);
    static const int lenProbeSignal;
};

const int Windowfunction::lenProbeSignal = 45;

void Windowfunction::funcHammingTukey(std::vector<double> *winFunction){
    assert(winFunction);
    (*winFunction).reserve(lenProbeSignal);
    int border = 0.1 * lenProbeSignal,
            sizeW = border * 2;
    double arr[sizeW];
    for(int i = 0; i < sizeW; i++){
        arr[i] = 0.54 + 0.46 * cosf(2 * M_PI * (0.5+i-sizeW*0.5)/sizeW);
    }
    for(int i = 0; i < lenProbeSignal; i++){
        winFunction->push_back(1);
    }
    for(int i = 0; i < sizeW/2; i++){
        winFunction->at(i) = winFunction->at(lenProbeSignal - 1 - i) = arr[i];
    }
}

#endif // WINDOWFUNCTION_H
