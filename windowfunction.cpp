#include "windowfunction.h"
#include "configure.h"

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


void Windowfunction::funcHamming(std::vector<double> *winFunction, int size){
    assert(winFunction);
    for(int i = 0; i < size; i++){
        winFunction->push_back(0.54 + 0.46 * cosf(2 * M_PI * (0.5+i-size*0.5)/size));
//        array[i] = 0.54 + 0.46 * cosf(2 * M_PI * (0.5+i-size*0.5)/size);
    }
}
