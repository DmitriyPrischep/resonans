#ifndef EMISSION_H
#define EMISSION_H
#include "configure.h"
#include "stddef.h"

struct azimuth{
    float avg;     //Среднее арифметическое по азимуту
    float sigma;   //Среднеквадратичное отклонение
    float border;  //Порог
};

struct mark{
    int i;
    int j;
    int k;
    float value;
    int size;
};


struct _signal{
    float x;
    float y;
};
struct reciver{
    struct _signal signalsArr[countChannels];
};

struct radioSignal{
    struct reciver recivers[countRecivers];
};

class Emission
{
public:
    Emission();
    Emission(int sizeEmission);
    ~Emission();
    struct radioSignal* allocateEmission(int size);
    struct radioSignal* data;
};

#endif // EMISSION_H
