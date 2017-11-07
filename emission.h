#ifndef EMISSION_H
#define EMISSION_H
#include "configure.h"
#include "stddef.h"

struct azimuth{
    double avg;     //Среднее арифметическое по азимуту
    double sigma;   //Среднеквадратичное отклонение
    double border;  //Порог
};


struct signal{
    double x;
    double y;
};
struct reciver{
    struct signal signalsArr[countChannels];
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
