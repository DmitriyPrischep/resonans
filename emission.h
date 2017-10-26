#ifndef EMISSION_H
#define EMISSION_H
#include <stddef.h>

static const int countChannels = 256;
static const int countRecivers = 16;
static const int countEmission = 512;

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
