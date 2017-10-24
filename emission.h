#ifndef EMISSION_H
#define EMISSION_H

static const int countChannels = 256;
static const int countRecivers = 16;
static const int countEmission = 64;

class Emission
{
public:
    Emission();
    struct signal{
        double x;
        double y;
    };
    struct signal sig[countChannels];

//    struct reciver{
//        struct signal signalsArr[countChannels];
//    };

//    struct radioSignal{
//        struct reciver recivers[countRecivers];
//    } ;

//    radioSignal data[countEmission];
};

#endif // EMISSION_H
