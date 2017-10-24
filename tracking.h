#ifndef TRACKING_H
#define TRACKING_H
//#define CountEmission 64        // Количество излучений радиосигналов
//#define CountReceiver 16        // Количество приемников (антенн)
//#define CountChannels 256        // Количество каналов дальности


//struct signal{
//    double x;
//    double y;
//};

//struct reciver{
//    struct signal signalsArr[CountChannels];
//};

//typedef struct tracking{
//    struct reciver recivers[CountReceiver];
//}Tracking;

static const int countChannels = 256;
static const int countRecivers = 16;
static const int countEmission = 64;


struct signal{
    double x;
    double y;
};
struct reciver{
    struct signal signalsArr[countChannels];
};

typedef struct emission{
    struct reciver recivers[countRecivers];
} Emission;


class Tracking
{
private:

public:
    Tracking();








//    class Temp{
//    public:
//        float x;
//        float y;
//        Temp();
//    };
//    Temp* t = new Temp[10];



//    class Reciver
//    {
//    public:
//        Reciver();
//        class Signal
//        {
//        public:
//            Signal();
//            double x;
//            double y;
//        };
//        Signal* signalsArr = new Signal[countChannels];
//    };
//    Reciver* recivers = new Reciver[countRecivers];


//    struct signal{
//        double x;
//        double y;
//    };
//    struct reciver{
//        struct signal signalsArr[countChannels];
//    };

//    struct reciver recivers[countRecivers];
};

#endif // TRACKING_H
