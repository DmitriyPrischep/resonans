#ifndef CONFIGURE
#define CONFIGURE
#include <QSettings>

static const int countChannels = 256;       // Количество каналов дальности
static const int countRecivers = 16;        // Количество приемников (антенн)
static const int countEmission = 64;        // Количество излучений радиосигналов
static const int cntAntennas = 16;           // Количество антенн

static const int countTargets = 1;          // Количество целей
static const int countParametrs = 7;        // Количество параметров цели
static const int heightSea = 0;             // Высота станции относительно уровня моря
static const float permittivity = 20.0;     // Диэлектрическая проницаемость влажной почвы
static const float conduct = 0.01;          // Проводимость влажной почвы
static const float heightAntenna = 6.25;    // Высота антенн приемников
static const float distanceAntenna = 2.3;   // Расстояние между антеннами
static const float noise = 1.0;             // Уровень шума
static const int frequencyDeviation = 100;  // Девиация частоты
static const float durationWave = 0.01;     // Период повторения зондирующего сигнала
static const int lenProbeSignal = 45;       // Длина зондирующего сигнала
static const float alphaBorder = 4.5;       // Альфа коэффициент порога
static const float gammaBorder = 10;        // Гамма коэффициент порога
static const int countVertRecivers = 8;     // Количество антенных элементов на вертикальной антенной решетке
static const float anten[countVertRecivers] // Высоты каждого элемента
    = {2.5, 5, 7.5, 10, 13.8, 17.2, 20.6, 24};

class Settings{
public:
    static void setValues(QSettings *settings);
    static void getValues(QSettings *settings);
};

#endif // CONFIGURE

