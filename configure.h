#ifndef CONFIGURE
#define CONFIGURE

static const int countTargets = 4;          // Количество целей

static const int countChannels = 256;       // Количество излучений радиосигналов
static const int countRecivers = 16;        // Количество приемников (антенн)
static const int countEmission = 512;       // Количество каналов дальности

static const int countParametrs = 7;        // Количество параметров цели
static const int heightSea = 0;             // Высота станции относительно уровня моря
static const int permittivity = 20.0;       // Диэлектрическая проницаемость влажной почвы
static const int conduct = 0.01;            // Проводимость влажной почвы
static const int heightAntenna = 6.25;      // Высота антенн приемников
static const int distanceAntenna = 2.3;     // Расстояние между антеннами
static const int cntAntennas = 16;          // Количество антенн
static const int noise = 1.0;               // Уровень шума
static const int frequencyDeviation = 100;  // Девиация частоты
static const int durationWave = 0.01;       // Период повторения зондирующего сигнала
static const int lenProbeSignal = 45;       // Длина зондирующего сигнала
static const int alphaBorder = 2.1;         // Альфа коэффициент порога
static const int gammaBorder = 5;           // Гамма коэффициент порога

#endif // CONFIGURE

