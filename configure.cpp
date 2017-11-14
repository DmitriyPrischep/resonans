#include "configure.h"

void Settings::setValues(QSettings *settings){
    settings->setValue("emission/countChannels", 256);
    settings->setValue("emission/countRecivers", 16);
    settings->setValue("emission/countEmission", 512);
    settings->setValue("targets/countTargets", 4);
    settings->setValue("targets/countParametrs", 7);
    settings->setValue("radar/heightSea", 0);
    settings->setValue("radar/permittivity", 20.0);
    settings->setValue("radar/conduct", 0.01);
    settings->setValue("radar/heightAntenna", 6.25);
    settings->setValue("radar/distanceAntenna", 2.3);
    settings->setValue("radar/cntAntennas", 16);
    settings->setValue("radar/noise", 1.0);
    settings->setValue("radar/frequencyDeviation", 100);
    settings->setValue("radar/durationWave", 0.01);
    settings->setValue("radar/lenProbeSignal", 45);
    settings->setValue("radar/alphaBorder", 2.1);
    settings->setValue("radar/gammaBorder", 5);
    settings->sync();
}
