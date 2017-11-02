#ifndef FUNCTIONAL
#define FUNCTIONAL
#include <QFile>
#include "target.h"
//#include "emission.h"

int readData(std::vector<Target> *targets);
int readDataQt(QString path, std::vector<Target> *targets);
void imitationTargets(Emission *emis, std::vector<double> windowFun, double DF, double Tn, int lenProbeSignal, double noise, std::vector<Target> targets);
#endif // FUNCTIONAL

