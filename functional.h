#ifndef FUNCTIONAL
#define FUNCTIONAL
#include <QFile>
#include "target.h"

int readData(std::vector<Target> *targets);
int readDataQt(QString path, std::vector<Target> *targets);
void reflectionCoefficient(double permittivity, double conduct, double elevat, double lamda, double *pol_h_x, double *pol_h_y);
void vio(Emission *emission, int indexTn, int indexReciver, int cntChannels, std::vector<double> windowFun);
void horizontalImitation(Emission *emission, std::vector<double> *windowFun, std::vector<Target> *targets);
void fft(struct signal in[], int n, int p);
int CFDN(Emission* emission);
int dopplerFiltration(Emission* emission);
int detection(Emission* emission);
int calcBorder(Emission* emission, std::vector<struct azimuth> *azimuths);
int marksSelection(Emission *emission, std::vector<struct azimuth> *azimuths, std::vector<struct mark>* marks);
int evalCoordinatesMarks(Emission *emission, std::vector<struct mark>* marks, std::vector<Target> *targets);
void verticalImitation(Emission *emission, std::vector<double> *windowFun, std::vector<Target> *targets);
void CFDN(std::vector<struct _signal>* array, double *elevat, double *amplitude);
void MX(std::vector<struct _signal>* arrayRealSignal, double *elevat );

#endif // FUNCTIONAL

