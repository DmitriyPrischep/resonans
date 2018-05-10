#ifndef FUNCTIONAL
#define FUNCTIONAL
#include <QFile>
#include "target.h"
#include "fftw3.h"

int writeDataQt(QString path, Emission *emission);
int readDataQt(QString path, std::vector<Target> *targets);
void reflectionCoefficient(float permittivity, float conduct, float elevat, float lamda, float *pol_h_x, float *pol_h_y, float *pol_v_x, float *pol_v_y);
void vio(Emission *emission, int indexTn, int indexReciver, int cntChannels, std::vector<float> windowFun);
void horizontalImitation(Emission *emission, std::vector<float> *windowFun, std::vector<Target> *targets);
void fft(struct signal in[], int n, int p);
int CFDN(Emission* emission);
int dopplerFiltration(Emission* emission, int countAE);
int detection(Emission* emission);
int calcBorder(Emission* emission, std::vector<struct azimuth> *azimuths);
int marksSelection(Emission *emission, std::vector<struct azimuth> *azimuths, std::vector<struct mark>* marks);
int evalCoordinatesMarks(Emission *emission, std::vector<struct mark>* marks, std::vector<Target> *targets);
void verticalImitation(Emission *emission, std::vector<float> *windowFun, std::vector<Target> *targets);
//void CFDN(std::vector<struct _signal>* array, float *elevat, float *amplitude);
void MX(std::vector<struct _signal>* arrayRealSignal, float *elevat );
void calculateElevation(Emission *emission, std::vector<struct mark>* marks, std::vector<Target> *targets);

#endif // FUNCTIONAL

