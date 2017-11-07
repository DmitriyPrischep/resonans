#include <iostream>
#include <fstream>
#include <QFile>
#include <QDebug>
#include <QDataStream>
#include "assert.h"
#include "target.h"
#include "emission.h"
#include "configure.h"
//#include "mainwindow.h"
#include "windowfunction.h"


int readData(std::vector<Target> *targets){
    const int lenString = 70;
    const int cntStrings = 4;
    const char ch = '\n';
    char mass[lenString][cntStrings];

    std::ifstream fs("targets1.txt", std::ios::in | std::ios::binary);
    if(!fs) return 1;

    for(int i = 0; i < cntStrings; i++){
        fs.getline(mass[i], lenString-1, ch);
        std::cerr << "String" << i + 1 << " = " << mass[i] <<endl;
    }
    qDebug() << "File is read successfully\n";
    fs.close();
    return 0;
}

int readDataQt(QString path, std::vector<Target> *targets){
    QFile file(path);
    if(!file.exists()){
        qDebug() << "File '" + path +  "' does not exist";
        return 2;
    }
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)){
        qDebug() << "Error reading file " + path;
        return 1;
    }
    file.readLine();    // Пропускаем первую строку
    while(!file.atEnd()){
        Target tar;
        QString line = file.readLine();
        QStringList list = line.split('\t');
        tar.setA(list.at(0).trimmed().toFloat());
        tar.setR(list.at(1).trimmed().toFloat());
        tar.setB(list.at(2).trimmed().toFloat());
        tar.setV(list.at(3).trimmed().toFloat());
        tar.setF(list.at(4).trimmed().toFloat());
        tar.setG(list.at(5).trimmed().toFloat());
        tar.setU(list.at(6).trimmed().toFloat());
        qDebug() << list[0] << list[1] << list[2] << list[3] << list[4] << list[5] << list[6] << endl;
        (*targets).push_back(tar);
    }

    qDebug() << "File is read successfully\n";
    file.close();
    return 0;
}


/// \brief - Генератор нормального шума
/// \param [in, out] *nu1 - массив коэффициентов
/// \param [in, out] *nu2 - уровень боковиков
void  noiseGeneration(double *x, double *y){
    assert(x);
    assert(y);
    *x = (double)(rand())/RAND_MAX*(3 - (-3))  -3;
    *y = (double)(rand())/RAND_MAX*(3 - (-3))  -3;
}



/// \brief - Коэффициент отражения от поверхности
/// \param [in] permittivity - Диэлектрическая проницаемость
/// \param [in] conduct - Проводимость
/// \param [in] elevat - Улог места
/// \param [in] lamda - Длина волны
/// \param [in] pol_v_x - Вертикальная поляризация по Х
/// \param [in] pol_v_y - Вертикальная поляризация по Y
/// \param [in] pol_h_x - Горизонтальная поляризация по Х
/// \param [in] pol_h_y - Горизонтальная поляризация по Y
void reflectionCoefficient(double permittivity, double conduct, double elevat, double lamda, double *pol_h_x, double *pol_h_y)/*float *pol_v_x, float *pol_v_y,*/ {
    double angle = elevat * M_PI / 180;
    double sin_angle = sin(angle);

    double par1 = permittivity - pow(cos(angle), 2);

    double par2 = -60 * conduct * lamda;
    double tmp = sqrt(sqrt(pow(par1, 2) + pow(par2, 2)));

//    float arctan = cos(atan2(par2, par1) / 2.0);
    double arctan = atan2(par2, par1);
    double px = tmp * cos(arctan),
            py = tmp * sin(arctan);

//    double f10x = permittivity * sin_angle - px;
//    double f10y = par2 * sin_angle - py;

    double pxx = permittivity * sin_angle + px,
            pyy = par2 * sin_angle + py;

    // Вертикаль
    double zn = pow(pxx, 2) + pow(pyy, 2);
/*
    double f1x = (f10x * pxx + f10y * pyy)/zn,
            f1y = (-f10x * pyy + f10y * pxx)/zn;
    *pol_v_x = sqrt(pow(f1x, 2) + pow(f1y, 2));
    *pol_v_y = - atan2(f1y, f1x) * 180.0 / M_PI;
*/
    // Горизонталь
    double f20x = sin_angle - px,
            f20y = - py;

    pxx = sin_angle + px;
    pyy = py;

    zn = pow(pxx, 2) + pow(pyy, 2);

    double f2x = (f20x * pxx + f20y * pyy) / zn,
            f2y = (-f20x * pyy + f20y * pxx) / zn;

    *pol_h_x = sqrt(pow(f2x, 2) + pow(f2y, 2));
    *pol_h_y = atan2(f2y, f2x) * 180.0 / M_PI;
    return;
}



/// \brief - Внутриимпульсная обработка сигнала
/// \param [in, out] sw[] - Свертка сигнала
/// \param [in] cnt_channels - Количество каналов дальностей
/// \param [in] window_fun[] - Оконная функция
/// \param [in] F_Deviation - Девиация частоты
/// \param [in] tauZS - длительность зондирующего сигнала
void vio(Emission *emission, int indexTn, int indexReciver, int cntChannels, std::vector<double> windowFun, int lenProbeSignal){
    assert(lenProbeSignal > 0);
    float arg = 0;
//    int kl = (int) (F_Deviation * tauZS + 0.1);
    int kl = (int) (lenProbeSignal + 0.1);

    for (int i = 0; i < cntChannels; i++){
        double c = 0, d = 0;
        for (int j = 0; j < lenProbeSignal; j++){
            int k = i + j;
            if (k >= cntChannels)
                continue;

            arg = (M_PI / kl) * pow((1 + j - (kl + 1) / 2.0), 2);

            c += windowFun.at(j) * (  emission->data[indexTn].recivers[indexReciver].signalsArr[k].x * cos(arg)
                                    + emission->data[indexTn].recivers[indexReciver].signalsArr[k].y * sin(arg));

            d += windowFun.at(j) * (- emission->data[indexTn].recivers[indexReciver].signalsArr[k].x * sin(arg)
                                    + emission->data[indexTn].recivers[indexReciver].signalsArr[k].y * cos(arg));
//            c += window_fun[j] * (arr[indexTn].reciver[indexReciver].signals[k].x * cos(arg)
//                                  + arr[indexTn].reciver[indexReciver].signals[k].y * sin(arg));
//            d += window_fun[j] * (-arr[indexTn].reciver[indexReciver].signals[k].x * sin(arg)
//                                  + arr[indexTn].reciver[indexReciver].signals[k].y * cos(arg));
        }
        emission->data[indexTn].recivers[indexReciver].signalsArr[i].x = c;
        emission->data[indexTn].recivers[indexReciver].signalsArr[i].y = d;
//        arr[indexTn].reciver[indexReciver].signals[i].x = c;
//        arr[indexTn].reciver[indexReciver].signals[i].y = d;
    }
}

void imitationTargets(Emission *emis, std::vector<double> windowFun, double DF, double Tn, int lenProbeSignal, double noise, std::vector<Target> targets){
    double intervals[countRecivers];
    double heights[countRecivers];

    for(int j = 0; j < countRecivers; j++){
        intervals[j] = distanceAntenna * j;
        heights[j] = heightAntenna;
    }

    double lambda = (double)300/60;
    double dt = (double)1/DF;
    const double piOn180 = (double)M_PI/180;

    for(int i = 0; i < countEmission; i++){
        for(int j = 0; j < countRecivers; j++){
            for(int k = 0; k < countChannels; k++){
                double x, y; // X,Y шума по нормальному распределению
                noiseGeneration(&x, &y);
                emis->data[i].recivers[j].signalsArr[k].x = x * noise;
                emis->data[i].recivers[j].signalsArr[k].y = y * noise;

                for(int t = 0; t < countTargets; t++){
                    double a = (targets.at(t).getA() <= -100) ? 0 : pow(10, (double)targets.at(t).getA()/20);
                    double b = targets.at(t).getB();
                    double g = targets.at(t).getG();
                    double r = targets.at(t).getR();
                    double f = targets.at(t).getF();
                    double v = targets.at(t).getV();
                    double u = targets.at(t).getU();

                    double vt = r + (v * i * Tn) /1000.0 + 0.5 * (u * pow(i,2) * pow(Tn, 2)) /1000.0;
                    double time = vt / (150 * dt);
                    int begin = (int)floor(time + 0.49);
                    int end = begin + lenProbeSignal - 1;

                    if(k >= begin && k <= end && a > 0){
                        double phase = (2 * M_PI/lambda) * (heights[j] + heightSea) * sin(g * piOn180);
//                        double sumPhase = (2 * M_PI/lambda) * (2 * v * (i * Tn + (k - begin) * (dt/1000)) * intervals[j] * sin(b * piOn180)
//                                                               + (double)u/2 * (pow(i,2) * pow(Tn,2) + pow(k - begin, 2) * pow((dt/1000),2)))
//                                + M_PI / lenProbeSignal * pow(0.5 + k - time - (double)lenProbeSignal/2, 2) + f;


//                        double arg0 = 2 * M_PI/lambda * (heights[j] + HeightSea) * sin(g * piOn180);
                        double arg = 2 * M_PI/lambda * intervals[j] * sin(b * piOn180);
                        arg += 2 * M_PI/lambda * 2 * v * i * Tn;
                        arg += 2 * M_PI/lambda * 2 * v * (k - begin) * dt/1000.0;
                        arg += 2 * M_PI/lambda * (u * pow(i, 2) * pow(Tn, 2)/2);
                        arg += 2 * M_PI/lambda * (u * pow((k - begin), 2) * pow((dt/1000.0),2));
                        arg += M_PI / lenProbeSignal * pow((0.5 + k - time - (double)lenProbeSignal/2), 2) + f;

                        double polarPhase = 0;
                        double polarAmplitude = 0;
                        reflectionCoefficient(permittivity, conduct, g, lambda, &polarAmplitude, &polarPhase);
                        double phase_x = arg + phase,
                               phase_y = arg - phase + polarPhase * piOn180;

                        double m1 = cos(phase_x);
                        m1 += polarAmplitude * cos(phase_y);
                        double m2 = (sin(phase_x));
                        m2 += polarAmplitude * sin(phase_y);

                        emis->data[i].recivers[j].signalsArr[k].x += m1 * windowFun.at(k - begin) * a;
                        emis->data[i].recivers[j].signalsArr[k].y += m2 * windowFun.at(k - begin) * a;
//                        data[i].reciver[j].signals[k].x += m1 * window_fun[k - begin] * a;
//                        data[i].reciver[j].signals[k].y += m2 * window_fun[k - begin] * a;
                    }

                }
            }
        }
        for(int j = 0; j < countRecivers; j++){
            vio(emis, i, j, countChannels, windowFun, lenProbeSignal);
        }

    }
}


void swapping(struct signal array[], int n){
    //assert()
    // !in  n>1
    int part = n / 2;
    for(int i = 0; i < part; i++){
        struct signal temp;
        temp.x = array[i].x;
        array[i].x = array[part + i].x;
        array[part + i].x = temp.x;

        temp.y = array[i].y;
        array[i].y = array[part + i].y;
        array[part + i].y = temp.y;
    }
}

/// \brief - Быстрое преобразование Фурье
/// \param [in, out] in[] - входной массив комплексных чисел
/// \param [in] n - Размерность массива
/// \param [in] p - Направление (p > 0 - Прямое, p == 0 - Обратное)
void fft(struct signal in[], int n, int p){
    int half = n / 2;
    int j = 2, i = 1, k = 1, m, c;
    while (j <= half){
        i = k + half;

        double temp;
        temp = in[j-1].x;
        in[j-1].x = in[i-1].x;
        in[i-1].x = temp;

        temp = in[j-1].y;
        in[j-1].y = in[i-1].y;
        in[i-1].y = temp;

        j++;
        m = half / 2;

        while (k > m){
            k -= m;
            m /= 2;
        }
        k += m;

        if (k > j){
            temp = in[j-1].x;
            in[j-1].x = in[k-1].x;
            in[k-1].x = temp;

            temp = in[j-1].y;
            in[j-1].y = in[k-1].y;
            in[k-1].y = temp;

            c = j + half + 1;
            i = k + half + 1;

            temp = in[c-1].x;
            in[c-1].x = in[i-1].x;
            in[i-1].x = temp;

            temp = in[c-1].y;
            in[c-1].y = in[i-1].y;
            in[i-1].y = temp;
        }
        j++;
    }
    i = 1;
    double t = M_PI;
    if (p > 0)
        t = -t;

    while (i < n){
        double im = 0,
                re = 1,
                Wpi = sin(t),
                Wpr = cos(t),
                temp;
        t *= 0.5;
        int cnt = i;
        i+=i;
        for (int h = 1; h <= cnt; h++){
            for(int j = h; j <= n; j += i){
               k = j + cnt;
               double a = in[k-1].x;
               double b = in[k-1].y;

               temp = a * re - b * im;
               in[k-1].x = in[j-1].x - temp;
               in[j-1].x += temp;

               temp = b * re + a * im;
               in[k-1].y = in[j-1].y - temp;
               in[j-1].y += temp;

            }
            temp = Wpr * re - Wpi * im;
            im = Wpr * im + Wpi * re;
            re = temp;
        }
    }
    swapping(in, n);
}


int CFDN(Emission* emission){
    if(!emission)
        return 1;
    std::vector<double> windowFun;
    Windowfunction::funcHamming(&windowFun, countRecivers);
    for(int i = 0; i < countEmission; i++){
        for(int k = 0; k < countChannels; k++){
            struct signal temp[countRecivers];

            for(int j = 0; j < countRecivers; j++){
                temp[j].x = emission->data[i].recivers[j].signalsArr[k].x * windowFun.at(j);
                temp[j].y = emission->data[i].recivers[j].signalsArr[k].y * windowFun.at(j);
//                temp[j].x = data[i].reciver[j].signals[k].x * winCFDN[j];
//                temp[j].y = data[i].reciver[j].signals[k].y * winCFDN[j];
            }

            fft(temp, countRecivers, 1);

            for(int j = 0; j < countRecivers; j++){
                emission->data[i].recivers[j].signalsArr[k].x = temp[j].x;
                emission->data[i].recivers[j].signalsArr[k].y = temp[j].y;
//                data[i].reciver[j].signals[k].x = temp[j].x;
//                data[i].reciver[j].signals[k].y = temp[j].y;
            }
        }
    }
    return 0;
}

int dopplerFiltration(Emission* emission){
    if(!emission)
        return 1;
    std::vector<double> windowFun;
    Windowfunction::funcHamming(&windowFun, countEmission);

    for(int j = 0; j < countRecivers; j++){
        for(int k = 0; k < countChannels; k++){
            struct signal temp[countEmission];

            for(int i = 0; i < countEmission; i++){
                temp[i].x = emission->data[i].recivers[j].signalsArr[k].x * windowFun.at(i);
                temp[i].y = emission->data[i].recivers[j].signalsArr[k].y * windowFun.at(i);
//                temp[i].x = data[i].reciver[j].signals[k].x * winDopler[i];
//                temp[i].y = data[i].reciver[j].signals[k].y * winDopler[i];
            }

            fft(temp, countEmission, 1);

            for(int i = 0; i < countEmission; i++){
                emission->data[i].recivers[j].signalsArr[k].x = temp[i].x;
                emission->data[i].recivers[j].signalsArr[k].y = temp[i].y;
//                data[i].reciver[j].signals[k].x = temp[i].x;
//                data[i].reciver[j].signals[k].y = temp[i].y;
            }
        }
    }
    return 0;
}


/// \brief - Вычисление порога
/// \param [in] data - входной массив амплитуд
/// \param [in, out] azimuths  - массив СКО и СА по каждому азимуту
/// \param [in] cntDoplers - количество доплеров
/// \param [in] cntAzimuth - количество направлений азимутов
/// \param [in] cntChannels - количество каналов дальности
int detection(Emission* emission){
    if(!emission)
        return 1;
    for(int i = 0; i < countEmission; i++){
        for(int j = 0; j < countRecivers; j++){
            for(int k = 0; k < countChannels; k++){
                emission->data[i].recivers[j].signalsArr[k].x = sqrt(pow(emission->data[i].recivers[j].signalsArr[k].x, 2) +
                                                                     pow(emission->data[i].recivers[j].signalsArr[k].y, 2));
                emission->data[i].recivers[j].signalsArr[k].y = 0;
//                data[i].reciver[j].signals[k].x = sqrt(pow(data[i].reciver[j].signals[k].x, 2)
//                                                       + pow(data[i].reciver[j].signals[k].y, 2));
//                data[i].reciver[j].signals[k].y = 0;
            }
        }
    }
    return 0;
}

int calcBorder(Emission emission, struct azimuth* azimuths){
//    cntDopllers = countEmission;
//    cntAzimuths = countRecivers;
//    cntChannels = countChannels;
    if (!azimuths)
        return 1;
    for(int j = 0; j < countRecivers; j++){
        double avg = 0;
        for(int k = 0; k < countChannels; k++){
            for(int i = 0; i < countEmission; i++){
                avg += emission.data[i].recivers[j].signalsArr[k].x;
//                avg += data[i].reciver[j].signals[k].x;
            }
        }
        avg /= (countChannels * countEmission);
        azimuths[j].avg = avg;

        double sigma = 0;
        for(int k = 0; k < countChannels; k++){
            for(int i = 0; i < countEmission; i++){
                sigma += pow(emission.data[i].recivers[j].signalsArr[k].x - avg, 2);
//                sigma += pow(data[i].reciver[j].signals[k].x - avg, 2);
            }
        }
        azimuths[j].sigma = sqrt(sigma / (countChannels * countEmission));
        azimuths[j].border = alphaBorder * azimuths[j].avg + gammaBorder * azimuths[j].sigma;
    }
    return 0;
}











