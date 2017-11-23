#include <iostream>
#include <fstream>
#include <QFile>
#include <QDebug>
#include <QDataStream>
#include "assert.h"
#include "target.h"
#include "emission.h"
#include "configure.h"
#include "windowfunction.h"


//int readData(std::vector<Target> *targets){
//    const int lenString = 70;
//    const int cntStrings = 4;
//    const char ch = '\n';
//    char mass[lenString][cntStrings];

//    std::ifstream fs("targets1.txt", std::ios::in | std::ios::binary);
//    if(!fs) return 1;

//    for(int i = 0; i < cntStrings; i++){
//        fs.getline(mass[i], lenString-1, ch);
//        std::cerr << "String" << i + 1 << " = " << mass[i] <<endl;
//    }
//    qDebug() << "File is read successfully\n";
//    fs.close();
//    return 0;
//}

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
/// \param [in, out] *emission - Свертка сигнала
/// \param [in] indexTn - Индекс канала дальности
/// \param [in] indexReciver - Индекс антенного элемента
/// \param [in] windowFun - Оконная функция
void vio(Emission *emission, int indexTn, int indexReciver, int cntChannels, std::vector<double> *windowFun){
    /////////////////////////////////// ДОПИСАТЬ ASSERT
//    assert();
    float arg = 0;
    int kl = (int) (lenProbeSignal + 0.1);

    for (int i = 0; i < cntChannels; i++){
        double c = 0, d = 0;
        for (int j = 0; j < lenProbeSignal; j++){
            int k = i + j;
            if (k >= cntChannels)
                continue;

            arg = (M_PI / kl) * pow((1 + j - (kl + 1) / 2.0), 2);

            c += windowFun->at(j) * (  emission->data[indexTn].recivers[indexReciver].signalsArr[k].x * cos(arg)
                                    + emission->data[indexTn].recivers[indexReciver].signalsArr[k].y * sin(arg));

            d += windowFun->at(j) * (- emission->data[indexTn].recivers[indexReciver].signalsArr[k].x * sin(arg)
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


/// \brief - Имитатор для вертикальной антенной решетки
/// \param [in, out] *emission - Свертка сигнала
/// \param [in] windowFun - Оконная функция
/// \param [in] DF - Частота девиации
/// \param [in] Tn - Период повторения сигнала
/// \param [in] noise - Уровень шума (Sigma)
/// \param [in] targets - Массив целей
void imitationTargets(Emission *emission, std::vector<double> *windowFun, std::vector<Target> *targets){
    double intervals[countRecivers];
    double heights[countRecivers];

    for(int j = 0; j < countRecivers; j++){
        intervals[j] = distanceAntenna * j;
        heights[j] = heightAntenna;
    }

    double lambda = (double)300/60;
    double dt = (double)1/frequencyDeviation;
    const double piOn180 = (double)M_PI/180;

    for(int i = 0; i < countEmission; i++){
        for(int j = 0; j < countRecivers; j++){
            for(int k = 0; k < countChannels; k++){
                double x, y; // X,Y шума по нормальному распределению
                noiseGeneration(&x, &y);
                emission->data[i].recivers[j].signalsArr[k].x = x * noise;
                emission->data[i].recivers[j].signalsArr[k].y = y * noise;

                for(int t = 0; t < countTargets; t++){
                    double a = (targets->at(t).getA() <= -100) ? 0 : pow(10, (double)targets->at(t).getA()/20);
                    double b = targets->at(t).getB();
                    double g = targets->at(t).getG();
                    double r = targets->at(t).getR();
                    double f = targets->at(t).getF();
                    double v = targets->at(t).getV();
                    double u = targets->at(t).getU();


                    double vt = r + (v * i * durationWave) /1000.0 + 0.5 * (u * pow(i,2) * pow(durationWave, 2)) /1000.0;
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
                        arg += 2 * M_PI/lambda * 2 * v * i * durationWave;
                        arg += 2 * M_PI/lambda * 2 * v * (k - begin) * dt/1000.0;
                        arg += 2 * M_PI/lambda * (u * pow(i, 2) * pow(durationWave, 2)/2);
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

                        emission->data[i].recivers[j].signalsArr[k].x += m1 * windowFun->at(k - begin) * a;
                        emission->data[i].recivers[j].signalsArr[k].y += m2 * windowFun->at(k - begin) * a;
//                        data[i].reciver[j].signals[k].x += m1 * window_fun[k - begin] * a;
//                        data[i].reciver[j].signals[k].y += m2 * window_fun[k - begin] * a;
                    }

                }
            }
        }
        for(int j = 0; j < countRecivers; j++){
            vio(emission, i, j, countChannels, windowFun);
        }

    }
}

///\brief - перестановка половинок массива
///\param [in, out] array[] - массив комплексных чисел
///\param [in] n - размер массива
void swapping(struct _signal array[], int n){
    assert(n > 0);
    assert(!array);
    int part = n / 2;
    for(int i = 0; i < part; i++){
        struct _signal temp;
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
void fft(struct _signal in[], int n, int p){
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


/// \brief - Центральное формирование диаграммы направленности
/// \param [in, out] *emission - Свертка сигнала
int CFDN(Emission* emission){
    if(!emission)
        return 1;
    std::vector<double> windowFun;
    Windowfunction::funcHamming(&windowFun, countRecivers);
    for(int i = 0; i < countEmission; i++){
        for(int k = 0; k < countChannels; k++){
            struct _signal temp[countRecivers];

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

/// \brief - Доплеровская фильтрация
/// \param [in, out] *emission - Свертка сигнала
int dopplerFiltration(Emission* emission){
    if(!emission)
        return 1;
    std::vector<double> windowFun;
    Windowfunction::funcHamming(&windowFun, countEmission);

    for(int j = 0; j < countRecivers; j++){
        for(int k = 0; k < countChannels; k++){
            struct _signal temp[countEmission];

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


/// \brief - Детекция сигнала
/// \param [in, out] *emission - Свертка сигнала
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

/// \brief - Вычисление порога
/// \param [in] data - Входной массив амплитуд
/// \param [in, out] azimuths  - Массив СКО и СА по каждому азимуту
int calcBorder(Emission emission, struct azimuth* azimuths){
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

/// \brief - Выделение координат отметок
/// \param [in, out] emission  - Свертка сигнала
/// \param [in] azimuths  - Массив СКО и СА по каждому азимуту
/// \param [in, out] marks - Выходной массив отметок, значение которого превышает порог border
int marksSelection(Emission *emission, struct azimuth* azimuths, std::vector<struct mark>* marks){
    if (!emission || !azimuths || !marks)
        return 1;
    for(int j = 1; j < countRecivers - 1; j++){
        for(int k = 1; k < countChannels - 1; k++){
            for(int i = 1; i < countEmission - 1; i++){
                if((emission->data[i].recivers[j].signalsArr[k].x > azimuths[j].border)  //Проверка с порогом
                        && (emission->data[i].recivers[j].signalsArr[k].x > emission->data[i].recivers[j+1].signalsArr[k].x)    // Проверка сверху
                        && (emission->data[i].recivers[j].signalsArr[k].x > emission->data[i].recivers[j].signalsArr[k+1].x)    // Проверка сзади
                        && (emission->data[i].recivers[j].signalsArr[k].x > emission->data[i+1].recivers[j].signalsArr[k].x)    // Проверка справа
                        && (emission->data[i].recivers[j].signalsArr[k].x > emission->data[i].recivers[j-1].signalsArr[k].x)    // Проверка сверху
                        && (emission->data[i].recivers[j].signalsArr[k].x > emission->data[i].recivers[j].signalsArr[k-1].x)    // Проверка спереди
                        && (emission->data[i].recivers[j].signalsArr[k].x > emission->data[i-1].recivers[j].signalsArr[k].x)){  // Проверка слева
                    struct mark temp;
                    temp.i = i;
                    temp.j = j;
                    temp.k = k;
                    temp.value = emission->data[i].recivers[j].signalsArr[k].x;
                    marks->push_back(temp);
                }
            }
        }
    }
    return 0;
}


/// \brief - Выделение координат отметок
/// \param [in] emission  - Свертка сигнала
/// \param [in] marks - Массив отметок
/// \param [in, out] targets - Выходной массив целей с характеристиками
int evalCoordinatesMarks(Emission emission, std::vector<struct mark>* marks, std::vector<Target> *targets){
    if(!marks || !targets)
        return 1;

    unsigned int n = 0;
    double lambda = (double)300/60;
    double velocity;
    double azimuth;
    double distance;

    while (n < marks->size()){
        double sum1 = 0, sum2 = 0;
        for(int p = -1; p < 2; p++){
            sum1 += emission.data[marks->at(n).i + p].recivers[marks->at(n).j].signalsArr[marks->at(n).k].x * (marks->at(n).i + p);
            sum2 += emission.data[marks->at(n).i + p].recivers[marks->at(n).j].signalsArr[marks->at(n).k].x;
//            sum1 += data[marks[n].i + p].reciver[marks[n].j].signals[marks[n].k].x * (marks[n].i + p);
//            sum2 += data[marks[n].i + p].reciver[marks[n].j].signals[marks[n].k].x;
//            printf("%d)\t%d\t%d\t%d\t%d\t%lf\n", n, p, marks[n].i, marks[n].j, marks[n].k, data[marks[n].i + p].reciver[marks[n].j].signals[marks[n].k].x);
        }
//        printf("\n");
        velocity = sum1 / sum2;
        targets->at(n).setV(lambda / 2 * velocity *((1.0/durationWave)/countEmission));
//        targets[n].V = ;

        sum1 = sum2 = 0;
        for(int p = -1; p < 2; p++){
            sum1 += emission.data[marks->at(n).i].recivers[marks->at(n).j + p].signalsArr[marks->at(n).k].x * (marks->at(n).j + p);
            sum2 += emission.data[marks->at(n).i].recivers[marks->at(n).j + p].signalsArr[marks->at(n).k].x;
//            sum1 += data[marks[n].i].reciver[marks[n].j + p].signals[marks[n].k].x * (marks[n].j + p);
//            sum2 += data[marks[n].i].reciver[marks[n].j + p].signals[marks[n].k].x;
        }
        azimuth = sum1 / sum2;
        int coef = 1; // коэффициент для сужения выборки, если было расширение при БПФ до 2 ^ n
        double arg = ((azimuth - coef * countRecivers/2)/ coef) * lambda/(distanceAntenna * countRecivers);
        targets->at(n).setB((double)180/M_PI * asin(arg));
//        targets[n].B = (double)180/M_PI * asin(arg);

        sum1 = sum2 = 0;
        for(int p = -1; p < 2; p++){
            sum1 += emission.data[marks->at(n).i].recivers[marks->at(n).j].signalsArr[marks->at(n).k + p].x * (marks->at(n).k + p);
            sum2 += emission.data[marks->at(n).i].recivers[marks->at(n).j].signalsArr[marks->at(n).k + p].x;
//            sum1 += data[marks[n].i].reciver[marks[n].j].signals[marks[n].k + p].x * (marks[n].k + p);
//            sum2 += data[marks[n].i].reciver[marks[n].j].signals[marks[n].k + p].x;
        }
        distance = sum1 / sum2;
        targets->at(n).setR(distance * ((double)150 / frequencyDeviation));
        targets->at(n).setA(marks->at(n).value);

//        targets[n].R = distance * ((double)150 / FrequencyDeviation);
//        targets[n].A = marks[n].value;
        n++;
    }
    return 0;
}


void verticalImitation(std::vector<struct _signal>* array, double elevat){
    double phase = 45;
    double E0 = 50;
    double lambda = (double)300/60;
    double koef = 2 * M_PI / lambda;
    double param = 180./M_PI;

    double polarPhase = 0;
    double polarAmplitude = 0;
    reflectionCoefficient(permittivity, conduct, elevat, lambda, &polarAmplitude, &polarPhase);

    for(int i = 0; i < cntVertAE; i++){
        double arg = - koef * anten[i] * sin(elevat) * phase/param;
        double arg0 = koef * anten[i] * sin(elevat) * phase/param;
        struct _signal temp;
        temp.x = E0 * (cos(arg) + polarAmplitude * cos(arg0 + polarPhase/param));
        temp.y = E0 * (sin(arg) + polarAmplitude * sin(arg0 + polarPhase/param));
        array->push_back(temp);
    }
}

void CFDN(std::vector<struct _signal>* array, double *elevat, double *amplitude){
    double max = 0;
    double lambda = 300./60.;
    for(double i = 0; i < 90; i += 0.1){
        double p = 0;
        double pp = 0;
        for(int j = 0; i < cntVertAE; i++){
            double arg = (2 * M_PI/lambda) * anten[j]*sin(i);
            double co = cos(arg);
            double si = sin(arg);

            p = p - array->at(j).x * co + array->at(j).y * si;
            pp = pp + array->at(j).x * si + array->at(j).y * co;
        }
        double tempAmpl = sqrt(pow(p, 2) + pow(pp, 2));
        if(tempAmpl > max){
            *elevat = i;
            *amplitude = max = tempAmpl;
        }
    }
}

void MX(std::vector<struct _signal>* array, double *elevat, double *amplitude){
    for(double angle = 0; angle < 90; angle += 0.1){
        double tt = 0;
        for(int i = 0; i < cntVertAE; i++){
            double pp = fabs(swh[].x - sw[i].x);
            p = pp * pp;
            tt = tt + p;
        }
        for (int i = 0; i < cntVertAE; i++){
            double pp = fabs(swh[].y - sw[i].y);
            p = pp * pp;
            tt = tt + p;
        }

        if (tt > 0)
            p = 1./tt;
        else
            p = 0;

        if(p > AMX){
            Amx = AMX = p;
            *elevat = angle;
        }
    }
}










