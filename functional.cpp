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

float module(float x, float y){
    float res = sqrtf(powf(x, 2) + powf(y, 2));
    return 20 * log10f(res);
}

int writeDataQt(QString path, Emission *emission){
    QFile file(path);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)){
        qDebug() << "Error opening file";
        return 1;
    }
    QTextStream writeStream(&file);

    for(int i = 0; i < countEmission; i++){
        for(int j = 0; j < countRecivers; j++){
            for(int k = 0; k < countChannels; k++){
                writeStream << k << "\t" << i << "\t"
                        << module(emission->data[i].recivers[j].signalsArr[k].x, emission->data[i].recivers[j].signalsArr[k].y) <<"\n" ;
            }
        }
    }
    file.close();
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
void  noiseGeneration(float *x, float *y){
    assert(x);
    assert(y);
    *x = (float)(rand())/RAND_MAX*(3 - (-3))  -3;
    *y = (float)(rand())/RAND_MAX*(3 - (-3))  -3;
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
void reflectionCoefficient(float permittivity, float conduct, float elevat, float lamda, float *pol_h_x, float *pol_h_y, float *pol_v_x, float *pol_v_y) {
    float angle = elevat * M_PI / 180;
    float sin_angle = sinf(angle);

    float par1 = permittivity - powf(cosf(angle), 2);

    float par2 = -60 * conduct * lamda;
    float tmp = sqrtf(sqrtf(powf(par1, 2) + powf(par2, 2)));

//    float arctan = cos(atan2(par2, par1) / 2.0);
    float arctan = atan2f(par2, par1);
    float px = tmp * cosf(arctan);
    float py = tmp * sinf(arctan);

    float f10x = permittivity * sin_angle - px;
    float f10y = par2 * sin_angle - py;

    float pxx = permittivity * sin_angle + px;
    float pyy = par2 * sin_angle + py;

    // Вертикаль
    float zn = powf(pxx, 2) + powf(pyy, 2);

    float f1x = (f10x * pxx + f10y * pyy)/zn;
    float f1y = (-f10x * pyy + f10y * pxx)/zn;
    *pol_v_x = sqrtf(powf(f1x, 2) + powf(f1y, 2));
    *pol_v_y = - atan2f(f1y, f1x) * 180.0 / M_PI;

    // Горизонталь
    float f20x = sin_angle - px;
    float f20y = - py;

    pxx = sin_angle + px;
    pyy = py;

    zn = powf(pxx, 2) + powf(pyy, 2);

    float f2x = (f20x * pxx + f20y * pyy) / zn;
    float f2y = (-f20x * pyy + f20y * pxx) / zn;

    *pol_h_x = sqrtf(powf(f2x, 2) + powf(f2y, 2));
    *pol_h_y = atan2f(f2y, f2x) * 180.0 / M_PI;
    return;
}



/// \brief - Внутриимпульсная обработка сигнала
/// \param [in, out] *emission - Свертка сигнала
/// \param [in] indexTn - Индекс канала дальности
/// \param [in] indexReciver - Индекс антенного элемента
/// \param [in] windowFun - Оконная функция
void vio(Emission *emission, int indexTn, int indexReciver, int cntChannels, std::vector<float> *windowFun){
    /////////////////////////////////// ДОПИСАТЬ ASSERT
    float arg = 0;
    int kl = (int) (lenProbeSignal + 0.1);

    for (int i = 0; i < cntChannels; i++){
        float c = 0, d = 0;
        for (int j = 0; j < lenProbeSignal; j++){
            int k = i + j;
            if (k >= cntChannels)
                continue;

            arg = (M_PI / kl) * powf((1 + j - (kl + 1) / 2.0), 2);

            c += windowFun->at(j)*  (  emission->data[indexTn].recivers[indexReciver].signalsArr[k].x * cosf(arg)
                                    + emission->data[indexTn].recivers[indexReciver].signalsArr[k].y * sinf(arg));

            d += windowFun->at(j)* (- emission->data[indexTn].recivers[indexReciver].signalsArr[k].x * sinf(arg)
                                    + emission->data[indexTn].recivers[indexReciver].signalsArr[k].y * cosf(arg));
        }
        emission->data[indexTn].recivers[indexReciver].signalsArr[i].x = c;
        emission->data[indexTn].recivers[indexReciver].signalsArr[i].y = d;
    }
}


/// \brief - Имитатор для вертикальной антенной решетки
/// \param [in, out] *emission - Свертка сигнала
/// \param [in] windowFun - Оконная функция
/// \param [in] DF - Частота девиации
/// \param [in] Tn - Период повторения сигнала
/// \param [in] noise - Уровень шума (Sigma)
/// \param [in] targets - Массив целей
void horizontalImitation(Emission *emission, std::vector<float> *windowFun, std::vector<Target> *targets){
    float intervals[countRecivers];
    float heights[countRecivers];

    for(int j = 0; j < countRecivers; j++){
        intervals[j] = distanceAntenna * j;
        heights[j] = heightAntenna;
    }

    float lambda = (float)300/60;
    float dt = (float)1/frequencyDeviation;
    const float piOn180 = (float)M_PI/180;

    for(int i = 0; i < countEmission; i++){
        int st, fin = 0;
        for(int j = 0; j < countRecivers; j++){
            for(int k = 0; k < countChannels; k++){
                float x, y; // X,Y шума по нормальному распределению
                noiseGeneration(&x, &y);
                emission->data[i].recivers[j].signalsArr[k].x = x * noise;
                emission->data[i].recivers[j].signalsArr[k].y = y * noise;

                for(int t = 0; t < countTargets; t++){
                    float a = (targets->at(t).getA() <= -100) ? 0 : powf(10, (float)targets->at(t).getA()/20);
                    float b = targets->at(t).getB();
                    float g = targets->at(t).getG();
                    float r = targets->at(t).getR();
                    float f = targets->at(t).getF();
                    float v = targets->at(t).getV();
                    float u = targets->at(t).getU();


                    float vt = r + (v * i * durationWave) /1000.0 + 0.5 * (u * powf(i,2) * powf(durationWave, 2)) /1000.0;
                    float time = vt / (150 * dt);
                    int begin = (int)floorf(time + 0.49);
                    int end = begin + lenProbeSignal - 1;

                    st = begin; fin = end;

                    if(k >= begin && k <= end && a > 0){
                        float phase = (2 * M_PI/lambda) * (heights[j] + heightSea) * sinf(g * piOn180);
                        float arg = 2 * M_PI/lambda * intervals[j] * sinf(b * piOn180);
                        arg += 2 * M_PI/lambda * 2 * v * i * durationWave;
                        arg += 2 * M_PI/lambda * 2 * v * (k - begin) * dt/1000.0;
                        arg += 2 * M_PI/lambda * (u * powf(i, 2) * powf(durationWave, 2)/2);
                        arg += 2 * M_PI/lambda * (u * powf((k - begin), 2) * powf((dt/1000.0),2));
                        arg += M_PI / lenProbeSignal * powf((0.5 + k - time - (float)lenProbeSignal/2), 2) + f;

                        float polarPhase = 0;
                        float polarAmplitude = 0;
                        float verticalPhase = 0;
                        float verticalAmplitude = 0;
                        reflectionCoefficient(permittivity, conduct, g, lambda, &polarAmplitude, &polarPhase, &verticalAmplitude, &verticalPhase);
                        float phase_x = arg + phase;
                        float phase_y = arg - phase + polarPhase * piOn180;

                        float m1 = cosf(phase_x);
                        m1 += polarAmplitude * cosf(phase_y);
                        float m2 = (sinf(phase_x));
                        m2 += polarAmplitude * sinf(phase_y);

                        emission->data[i].recivers[j].signalsArr[k].x += m1 * windowFun->at(k - begin) * a;
                        emission->data[i].recivers[j].signalsArr[k].y += m2 * windowFun->at(k - begin) * a;
                    }
                }
            }
//            for(int ind = st; ind < (int)(st + (fin - st)/2); ind++){
//                emission->data[i].recivers[j].signalsArr[ind].x = 0.;
//                emission->data[i].recivers[j].signalsArr[ind].y = 0.;
//            }
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
/// \param [in] n - Размерность массива (2^k)
/// \param [in] p - Направление (p > 0 - Прямое, p == 0 - Обратное)
void fft(struct _signal in[], int n, int p){
    int half = n / 2;
    int j = 2, i = 1, k = 1, m, c;
    while (j <= half){
        i = k + half;

        float temp;
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
    float t = M_PI;
    if (p > 0)
        t = -t;

    while (i < n){
        float im = 0,
                re = 1,
                Wpi = sinf(t),
                Wpr = cosf(t),
                temp;
        t *= 0.5;
        int cnt = i;
        i+=i;
        for (int h = 1; h <= cnt; h++){
            for(int j = h; j <= n; j += i){
               k = j + cnt;
               float a = in[k-1].x;
               float b = in[k-1].y;

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
    std::vector<float> windowFun;
    Windowfunction::funcHamming(&windowFun, countRecivers);
    for(int i = 0; i < countEmission; i++){
        for(int k = 0; k < countChannels; k++){
            struct _signal temp[countRecivers];

            for(int j = 0; j < countRecivers; j++){
                temp[j].x = emission->data[i].recivers[j].signalsArr[k].x * windowFun.at(j);
                temp[j].y = emission->data[i].recivers[j].signalsArr[k].y * windowFun.at(j);
            }

            fft(temp, countRecivers, 1);

            for(int j = 0; j < countRecivers; j++){
                emission->data[i].recivers[j].signalsArr[k].x = temp[j].x;
                emission->data[i].recivers[j].signalsArr[k].y = temp[j].y;
            }
        }
    }
    return 0;
}

/// \brief - Доплеровская фильтрация
/// \param [in, out] *emission - Свертка сигнала
/// \param [in, out] countAE - Количество антенных элементов
int dopplerFiltration(Emission* emission, int countAE){
    if(!emission)
        return 1;
    std::vector<float> windowFun;
    Windowfunction::funcHamming(&windowFun, countEmission);

    for(int j = 0; j < countAE; j++){
        for(int k = 0; k < countChannels; k++){
            struct _signal temp[countEmission];

            for(int i = 0; i < countEmission; i++){
                temp[i].x = emission->data[i].recivers[j].signalsArr[k].x * windowFun.at(i);
                temp[i].y = emission->data[i].recivers[j].signalsArr[k].y * windowFun.at(i);
            }

            fft(temp, countEmission, 1);

            for(int i = 0; i < countEmission; i++){
                emission->data[i].recivers[j].signalsArr[k].x = temp[i].x;
                emission->data[i].recivers[j].signalsArr[k].y = temp[i].y;
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
                emission->data[i].recivers[j].signalsArr[k].x = sqrtf(powf(emission->data[i].recivers[j].signalsArr[k].x, 2) +
                                                                     powf(emission->data[i].recivers[j].signalsArr[k].y, 2));
                emission->data[i].recivers[j].signalsArr[k].y = 0;
            }
        }
    }
    return 0;
}

/// \brief - Вычисление порога
/// \param [in] data - Входной массив амплитуд
/// \param [in, out] azimuths  - Массив СКО и СА по каждому азимуту
int calcBorder(Emission *emission, std::vector<struct azimuth> *azimuths){
    if (!azimuths)
        return 1;
    for(int j = 0; j < countRecivers; j++){
        float avg = 0;
        for(int k = 0; k < countChannels; k++){
            for(int i = 0; i < countEmission; i++){
                avg += emission->data[i].recivers[j].signalsArr[k].x;
            }
        }
        struct azimuth temp;
        avg /= (countChannels * countEmission);
        temp.avg = avg;

        float sigma = 0;
        for(int k = 0; k < countChannels; k++){
            for(int i = 0; i < countEmission; i++){
                sigma += powf(emission->data[i].recivers[j].signalsArr[k].x - avg, 2);
            }
        }
        temp.sigma = sqrtf(sigma / (countChannels * countEmission));
        temp.border = alphaBorder * temp.avg + gammaBorder * temp.sigma;
        azimuths->push_back(temp);
    }
    return 0;
}

/// \brief - Выделение координат отметок
/// \param [in, out] emission  - Свертка сигнала
/// \param [in] azimuths  - Массив СКО и СА по каждому азимуту
/// \param [in, out] marks - Выходной массив отметок, значение которого превышает порог border
int marksSelection(Emission *emission, std::vector<struct azimuth> *azimuths, std::vector<struct mark>* marks){
    if (!emission || !azimuths || !marks)
        return 1;
    for(int j = 1; j < countRecivers - 1; j++){
        for(int k = 1; k < countChannels - 1; k++){
            for(int i = 1; i < countEmission - 1; i++){
                if((emission->data[i].recivers[j].signalsArr[k].x > azimuths->at(j).border)  //Проверка с порогом
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
int evalCoordinatesMarks(Emission *emission, std::vector<struct mark>* marks, std::vector<Target> *targets){
    if(!marks || !targets || !emission)
        return 1;

    unsigned int n = 0;
    float lambda = (float)300/60;
    float velocity;
    float azimuth;
    float distance;

    while (n < marks->size()){
        float sum1 = 0, sum2 = 0;
        Target temp;
        int i;
        if(marks->at(n).i < countEmission/2)
            i = marks->at(n).i;
        else
            i = marks->at(n).i - countEmission/2;

        for(int p = -1; p < 2; p++){
            sum1 += emission->data[i + p].recivers[marks->at(n).j].signalsArr[marks->at(n).k].x * (i + p);
            sum2 += emission->data[i + p].recivers[marks->at(n).j].signalsArr[marks->at(n).k].x;
        }
        velocity = sum1 / sum2;
        float v  = lambda / 2 * velocity *(100./countEmission);
        temp.setV(v);

        sum1 = sum2 = 0;
        for(int p = -1; p < 2; p++){
            sum1 += emission->data[i].recivers[marks->at(n).j + p].signalsArr[marks->at(n).k].x * (marks->at(n).j + p);
            sum2 += emission->data[i].recivers[marks->at(n).j + p].signalsArr[marks->at(n).k].x;
        }
        azimuth = sum1 / sum2;
        int coef = 1; // коэффициент для сужения выборки, если было расширение при БПФ до 2 ^ n
        float arg = ((azimuth - coef * countRecivers/2)/ coef) * lambda/(distanceAntenna * countRecivers);
        temp.setB((float)180/M_PI * asinf(arg));

        sum1 = sum2 = 0;
        for(int p = -1; p < 2; p++){
            sum1 += emission->data[i].recivers[marks->at(n).j].signalsArr[marks->at(n).k + p].x * (marks->at(n).k + p);
            sum2 += emission->data[i].recivers[marks->at(n).j].signalsArr[marks->at(n).k + p].x;
        }
        distance = sum1 / sum2;
        temp.setR(distance * ((float)150 / frequencyDeviation));
        temp.setA(marks->at(n).value);
        targets->push_back(temp);
        n++;
    }
    return 0;
}

void generateSignal(std::vector<struct _signal>* array, float elevat){
    float lambda = (float)300./60;
    float koef = 2 * M_PI / lambda;
    float param = 180./M_PI;

    float polarPhase = 0;
    float polarAmplitude = 0;
    float verticalPhase = 0;
    float verticalAmplitude = 0;
    reflectionCoefficient(permittivity, conduct, elevat, lambda, &polarAmplitude, &polarPhase, &verticalAmplitude, &verticalPhase);

    for(int i = 0; i < countVertRecivers; i++){
        float arg = - koef * anten[i] * sinf(elevat/param);
        float arg0 = koef * anten[i] * sinf(elevat/param);
        struct _signal temp;
        temp.x = cosf(arg) + polarAmplitude * cosf(arg0 + polarPhase/param);
        temp.y = sinf(arg) + polarAmplitude * sinf(arg0 + polarPhase/param);
        array->push_back(temp);
    }
}

/// \brief - Имитатор для вертикальной антенной решетки
/// \param [in, out] *emission - Свертка сигнала
/// \param [in] windowFun - Оконная функция
/// \param [in] DF - Частота девиации
/// \param [in] Tn - Период повторения сигнала
/// \param [in] noise - Уровень шума (Sigma)
/// \param [in] targets - Массив целей
void verticalImitation(Emission *emission, std::vector<float> *windowFun, std::vector<Target> *targets){
    float lambda = (float)300/60;   /// ВЫНЕСТИ в const
    float dt = (float)1/frequencyDeviation;
    const float piOn180 = (float)M_PI/180;

    for(int i = 0; i < countEmission; i++){
        for(int j = 0; j < cntAntennas; j++){
            for(int k = 0; k < countChannels; k++){
                float x, y; // X,Y шума по нор мальному распределению
                noiseGeneration(&x, &y);
                emission->data[i].recivers[j].signalsArr[k].x = x * noise;
                emission->data[i].recivers[j].signalsArr[k].y = y * noise;

                for(int t = 0; t < countTargets; t++){
                    float a = (targets->at(t).getA() <= -100) ? 0 : powf(10, (float)targets->at(t).getA()/20);
                    float cos_b = cosf(targets->at(t).getB() * piOn180);
                    float g = targets->at(t).getG();
                    float r = targets->at(t).getR();
                    float f = targets->at(t).getF();
                    float v = targets->at(t).getV();
                    float u = targets->at(t).getU();


                    float vt = r + (v * i * durationWave) /1000.0 + 0.5 * (u * powf(i,2) * powf(durationWave, 2)) /1000.0;
                    float time = vt / (150 * dt);
                    int begin = (int)floor(time + 0.49);
                    int end = begin + lenProbeSignal - 1;

                    if(k >= begin && k <= end && a > 0){
                        float phase = 0.0;
                        if (j >= cntAntennas / 2)
                            phase = 2 * M_PI/lambda * anten[j - cntAntennas / 2] * sinf(g * piOn180);
                        else
                            phase = 2 * M_PI/lambda * anten[j] * sinf(g * piOn180);
//
                        float arg = 2 * M_PI/lambda * 2 * v * i * durationWave;
                        arg += 2 * M_PI/lambda * 2 * v * (k - begin) * dt/1000.0;
                        arg += 2 * M_PI/lambda * (u * powf(i, 2) * powf(durationWave, 2)/2);
                        arg += 2 * M_PI/lambda * (u * powf((k - begin), 2) * powf((dt/1000.0),2));
                        arg += M_PI / lenProbeSignal * powf((0.5 + k - time - (float)lenProbeSignal/2), 2) + f;

                        float gorizontalPhase = 0;
                        float polarAmplitude = 0;
                        float verticalPhase = 0;
                        float verticalAmplitude = 0;
                        reflectionCoefficient(permittivity, conduct, g, lambda, &polarAmplitude, &gorizontalPhase, &verticalAmplitude, &verticalPhase);

                        float m1 = 0;
                        float m2 = 0;
                        if(j >= 8){
//                            double phase_x = arg;
//                            double phase_y = arg + verticalPhase * piOn180;
//                            m1 = cos(phase_x);
//                            m1 += verticalAmplitude * cos(phase_y);
//                            m2 = (sin(phase_x));
//                            m2 += polarAmplitude * sin(phase_y);
                        } else {
                            float phase_x = arg + phase;
                            float phase_y = arg - phase + gorizontalPhase * piOn180;
                            m1 = cosf(phase_x);
                            m1 += polarAmplitude * cosf(phase_y);
                            m2 = (sinf(phase_x));
                            m2 += polarAmplitude * sinf(phase_y);
                        }
                        emission->data[i].recivers[j].signalsArr[k].x += m1 * a * cos_b * windowFun->at(k - begin);
                        emission->data[i].recivers[j].signalsArr[k].y += m2 * a * cos_b * windowFun->at(k - begin);
                    }
                }
            }
        }
        for(int j = 0; j < cntAntennas; j++){
            vio(emission, i, j, countChannels, windowFun);
        }

    }
}


void CFDN(std::vector<struct _signal>* array, float *elevat, float *amplitude){
    float max = 0;
    float lambda = 300./60.;
    for(float i = 0; i < 90; i += 0.1){
        float p = 0;
        float pp = 0;
        for(int j = 0; i < countVertRecivers; i++){
            float arg = (2 * M_PI/lambda) * anten[j] * sinf(i);
            float co = cosf(arg);
            float si = sinf(arg);

            p = p - array->at(j).x * co + array->at(j).y * si;
            pp = pp + array->at(j).x * si + array->at(j).y * co;
        }
        float tempAmpl = sqrtf(powf(p, 2) + powf(pp, 2));
        if(tempAmpl > max){
            *elevat = i;
            *amplitude = max = tempAmpl;
        }
    }
}

void MX(std::vector<struct _signal>* arrayRealSignal, float *elevat ){
    float maxAngle = 0;
    for(float angle = 0; angle < 90; angle += 0.1){
        std::vector<struct _signal> arrayModelSignal;
        generateSignal(&arrayModelSignal, angle);

        //Signal normalization by amplitude
        float normRealSg = 0;
        float normModelSg = 0;
        for(int i = 0; i < countVertRecivers; i++){
            normModelSg += sqrtf(powf(arrayModelSignal.at(i).x,2) + powf(arrayModelSignal.at(i).y,2));
            normRealSg += sqrtf(powf(arrayRealSignal->at(i).x ,2) + powf(arrayRealSignal->at(i).y,2));
        }

        for(int i = 0; i < countVertRecivers; i++){
            arrayModelSignal.at(i).x /= normModelSg;
            arrayModelSignal.at(i).y /= normModelSg;
            arrayRealSignal->at(i).x /= normRealSg;
            arrayRealSignal->at(i).y /= normRealSg;
        }

        //Signal normalization by phase
        float phaseReferenceModel = atan2f(arrayModelSignal.at(1).y, arrayModelSignal.at(1).x);
        float phaseReferenceReal = atan2f(arrayRealSignal->at(1).y, arrayRealSignal->at(1).x);
        float cosPhaseRefModel = cosf(phaseReferenceModel);
        float sinPhaseRefModel = sinf(phaseReferenceModel);
        float cosPhaseRefReal = cosf(phaseReferenceReal);
        float sinPhaseRefReal = sinf(phaseReferenceReal);

        for(int i = 0; i < countVertRecivers; i++){
            float x = arrayModelSignal.at(i).x * cosPhaseRefModel + arrayModelSignal.at(i).y * sinPhaseRefModel;
            float y = - arrayModelSignal.at(i).x * sinPhaseRefModel + arrayModelSignal.at(i).y * cosPhaseRefModel;
            arrayModelSignal.at(i).x = x;
            arrayModelSignal.at(i).y = y;

            x = arrayRealSignal->at(i).x * cosPhaseRefReal + arrayRealSignal->at(i).y * sinPhaseRefReal;
            y = - arrayRealSignal->at(i).x * sinPhaseRefReal + arrayRealSignal->at(i).y * cosPhaseRefReal;
            arrayRealSignal->at(i).x = x;
            arrayRealSignal->at(i).y = y;
        }

        float sumDifference = 0;
        for(int i = 0; i < countVertRecivers; i++){
            float amplitudeRealSignal = sqrtf(powf(arrayRealSignal->at(i).x, 2) + powf(arrayRealSignal->at(i).y, 2));
            float amplitudeModalSignal = sqrtf(powf(arrayModelSignal.at(i).x, 2) + powf(arrayModelSignal.at(i).y, 2));
            float diffAmplitude = powf(fabsf(amplitudeModalSignal - amplitudeRealSignal), 2);
            sumDifference += diffAmplitude;
        }
        float sumX = 0;
        float sumY = 0;
        for(int i = 0; i < countVertRecivers; i++){
            float differenceSignalX = powf(fabsf(arrayModelSignal.at(i).x - arrayRealSignal->at(i).x), 2);
            float differenceSignalY = powf(fabsf(arrayModelSignal.at(i).y - arrayRealSignal->at(i).y), 2);
            sumX += differenceSignalX;
            sumY += differenceSignalY;
        }

        float resultAngle = 0;
        if(sumDifference > 0 && sumX > 0 && sumY > 0)
            resultAngle = 1. / sumDifference + 1. / sumX + 1. / sumY;
        else
            resultAngle = 0;

        if (resultAngle > maxAngle){
            maxAngle = resultAngle;
            *elevat = angle;
        }
    }
}


void calculateElevation(Emission *emission, std::vector<struct mark>* marks, std::vector<Target> *targets){
    if(!targets || !marks || !emission)
        return;

    for(unsigned int i = 0; i < marks->size(); i++){
        std::vector<struct _signal> array;
        for(int j = 0; j < countVertRecivers; j++){
            struct _signal temp = {0, 0};
            temp.x = emission->data[marks->at(i).i].recivers[j].signalsArr[marks->at(i).k].x;
            temp.y = emission->data[marks->at(i).i].recivers[j].signalsArr[marks->at(i).k].y;
            array.push_back(temp);
        }

        for(unsigned int ind = 0; ind < array.size(); ind++){
            qDebug() << ind << "    " <<array.at(ind).x << "    " << array.at(ind).y;
        }

        float elevat;
        MX(&array, &elevat);
        targets->at(i).setG(elevat);
    }
}









