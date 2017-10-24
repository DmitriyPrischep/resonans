#include <iostream>
#include <fstream>
#include <QFile>
#include <QDebug>
#include <QDataStream>
#include "target.h"


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

int readDataQt(std::vector<Target> *targets){
    QFile file("targets1.txt");
    if(!file.exists()){
        qDebug() << "File does not exist";
        return 2;
    }
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)){
        qDebug() << "Error reading file";
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
