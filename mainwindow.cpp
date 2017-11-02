#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "target.h"
#include "emission.h"
#include "functional.h"
#include "windowfunction.h"
#include <QDebug>
#include <QFileDialog>
#include <QMessageBox>
#include <ctime>

float function(int value){
    return 0.6 * pow(value, 2) + 100;
}


int rotatePoint(int* outX, int* outY, int inX, int inY, int centerX, int centerY, int angle){
    if (!outX || !outY){
        return -1;
    }
    *outX = (int)(centerX + (inX - centerX) * cos(angle*M_PI/180) - (inY - centerY) * sin(angle*M_PI/180));
    *outY = (int)(centerY + (inY - centerY) * cos(angle*M_PI/180) + (inX - centerX) * sin(angle*M_PI/180));
    return 0;
}

//X = x0 + (x - x0) * cos(a) - (y - y0) * sin(a);
//Y = y0 + (y - y0) * cos(a) + (x - x0) * sin(a);

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{    
    ui->setupUi(this);


    scene = new QGraphicsScene(this);
    ui->graphicsView->setScene(scene);

//    QBrush redBrush(Qt::red);
//    QBrush greenBrush(Qt::green);
//    QPen blackPen(Qt::black);
////    blackPen.setWidth(6);

//    int centerX = scene->width() / 2;
//    int centerY = scene->height() / 2;

//    scene->addEllipse(centerX- 350, centerY - 350, 700, 700, blackPen);
//    scene->addEllipse(centerX - 250, centerY - 250, 500, 500, blackPen);
//    scene->addEllipse(centerX - 150, centerY - 150, 300, 300, blackPen);

//    int x1 = 350;
//    int y1 = 0;
//    int x2 = -350;
//    int y2 = 0;
//    scene->addLine(x1, y1, x2, y2);
//    for(int i = 0; i < 5; i++){
//        rotatePoint(&x1, &y1, x1, y1, centerX, centerY, 30);
//        rotatePoint(&x2, &y2, x2, y2, centerX, centerY, 30);
//        scene->addLine(x1, y1, x2, y2);
//    }

//    scene->addEllipse(centerX - 4, centerY - 4, 8, 8, blackPen, greenBrush);

//    //Движение цели
//    int startX = -20;
//    int startY = -20;
//    for(int i = 0, j = 0; i < 100 && j < 100; i++, j+=3){
//        scene->addEllipse(startX+i, function(startY+i), 4, 4, blackPen, redBrush);
//    }


}


MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_btnFileDialog_clicked()
{
    setPath(QFileDialog::getOpenFileName(this, "Open file targets", "", "*.txt"));
    ui->editPath->setText(getPath());
}

void MainWindow::on_pushButton_clicked()
{
    srand(time(NULL));

    std::vector<Target> targets;
    int code = readDataQt(this->getPath(), &targets);
    if(code != 0){
        QMessageBox msgError;
         if(code == 2){
             msgError.setText(tr("Error open file"));
             msgError.setInformativeText(tr("File does not exist"));
             msgError.setStandardButtons(QMessageBox::Ok);
             msgError.setDefaultButton(QMessageBox::Ok);
             msgError.exec();
             return;
         }
         msgError.setText(tr("Error open file"));
         msgError.setInformativeText(tr("Error reading file"));
         msgError.setStandardButtons(QMessageBox::Ok);
         msgError.setDefaultButton(QMessageBox::Ok);
         msgError.exec();
         return;
    }

    Emission emis;
    emis.data[120].recivers[9].signalsArr[200].x = 99;
    emis.data[120].recivers[9].signalsArr[200].y = 33;
    ui->label->setText(QString::number(emis.data[120].recivers[9].signalsArr[200].x) +  QString::number(emis.data[120].recivers[9].signalsArr[200].y));

    std::vector<double> windowFunc;
    Windowfunction::funcHammingTukey(&windowFunc);


    int a;
    a = 5;

}
