#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "target.h"
#include "emission.h"
#include "functional.h"
#include "windowfunction.h"
#include <QDebug>
#include <QFileDialog>
#include <QMessageBox>
#include <QStandardItemModel>
#include <QStandardItem>
#include <ctime>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{    
    ui->setupUi(this);
    srand(time(NULL));
    this->setWindowTitle(QCoreApplication::applicationName());
    ui->lblPath->setText(trUtf8("File targets:"));
}


MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_btnFileDialog_clicked()
{
    setPath(QFileDialog::getOpenFileName(this, trUtf8("Open file targets"), "", "*.txt"));
    ui->editPath->setText(getPath());
}

void MainWindow::initializationTable(std::vector<Target> *targets){
    QStandardItemModel *model = new QStandardItemModel;
    QStandardItem *item;

    QStringList columnHeader;
    columnHeader.append(trUtf8("Amplitude"));
    columnHeader.append(trUtf8("Distanse"));
    columnHeader.append(trUtf8("Azimuth"));
    columnHeader.append(trUtf8("Speed"));
    columnHeader.append(trUtf8("Start phase"));
    columnHeader.append(trUtf8("Elevation"));
    columnHeader.append(trUtf8("Acceleration"));
    model->setHorizontalHeaderLabels(columnHeader);

    for(unsigned int i = 0; i < targets->size(); i++){
        item = new QStandardItem(QString::number(targets->at(i).getA()));
        model->setItem(i, 0, item);

        item = new QStandardItem(QString::number(targets->at(i).getR()));
        model->setItem(i, 1, item);

        item = new QStandardItem(QString::number(targets->at(i).getB()));
        model->setItem(i, 2, item);

        item = new QStandardItem(QString::number(targets->at(i).getV()));
        model->setItem(i, 3, item);

//        item = new QStandardItem(QString(targets->at(i).getF()));
//        model->setItem(i, 4, item);

        item = new QStandardItem(QString::number(targets->at(i).getG(), 'f', 2));
        model->setItem(i, 5, item);

//        item = new QStandardItem(QString(targets->at(i).getU()));
//        model->setItem(i, 6, item);
    }

    ui->tableView->setModel(model);

    ui->tableView->resizeRowsToContents();
    ui->tableView->resizeColumnsToContents();
}

void MainWindow::viewMessageBox(QString title, QString text){
    QMessageBox msgError;
    msgError.setWindowTitle(title);
    msgError.setInformativeText(text);
    msgError.setStandardButtons(QMessageBox::Ok);
    msgError.setDefaultButton(QMessageBox::Ok);
    //QSpacerItem* horizontalSpacer = new QSpacerItem(300, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    //QGridLayout* layout = (QGridLayout*)msgError.layout();
    //layout->addItem(horizontalSpacer, layout->rowCount(), 0, 1, layout->columnCount());
    msgError.exec();
}

void MainWindow::on_pushButton_clicked()
{
    std::vector<Target> targets;
    int code = readDataQt(this->getPath(), &targets);
    if(code != 0){
         if(code == 2){
             viewMessageBox(trUtf8("Error open file"), trUtf8("File does not exist"));
             return;
         }
         viewMessageBox(trUtf8("Error open file"), trUtf8("Error reading file"));
         return;
    }

    Emission* emis = new Emission(countEmission);

    std::vector<double> windowFunc;
    Windowfunction::funcHammingTukey(&windowFunc);

    horizontalImitation(emis, &windowFunc, &targets);

    CFDN(emis);
    dopplerFiltration(emis, countRecivers);
    detection(emis);

    std::vector<struct azimuth> azimuths;

    calcBorder(emis, &azimuths);

    std::vector<struct mark> marks;

    marksSelection(emis, &azimuths, &marks);

    qDebug() << "  Mark size is " + marks.size();

    std::vector<Target> resultTargets;
    evalCoordinatesMarks(emis, &marks, &resultTargets);

    Emission* emisVertical = new Emission(countEmission);
    verticalImitation(emisVertical, &windowFunc, &targets);

    dopplerFiltration(emisVertical, countVertRecivers);

    calculateElevation(emisVertical, &marks, &resultTargets);

//    for(int i = 0; i < countTargets; i++){
//        std::vector<struct _signal> signalTarget;
//        verticalImitation(&signalTarget, targets.at(i).getG());
//    }

    initializationTable(&resultTargets);

//    QSettings *settings = new QSettings("configure.conf", QSettings::IniFormat);
//    Settings::setValues(settings);
    delete[] emis->data;
    delete[] emisVertical->data;
}

void MainWindow::on_pushButton_2_clicked()
{

    Emission* emisVertical = new Emission(countEmission);
    std::vector<Target> targets;
    int code = readDataQt("targets.txt", &targets);
    std::vector<double> windowFunc;
    Windowfunction::funcHammingTukey(&windowFunc);
    verticalImitation(emisVertical, &windowFunc, &targets);

//    Emission* emis = new Emission(countEmission);
//    std::vector<Target> targets;
//    int code = readDataQt("targets.txt", &targets);
//    std::vector<double> windowFunc;
//    Windowfunction::funcHammingTukey(&windowFunc);
//    horizontalImitation(emis, &windowFunc, &targets);
    double angle = 1.80;

    writeDataQt("data_graphics.txt", emisVertical);




//    double a[] = {2.380, 4.532, 5.967, 6.763, 6.302, 4.444, 1.652, -1.659};
//    double b[] = {-1.168, -2.292, -3.242, -3.647, -3.655, -2.595, -1.113, 0.643};

//    float a[] = {-0.51, -0.82, -0.94, -0.79, -0.31, 0.38, 0.79, 0.91};
//    float b[] = {0.82, 1.49, 1.75, 1.58, 0.56, -0.58, -1.51, -1.71};
    /* TEST ON MX()
    //Угол 1.8 град
    float a[] = {0.18, 0.34, 0.52, 0.69, 0.89, 1.07, 1.26, 1.42};
    float b[] = {-0.10, -0.19, -0.32, -0.38, -0.54, -0.65, -0.74, -0.86};
    std::vector<struct _signal> array;
    for(int j = 0; j < countVertRecivers; j++){
        struct _signal temp = {0, 0};
        temp.x = a[j];
        temp.y = b[j];
        array.push_back(temp);
    }
    double e;
    MX(&array, &e);
    e += 1 - 1;
    */

//    picture = new MyGraphicsView();
//    //Движение цели
//    int startX = -20;
//    int startY = -20;
//    for(int i = 0, j = 0; i < 100 && j < 100; i++, j+=3){
//        scene->addEllipse(startX+i, function(startY+i), 4, 4, blackPen, redBrush);
//    }
}
