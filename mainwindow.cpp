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

//        item = new QStandardItem(QString(targets->at(i).getG()));
//        model->setItem(i, 5, item);

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
    QSpacerItem* horizontalSpacer = new QSpacerItem(300, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    QGridLayout* layout = (QGridLayout*)msgError.layout();
    layout->addItem(horizontalSpacer, layout->rowCount(), 0, 1, layout->columnCount());
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

    imitationTargets(emis, &windowFunc, &targets);

    CFDN(emis);
    dopplerFiltration(emis);
    detection(emis);

    std::vector<struct azimuth> azimuths;

    calcBorder(emis, &azimuths);

    std::vector<struct mark> marks;

    marksSelection(emis, &azimuths, &marks);

    qDebug() << "  Mark size is " + marks.size();

    std::vector<Target> resultTargets;
    evalCoordinatesMarks(emis, &marks, &resultTargets);

    int a;
    a = resultTargets.size();

    initializationTable(&resultTargets);


//    QSettings *settings = new QSettings("configure.conf", QSettings::IniFormat);
//    Settings::setValues(settings);

    ////////////////////////////////////////

    delete[] emis->data;
}

void MainWindow::on_pushButton_2_clicked()
{
    picture = new MyGraphicsView();
//    //Движение цели
//    int startX = -20;
//    int startY = -20;
//    for(int i = 0, j = 0; i < 100 && j < 100; i++, j+=3){
//        scene->addEllipse(startX+i, function(startY+i), 4, 4, blackPen, redBrush);
//    }
}
