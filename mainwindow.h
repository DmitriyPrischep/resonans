#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QWidget>
#include <QPainter>
#include <QtCore>
#include <QtGui>
#include "mygraphicsview.h"
#include "target.h"


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void setPath(QString str){this->pathFile = str;}
    QString getPath(){return this->pathFile;}

private slots:
    void on_btnFileDialog_clicked();

    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

private:
    void initializationTable(std::vector<Target> *targets);
    void viewMessageBox(QString title, QString text);
    Ui::MainWindow *ui;
    MyGraphicsView *picture;
//    QGraphicsScene *scene;
//    QGraphicsEllipseItem *ellipse;
//    QGraphicsRectItem *rectangle;

    QString pathFile;

};

#endif // MAINWINDOW_H
