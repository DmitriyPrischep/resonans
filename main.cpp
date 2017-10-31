#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QCoreApplication::setOrganizationName("Research Institute Rezonans");
    QCoreApplication::setApplicationName("Radar simulator");

    QApplication a(argc, argv);
    QTranslator traslator;
    traslator.load(":/languages/program_" + QLocale::system().name());
    a.installTranslator(&traslator);

    MainWindow w;
    w.show();

    return a.exec();
}
