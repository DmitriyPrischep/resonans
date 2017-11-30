#ifndef MYGRAPHICSVIEW_H
#define MYGRAPHICSVIEW_H

#include <QWidget>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QGraphicsItemGroup>
#include <QTimer>

class MyGraphicsView : public QGraphicsView
{
    Q_OBJECT
public:
    explicit MyGraphicsView(QWidget *parent = 0);
    ~MyGraphicsView();

private slots:
    void redrawTimer();

private:
    void resizeEvent(QResizeEvent *event);
    void deleteItemsFromGroup(QGraphicsItemGroup *group_1s);

    QGraphicsScene *scene;
    QGraphicsItemGroup *circles;
    QGraphicsItemGroup *sectors;

    QTimer *timer;
};

#endif // MYGRAPHICSVIEW_H
