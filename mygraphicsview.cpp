#include "mygraphicsview.h"
#include <cmath>

MyGraphicsView::MyGraphicsView(QWidget *parent)
    : QGraphicsView(parent){
    this->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    this->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    this->setAlignment(Qt::AlignCenter);

    this->setMinimumHeight(600);
    this->setMinimumWidth(600);

    scene = new QGraphicsScene();
    this->setScene(scene);

    circles = new QGraphicsItemGroup();
    sectors = new QGraphicsItemGroup();

    scene->addItem(circles);
    scene->addItem(sectors);

    timer = new QTimer();

    connect(timer, SIGNAL(timeout()), this, SLOT(redrawTimer()));
    timer->start(50);
}

MyGraphicsView::~MyGraphicsView()
{

}

float function(int value){
    return 0.6 * pow(value, 2) + 100;
}


int rotatePoint_1(QPoint *outPoint, QPoint *inPoint, QPoint *center, int angle){
    if (!outPoint){
        return -1;
    }
    outPoint->setX((int)(center->x() + (inPoint->x() - center->x()) * cos(angle * M_PI/180) - (inPoint->y() - center->y()) * sin(angle*M_PI/180)));
    outPoint->setY((int)(center->y() + (inPoint->y() - center->y()) * cos(angle * M_PI/180) + (inPoint->x() - center->y()) * sin(angle*M_PI/180)));
    return 0;
}

int rotatePoint(QPoint *outPoint, QPoint inPoint, QPoint center, int angle){
    if (!outPoint){
        return -1;
    }
    outPoint->setX((int)(center.x() + (inPoint.x() - center.x()) * cos(angle * M_PI/180) - (inPoint.y() - center.y()) * sin(angle*M_PI/180)));
    outPoint->setY((int)(center.y() + (inPoint.y() - center.y()) * cos(angle * M_PI/180) + (inPoint.x() - center.y()) * sin(angle*M_PI/180)));
    return 0;
}

void MyGraphicsView::redrawTimer()
{
    this->deleteItemsFromGroup(circles);
    this->deleteItemsFromGroup(sectors);

    int width = this->width();
    int height = this->height();
    scene->setSceneRect(0,0,width,height);

    int centerX = this->width() / 2;
    int centerY = this->height() / 2;

    QBrush redBrush(Qt::red);
    QBrush greenBrush(Qt::green);
    QPen blackPen(Qt::black);

    circles->addToGroup(scene->addEllipse(centerX-4, centerY-4, 8, 8, blackPen, greenBrush));

    int sideOfSquare = (height > width) ? (width-40) : (height-40);
    for(int i = 0; i < 400; i+=100){
        circles->addToGroup(scene->addEllipse(centerX - (sideOfSquare/2) + i, centerY - (sideOfSquare/2) + i, sideOfSquare - 2*i, sideOfSquare - 2*i, blackPen));
    }

//    scene->setSceneRect(20, 20, sideOfSquare, sideOfSquare);
//    centerX = scene->width()/2;
//    centerY = scene->height()/2;

    QPoint point_1(centerX - (sideOfSquare/2), centerY);
    QPoint point_2(centerX + (sideOfSquare/2), centerY);
    QPoint center(centerX, centerY);

    for(int i = 0; i < 6; i++){
        sectors->addToGroup(scene->addLine((int)(point_1.x()), (int)(point_1.y()), (int)(point_2.x()), (int)(point_2.y()), blackPen));
        rotatePoint(&point_1, point_1, center, 30);
        rotatePoint(&point_2, point_2, center, 30);
    }

//    int sideOfSquare = (height > width) ? (width - 60) : (height - 60);
//    int centerOfWidget_X = width/2;
//    int centerOfWidget_Y = height/2;
//    group_2->addToGroup(scene->addLine(centerOfWidget_X - (sideOfSquare/2),
//                                       centerOfWidget_Y - (sideOfSquare/2),
//                                       centerOfWidget_X + (sideOfSquare/2),
//                                       centerOfWidget_Y - (sideOfSquare/2),
//                                       penRed));

//    group_2->addToGroup(scene->addLine(centerOfWidget_X + (sideOfSquare/2),
//                                       centerOfWidget_Y - (sideOfSquare/2),
//                                       centerOfWidget_X + (sideOfSquare/2),
//                                       centerOfWidget_Y + (sideOfSquare/2),
//                                       penRed));

//    group_2->addToGroup(scene->addLine(centerOfWidget_X + (sideOfSquare/2),
//                                       centerOfWidget_Y + (sideOfSquare/2),
//                                       centerOfWidget_X - (sideOfSquare/2),
//                                       centerOfWidget_Y + (sideOfSquare/2),
//                                       penRed));

//    group_2->addToGroup(scene->addLine(centerOfWidget_X - (sideOfSquare/2),
//                                       centerOfWidget_Y + (sideOfSquare/2),
//                                       centerOfWidget_X - (sideOfSquare/2),
//                                       centerOfWidget_Y - (sideOfSquare/2),
//                                       penRed));

}


void MyGraphicsView::resizeEvent(QResizeEvent *event)
{
    timer->start(50);
    QGraphicsView::resizeEvent(event);
}


void MyGraphicsView::deleteItemsFromGroup(QGraphicsItemGroup *group)
{
    foreach( QGraphicsItem *item, scene->items(group->boundingRect())) {
       if(item->group() == group ) {
          delete item;
       }
    }
}
