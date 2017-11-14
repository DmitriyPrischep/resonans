#include "emission.h"
#include <QDebug>

Emission::Emission()
{

}

struct radioSignal* Emission::allocateEmission(int size){
    radioSignal* temp = new radioSignal[size];
    if(!temp)
        return NULL;
    else
        return temp;
}

Emission::Emission(int sizeEmission){
    this->data = allocateEmission(sizeEmission);
}

Emission::~Emission(){
    qDebug() << "Destructor worked\n" ;
    delete[] this->data;
}
