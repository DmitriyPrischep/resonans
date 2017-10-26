#include "emission.h"

//struct radio allocateEmission(int size)

Emission::Emission()
{
//    this->data = allocateEmission(int size);
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
    delete[] this->data;
}
