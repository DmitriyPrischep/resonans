#ifndef TARGET_H
#define TARGET_H


class Target
{
private:
    float A;  // Амплитуды целей      (dB)
    float R;  // Дальности до целей   (km)
    float B;  // Азимуты целей        (grad)
    float V;  // Скрорость целей      (m/s)
    float F;  // Начальная фаза целей (rad)
    float G;  // Угол места целей     (grad)
    float U;  // Ускорение целей      (m/(s*s))
public:
    Target();
    void setA(float value){ this->A = value; }
    float getA(){ return this->A;}

    void setR(float value){ this->R = value; }
    float getR(){ return this->R;}

    void setB(float value){ this->B = value; }
    float getB(){ return this->B;}

    void setV(float value){ this->V = value; }
    float getV(){ return this->V;}

    void setF(float value){ this->F = value; }
    float getF(){ return this->F;}

    void setG(float value){ this->G = value; }
    float getG(){ return this->G;}

    void setU(float value){ this->U = value; }
    float getU(){ return this->U;}
};

#endif // TARGET_H
