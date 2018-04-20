#ifndef TARGET_H
#define TARGET_H


class Target
{
private:
    double A;  // Амплитуды целей      (dB)
    double R;  // Дальности до целей   (km)
    double B;  // Азимуты целей        (grad)
    double V;  // Скрорость целей      (m/s)
    double F;  // Начальная фаза целей (rad)
    double G;  // Угол места целей     (grad)
    double U;  // Ускорение целей      (m/(s*s))
public:
    Target();
    void setA(double value){ this->A = value; }
    double getA(){ return this->A;}

    void setR(double value){ this->R = value; }
    double getR(){ return this->R;}

    void setB(double value){ this->B = value; }
    double getB(){ return this->B;}

    void setV(double value){ this->V = value; }
    double getV(){ return this->V;}

    void setF(double value){ this->F = value; }
    double getF(){ return this->F;}

    void setG(double value){ this->G = value; }
    double getG(){ return this->G;}

    void setU(double value){ this->U = value; }
    double getU(){ return this->U;}
};

#endif // TARGET_H
