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
    getA(){ return this->A;}

    void setR(double value){ this->R = value; }
    getR(){ return this->R;}

    void setB(double value){ this->B = value; }
    getB(){ return this->B;}

    void setV(double value){ this->V = value; }
    getV(){ return this->V;}

    void setF(double value){ this->F = value; }
    getF(){ return this->F;}

    void setG(double value){ this->G = value; }
    getG(){ return this->G;}

    void setU(double value){ this->U = value; }
    getU(){ return this->U;}

protected:
//    inline QDataStream &operator <<(QDataStream &stream,const Target &tar) //serialization
//    {
////        stream >> tar.setA(stream);
////        stream >> sC.b;
////        stream >> sC.c;
//        return stream;
//    }

};

#endif // TARGET_H
