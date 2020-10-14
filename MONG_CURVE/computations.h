#ifndef COMPUTATIONS_H
#define COMPUTATIONS_H
#include "mong_params.h"
class computations
{
public:
    computations()//http://hyperelliptic.org/EFD/g1p/auto-montgom-xz.html
    {}
    void double_point(point &point1){ //dbl-1987-m-3
        NTL::ZZ a24=point1.a,two(2),A,B,C,four(4);
        AddMod(A,point1.X,point1.Z,point1.p);
        PowerMod(A,A,two,point1.p);
        SubMod(B,point1.X,point1.Z,point1.p);
        PowerMod(B,B,two,point1.p);
        SubMod(C,A,B,point1.p);// (x+z)^2-(x-z)^2
        MulMod(point1.X,A,B,point1.p);// (x+z)^2*(x-z)^2
        InvMod(four,four,point1.p);// 1/4
        AddMod(a24,a24,two,point1.p);//a-2
        MulMod(a24,a24,four,point1.p);//a-2/4
        MulMod(point1.Z,a24,C,point1.p);//(a-2/4)*(x+z)^2-(x-z)^2
        AddMod(point1.Z,point1.Z,B,point1.p);//((a-2/4)*(x+z)^2-(x-z)^2)+(x-z)^2
        MulMod(point1.Z,point1.Z,C,point1.p);
    }
    void add_point(point &point3, point &point2, point &default_){// dadd-1987-m-2
        NTL::ZZ two(2), t1, t2, t3, t4;
        SubMod(t1,point3.X, point3.Z, point2.p);// X3-Z3
        AddMod(t2, point3.X, point3.Z, point2.p);//X3+Z3
        AddMod(t3, point2.X, point2.Z, point2.p);//X2+Z2
        SubMod(t4, point2.X, point2.Z, point2.p);//X2-Z2
        MulMod(t1, t1, t3, point2.p);//(X3-Z3)*(X2+Z2)
        MulMod(t2, t2, t4, point2.p);//(X3+Z3)*(X2-Z2)
        AddMod(point3.X, t1, t2, point2.p);//(X3-Z3)*(X2+Z2)+(X3+Z3)*(X2-Z2)
        PowerMod(point3.X, point3.X, two, point2.p);//((X3-Z3)*(X2+Z2)+(X3+Z3)*(X2-Z2))^2
        MulMod(point3.X, point3.X, default_.Z, point2.p);//Z0*((X3-Z3)*(X2+Z2)+(X3+Z3)*(X2-Z2))^2
        SubMod(point3.Z, t1, t2, point2.p);//(X3-Z3)*(X1+Z1)-(X1-Z1)*(X2+Z2)
        PowerMod(point3.Z, point3.Z, two, point2.p);
        MulMod(point3.Z, point3.Z, default_.X, point2.p);
    }
    void point_transform(point &point1){
        NTL::ZZ invert;
        if(point1.Z!=0){
            InvMod(invert, point1.Z, point1.p);
            MulMod(point1.X,point1.X,invert,point1.p);
            MulMod(point1.Z,point1.Z,invert,point1.p);
        }
        else{
            InvMod(invert,point1.X,point1.p);
            MulMod(point1.X,point1.X,invert,point1.p);
        }
    }
    void montgomery_ladder(point &point_old, NTL::ZZ q){ //https://drive.google.com/drive/u/1/folders/17_F1NM91KR-6HOnUG_pHWxmtlRf7auU1 slide 54-55
        long size = NumBits(q);
        point point_new(point_old);
        point point_zero;
        point_zero.X = 1;
        point_zero.Z = 0;
        for(long i = size-1; i>=0; i-- ){
            if(NTL::bit(q,i)){
                add_point(point_zero, point_new, point_old);
                double_point(point_new);
            }
            else{
                add_point(point_new, point_zero, point_old);
                double_point(point_zero);
            }
        }
        point_old.X = point_zero.X;
        point_old.Z = point_zero.Z;
        point_transform(point_old);
    }

};

static computations cmp;
#endif // COMPUTATIONS_H
