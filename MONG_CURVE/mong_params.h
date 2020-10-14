#ifndef MONG_PARAMS_H
#define MONG_PARAMS_H
#include <NTL/ZZ.h>
#include "parameters.h"

class point{
public:
    NTL::ZZ X,Y,Z,p,q,a,b;
    point(){
        NTL::ZZ v,u,e,d;
        NTL::ZZ one(1),two(2),four(4);
        Z=one;
        conv(v,V);
        conv(u,U);
        conv(e,E);
        conv(d,D);
        conv(p, MODULE);
        conv(q,Q);
        AddMod(a,e,d,p); // a + d
        MulMod(a,a,two,p); // 2(a+d)
        SubMod(b,e,d,p);// a - d
        InvMod(b,b,p);// (a - d)^(-1)
        MulMod(a,a,b,p);// 2(a+d)/(a-d)
        MulMod(b,b,four,p);// 4/(a-d)
        AddMod(X,one,v,p);// 1 + v
        SubMod(Y,one,v,p);// 1 - v
        InvMod(Y,Y,p);// (1 - v)^(-1)
        MulMod(X,X,Y,p);//(1 + v)/(1 - v)
        InvMod(Y,u,p);// u^(-1)
        MulMod(Y,X,Y,p);// (1+v)/((1-v)*u)
    }
    point(const point &p2):
           X(p2.X),Z(p2.Z),p(p2.p),q(p2.q),a(p2.a),b(p2.b)
    {}
};



#endif
