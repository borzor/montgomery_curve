#ifndef TESTS_H
#define TESTS_H
#include "computations.h"
class tests
{
public:
    void is_point_on_curve(point &point1){
        NTL::ZZ t1, t2, t3, two(2), three(3);
        PowerMod(t1, point1.Y, two, point1.p);
        MulMod(t1, t1, point1.b, point1.p);
        PowerMod(t2, point1.X, three, point1.p);
        PowerMod(t3, point1.X, two, point1.p);
        MulMod(t3, t3, point1.a, point1.p);
        AddMod(t2, t2, t3, point1.p);
        AddMod(t2, t2, point1.X, point1.p);
        std::cout<<"checking if point on curve\nX = "<<point1.X<<"\nZ = "<<point1.Z<<'\n';
        if(compare(t1,t2)==0){
            std::cout<<"point on curve\n\n";
        }
        else{
            std::cout<<"point not on curve\n\n";
        }
    }
    void distributivity(point &point1, NTL::ZZ k1, NTL::ZZ k2){
        point point2(point1);
        point point3(point1);
        point point4(point1);
        point point5(point1);
        if(k1<k2){
            swap(k1,k2);
        }
        cmp.montgomery_ladder(point2, k1);
        cmp.montgomery_ladder(point3, k2);
        cmp.montgomery_ladder(point4, k1+k2);
        cmp.montgomery_ladder(point5, k1-k2);
        cmp.add_point(point2,point3,point5);
        cmp.point_transform(point2);
        std::cout<<"checking K1*P + K2*P = (K1+K2)*P\n";
        if(!compare(point2.X, point4.X)){
            std::cout<<"Test passed\n";
        }
        else{
            std::cout<<"Test failed\n";
        }
    }
    void qP_(point &point1){
        point point_check(point1);
        cmp.montgomery_ladder(point_check, point1.q);
        std::cout<<"Point with Q multiplier\n";
        if(point_check.X==1 && point_check.Z==0){
            std::cout<<"new point is neutral\n\n";
        }
    }
    void qP_second(point &point1){
        point point_plus(point1);
        point point_minus(point1);
        NTL::ZZ one(1);
        std::cout<<"default point after transformations from GOST\nX = "<<point1.X<<"\nZ = "<<point1.Z<<"\n\n";
        AddMod(point_plus.q,point_plus.q, one, point1.p);
        cmp.montgomery_ladder(point_plus,point_plus.q);
        if (!compare(point_plus.X, point1.X)&&!compare(point_plus.Z, point1.Z)){
            std::cout<<"Point with Q+1 multiplier\nX = "<<point_minus.X<<"\nZ = "<<point_minus.Z<<'\n'<<"Test passed\n";
        }
        SubMod(point_minus.q,point_minus.q, one, point1.p);
        cmp.montgomery_ladder(point_minus,point_minus.q);
        if (!compare(point_minus.X, point1.X)&&!compare(point_minus.Z, point1.Z)){
            std::cout<<"Point with q-1 multiplier\nX = "<<point_minus.X<<"\nZ = "<<point_minus.Z<<'\n'<<"Test passed\n\n";
        }
    }


};
static tests test;
#endif // TESTS_H
