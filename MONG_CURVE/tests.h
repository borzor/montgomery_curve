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
        if(compare(t1,t2)==0){
            std::cout<<"point on curve\n";
        }
        else{
            std::cout<<"point not on curve\n";
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
        if(!compare(point2.X, point4.X)){
            std::cout<<"distributivity is fine\n";
        }
        else{
            std::cout<<"distributivity is not fine\n";
        }
    }
    void qP_(point &point1){
        point point_check(point1);
        cmp.montgomery_ladder(point_check, point1.q);
        if(point_check.X==1 && point_check.Z==0){
            std::cout<<"new point is neutral\n";
        }
    }
    void qP_second(point &point1){
        point point_plus(point1);
        point point_munus(point1);
        NTL::ZZ one(1);
        std::cout<<"default point after transformations from GOST\n"<<point1.X<<'\n'<<point1.Z<<'\n';
        AddMod(point_plus.q,point_plus.q, one, point1.p);
        cmp.montgomery_ladder(point_plus,point_plus.q);
        if (!compare(point_plus.X, point1.X)&&!compare(point_plus.Z, point1.Z)){
            std::cout<<"Point with Q+1 multiplier\nTest passed\n";
        }
        SubMod(point_munus.q,point_munus.q, one, point1.p);
        cmp.montgomery_ladder(point_munus,point_munus.q);
        std::cout<<"Point with q-1 multiplier\n"<<point_munus.X<<'\n'<<point_munus.Z<<'\n';
    }


};
static tests test;
#endif // TESTS_H
