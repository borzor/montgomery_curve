#include <iostream>
#include "tests.h"
int main()
{
point point1;
NTL::ZZ k1,k2;
RandomBnd(k1,point1.p);
RandomBnd(k2,point1.p);
test.qP_second(point1);
test.qP_(point1);
test.is_point_on_curve(point1);
test.distributivity(point1, k1, k2);
return 0;
}
