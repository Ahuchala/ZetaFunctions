needsPackage "Divisor";
n = 5;
k = 3;


K = ZZ/101;

pluckerIdeal = Grassmannian(k-1,n-1,CoefficientRing=>K);
R = ring pluckerIdeal;


f = -3/7*p_(0,1,2)^4 + 2*p_(0,1,3)^4 - 7/5*p_(0,2,3)^4 - 1/3*p_(1,2,3)^4 - 5/9*p_(0,1,4)^4 - 5/2*p_(0,2,4)^4 - 7/8*p_(1,2,4)^4 - 9/5*p_(0,3,4)^4 - 3/7*p_(1,3,4)^4 + 7/10*p_(2,3,4)^4

assert isSmooth(pluckerIdeal + ideal f,IsGraded=>true)
