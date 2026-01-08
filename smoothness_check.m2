-- Ring with Plücker coordinates for Gr(3,6) ⊂ P^19
-- q = 101;
n = 6;
k = 3;
d = 1;
-- K = ZZ/q -->
K = QQ;

prec = 1;

-- access generators like p_(0,2)

-- R = K[p_(0,0) .. p_(n-1,n-1)]

pluckerIdeal = Grassmannian(k-1,n-1,CoefficientRing=>K);
R = ring pluckerIdeal;


f = p_(0..k-1)-p_(k..2*k-1);
-- Variety X = Gr ∩ V(F)

J = R / (pluckerIdeal + ideal f)
dim singularLocus Proj J