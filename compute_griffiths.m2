-- call with      M2 --script test.m2 8


-- print value(scriptCommandLine#1)

-- sage code is like run_macualay2_program(function_name, args)

-- polyList is a list [f_1,...,f_c] of polynomials

computeGriffithsRing = (n,polynomials) -> (
k = QQ;
S = k[x_0..x_n]; -- useful for parsing input, not strictly necessary

polyList = value(polynomials);
numPoly = length polyList;
degs = for b in polyList list (degree b)#0;
m = n + 1 - sum(degs);


degreeInput = for i from 0 to n list {0,1};
degreeInput = degreeInput | for poly in polyList list {1,- first degree poly};

R = k[x_0..x_n,y_1..y_numPoly,Degrees=>degreeInput];
polyList = value(polynomials);


F = sum ( for i from 1 to numPoly list y_i * polyList_(i-1));
J = R/ideal jacobian F;
return for i from 0 to (n-numPoly) list toString basis({i,m}, J)
);


print computeGriffithsRing(2,"[x_0^2]");

--function_name = value(scriptCommandLine#1)

--capture function_name

