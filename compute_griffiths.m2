-- call with      M2 --script test.m2 8


-- todo: merge these into one function

load "parse_input.m2"

assertSmooth = (n,polynomials,p) -> (
    k = ZZ/p;
    R = k[x_0..x_n]; -- useful for parsing input, not strictly necessary

    polyList = value(polynomials);
    numPoly = length polyList;

    polyList = value(polynomials);

    I = ideal toList polyList;
    J = I + minors(numPoly, jacobian I);
    return saturate J == R;
);

computeGriffithsRing = (n,polynomials,p) -> (
    -- k = QQ;
    k = ZZ/p;
    R = k[x_0..x_n]; -- useful for parsing input, not strictly necessary

    polyList = value(polynomials);
    numPoly = length polyList;
    degs = for b in polyList list first degree b;
    m = sum(degs)- n - 1;


    degreeInput = for i from 0 to n list {0,1};
    degreeInput = degreeInput | for poly in polyList list {1,- first degree poly};

    R = k[x_0..x_n,y_1..y_numPoly,Degrees=>degreeInput];
    polyList = value(polynomials);


    F = sum ( for i from 1 to numPoly list y_i * polyList_(i-1));
    J = R/ideal jacobian F;
    return for i from 0 to (n-numPoly) list toString basis({i,m}, J)

    -- return for i from 0 to (n-numPoly) list toString basis({i,m}, J)
);

computeP1 = (n,polynomials,p) -> (
    k = QQ;
    -- k = ZZ/p;
    S = k[x_0..x_n]; -- useful for parsing input, not strictly necessary

    polyList = value(polynomials);
    numPoly = length polyList;
    degs = for b in polyList list (degree b)#0;
    m = sum(degs)-n - 1;
    
    degreeInput = for i from 0 to n list {0,1};
    degreeInput = degreeInput | for poly in polyList list {1,- first degree poly};

    R = k[x_0..x_n,y_1..y_numPoly,Degrees=>degreeInput];
    return toString basis({1,0}, R)
);

computePn = (n,polynomials,p) -> (
    k = QQ;
    -- k = ZZ/p;
    S = k[x_0..x_n]; -- useful for parsing input, not strictly necessary

    polyList = value(polynomials);
    numPoly = length polyList;
    degs = for b in polyList list (degree b)#0;
    m = sum(degs)-n - 1;


    degreeInput = for i from 0 to n list {0,1};
    degreeInput = degreeInput | for poly in polyList list {1,- first degree poly};

    R = k[x_0..x_n,y_1..y_numPoly,Degrees=>degreeInput];

    return toString basis({n,-n}, R)
);



print parseInput();