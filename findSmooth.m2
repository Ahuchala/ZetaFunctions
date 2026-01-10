n = 6;
k = 2;
K = QQ;
d = 2;
pluckerIdeal = Grassmannian(k-1,n-1,CoefficientRing=>K);
R = ring pluckerIdeal;
-- f = random(d,R);
f = sum(apply(gens R, i->random(ZZ)*(i)^d))

countInversions = (perm) -> (
    ans := 0;
    for i from 0 to #perm-2 do (
        for j from i+1 to #perm-1 do (
            if perm#j < perm#i then ans = ans + 1;
        );
    );
    ans
);

-- Build lookup table once
varLookup := hashTable apply(gens R, v -> (
    toList (baseName v)#1 => v
));

applyPermutationToPoly = (f, perm) -> (
    subList := apply(gens R, v -> (
        indices := toList (baseName v)#1;
        newIndices := apply(indices, i -> perm#i);
        sortedIndices := sort newIndices;
        sign := (-1)^(countInversions(newIndices));
        v => sign * varLookup#sortedIndices
    ));
    sub(f, subList)
);

sumOverPermutations = (f, permList) -> (
    sum apply(permList, perm -> applyPermutationToPoly(f, perm))
);

-- permList = {{0,1,2,3,4}, {1,0,2,3,4}};
permList = {{0,1,2,3,4,5},{1,0,5,4,3,2},{2,4,0,5,1,3},{3,5,4,0,2,1},{4,2,3,1,5,0},{5,3,1,2,0,4}};

symF = sumOverPermutations(f, permList);
print symF;
assert isSmooth(variety ideal symF);


-- {{p_(0,1)*p_(3,4)^2, p_(0,2)*p_(3,4)^2, p_(1,2)*p_(3,4)^2,p_(0,3)*p_(3,4)^2, p_(1,3)*p_(3,4)^2, p_(2,3)*p_(3,4)^2, p_(0,4)*p_(3,4)^2,p_(1,4)*p_(3,4)^2, p_(2,4)*p_(3,4)^2, p_(3,4)^3}}

-- Function to get the matrix representation of a permutation acting on a basis in J
-- permutationMatrixInJ = (perm, basisElems, J) -> (
--     -- basisElems should be a list of elements in J
--     m := #basisElems;
    
--     -- Lift to R, apply permutation, reduce back to J
--     images := apply(basisElems, b -> (
--         bLifted := lift(b, R);
--         imageInR := applyPermutationToPoly(bLifted, perm);
--         sub(imageInR, J)
--     ));
    
--     -- Express each image as a linear combination of basis elements
--     coeffMatrix := {};
--     for img in images do (
--         imgLifted := lift(img, R);
--         coeffs := apply(basisElems, b -> (
--             bLifted := lift(b, R);
--             coefficient(bLifted, imgLifted)
--         ));
--         coeffMatrix = append(coeffMatrix, coeffs);
--     );
    
--     transpose matrix coeffMatrix
-- );

-- -- Get basis in J at a given degree
-- getBasisInJ = (deg, J) -> (
--     B := basis(deg, J);
--     if B == 0 then return {};
--     flatten entries B
-- );

-- -- Example usage with your basis
-- myBasis = flatten entries basis(3, J);  -- degree 2d

-- perm = {1,0,2,3,4};
-- M = permutationMatrixInJ(perm, myBasis, J);

-- print "Matrix representation in J:";
-- print toString M;