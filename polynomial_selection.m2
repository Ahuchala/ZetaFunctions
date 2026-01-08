k = ZZ/101
G = Grassmannian(2,9,CoefficientRing => k);
R = ring G;

-- A = {0, 4, 5, 6, 1, 2, 3, 7, 8, 9} -- secretly (0 1)
A = {6, 0, 5, 4, 3, 9, 8, 2, 1, 7} -- secretly (0 1 4 2)
-- B = {6, 4, 0, 5, 8, 3, 9, 1, 7, 2} -- secretly (0 1 4 3)

-- A = { 0, 4, 3, 6, 1, 2, 3, 7, 8, 9}; -- secretly (0 1)
-- B = { 4, 5, 6, 0, 7, 8, 1, 9, 2, 3};

-- A = { 0, 5, 4, 6, 2, 1, 3, 7, 9, 8 } -- secretly (0 1) (2 3)
B = { 4, 8, 7, 1, 6, 5, 0, 9, 3, 2 } -- secretly (0 2 4)

fooA = I -> (
   J := for i in I list A_i;
   sgn := (if J_0 < J_1 then 1 else -1) *
          (if J_0 < J_2 then 1 else -1) *
          (if J_1 < J_2 then 1 else -1);
   return sgn*p_(toSequence sort J)
)

barA = map(R, R, subsets(10, 3) / fooA)

fooB = I -> (
   J := for i in I list B_i;
   sgn := (if J_0 < J_1 then 1 else -1) *
          (if J_0 < J_2 then 1 else -1) *
          (if J_1 < J_2 then 1 else -1);
   return sgn*p_(toSequence sort J)
)

barB = map(R, R, subsets(10, 3) / fooB)
L = {p_(0,1,5)}
-- L = { p_(0,1,2) } -- could also have tried p_(0,4,7) or p_(0,4,5) or p_(0,1,9)
L = unique sort(L | (L / barA) | (L / barB)) -- repeat until it stabilizes
