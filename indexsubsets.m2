-- Index all 3-subsets of {0..9}
S3 = subsets(toList(0..9), 3);
nS = #S3;
idxOf = hashTable apply(nS, k -> (S3#k, k));

-- Sign of permutation (given as a list)
signPerm = lst -> (
    inv := sum flatten apply(#lst, i -> apply(i, j -> if lst#j > lst#i then 1 else 0));
    if even inv then 1 else -1
);

-- Sparse action of E_{r,c} on basis element basisIdx
-- Returns list of (coefficient, resultIndex) pairs
actE = (r, c, bi) -> (
    abc := S3#bi;
    flatten apply(3, pos -> (
        if abc#pos == c then (
            nv := replace(pos, r, toList abc);
            if #(unique nv) == 3 then (
                s := sort nv;
                perm := apply(s, x -> position(nv, y -> y == x));
                {(signPerm perm, idxOf#s)}
            ) else {}
        ) else {}
    ))
);

-- Build 120x120 action matrix for E_{r,c}
matE = (r, c) -> (
    M := mutableMatrix(QQ, nS, nS);
    for j from 0 to nS-1 do
        for p in actE(r, c, j) do
            M_(p#1, j) = M_(p#1, j) + p#0;
    matrix M
);

-- Helper: sparse list of (coeff, index) -> column vector in QQ^120
pI = (i, j, k) -> idxOf#{i, j, k};
toVec = coeffs -> (
    v := mutableMatrix(QQ, nS, 1);
    for p in coeffs do v_(p#1, 0) = v_(p#1, 0) + p#0;
    matrix v
);

a1  = toVec {(1,pI(0,1,2)),(1,pI(0,1,3)),(-1,pI(0,2,6)),(1,pI(1,3,6)),(1,pI(2,3,6)),(-1,pI(4,5,8)),(1,pI(5,7,8)),(1,pI(4,7,9)),(1,pI(5,7,9)),(1,pI(4,8,9))};
a2  = toVec {(1,pI(0,1,4)),(1,pI(2,5,6)),(1,pI(3,6,7)),(-1,pI(5,6,7)),(-1,pI(0,2,8)),(1,pI(0,4,8)),(-1,pI(2,5,8)),(1,pI(1,3,9)),(-1,pI(1,4,9)),(1,pI(3,7,9))};
a3  = toVec {(1,pI(0,1,5)),(-1,pI(3,4,6)),(-1,pI(0,2,7)),(1,pI(0,5,7)),(1,pI(1,3,8)),(1,pI(3,4,8)),(-1,pI(1,5,8)),(-1,pI(2,6,9)),(1,pI(4,6,9)),(1,pI(2,7,9))};
a4  = toVec {(1,pI(0,2,3)),(1,pI(1,2,3)),(-1,pI(0,1,6)),(1,pI(1,2,6)),(-1,pI(0,3,6)),(-1,pI(4,5,7)),(1,pI(4,7,8)),(-1,pI(4,5,9)),(1,pI(5,8,9)),(1,pI(7,8,9))};
a5  = toVec {(1,pI(0,2,4)),(1,pI(3,5,6)),(-1,pI(1,3,7)),(-1,pI(3,5,7)),(-1,pI(2,4,8)),(1,pI(2,6,8)),(-1,pI(5,6,8)),(-1,pI(0,1,9)),(1,pI(0,4,9)),(-1,pI(1,7,9))};
a6  = toVec {(1,pI(0,3,4)),(1,pI(0,4,5)),(1,pI(0,5,6)),(-1,pI(2,3,7)),(-1,pI(3,4,7)),(-1,pI(1,2,8)),(-1,pI(2,7,8)),(1,pI(1,6,9)),(-1,pI(5,6,9)),(-1,pI(1,8,9))};
a7  = toVec {(1,pI(0,4,6)),(-1,pI(1,5,6)),(1,pI(4,5,6)),(1,pI(0,3,7)),(-1,pI(0,4,7)),(-1,pI(2,3,8)),(-1,pI(3,7,8)),(-1,pI(1,2,9)),(1,pI(1,5,9)),(-1,pI(2,8,9))};
a8  = toVec {(1,pI(1,2,4)),(1,pI(2,3,5)),(-1,pI(2,4,5)),(1,pI(1,4,7)),(-1,pI(1,6,7)),(1,pI(0,6,8)),(1,pI(6,7,8)),(-1,pI(0,3,9)),(-1,pI(3,5,9)),(1,pI(0,8,9))};
a9  = toVec {(1,pI(1,2,5)),(-1,pI(1,4,5)),(1,pI(1,4,6)),(1,pI(0,6,7)),(-1,pI(4,6,7)),(-1,pI(0,3,8)),(1,pI(0,7,8)),(1,pI(2,3,9)),(-1,pI(2,5,9)),(1,pI(3,8,9))};
a10 = toVec {(1,pI(1,3,4)),(-1,pI(0,2,5)),(1,pI(2,5,7)),(-1,pI(2,6,7)),(1,pI(0,1,8)),(1,pI(1,4,8)),(-1,pI(0,5,8)),(-1,pI(3,4,9)),(1,pI(3,6,9)),(1,pI(6,7,9))};
a11 = toVec {(1,pI(1,3,5)),(1,pI(2,4,6)),(1,pI(0,1,7)),(1,pI(1,5,7)),(-1,pI(3,5,8)),(1,pI(3,6,8)),(-1,pI(4,6,8)),(-1,pI(0,2,9)),(-1,pI(2,4,9)),(1,pI(0,7,9))};
a12 = toVec {(1,pI(2,3,4)),(-1,pI(0,3,5)),(-1,pI(3,4,5)),(1,pI(1,2,7)),(1,pI(2,4,7)),(-1,pI(1,6,8)),(1,pI(1,7,8)),(-1,pI(0,5,9)),(1,pI(0,6,9)),(1,pI(6,8,9))};

W = a1 | a2 | a3 | a4 | a5 | a6 | a7 | a8 | a9 | a10 | a11 | a12;
rW = rank W;
print("Rank of W: " | toString rW);

-- Check gl10-invariance: for all E_{r,c}, image of W must stay in col(W)
isInvariant := true;
for r from 0 to 9 do
    for c from 0 to 9 do (
        ME := matE(r, c);
        print(rank (W | ME*W));
        if rank(W | ME * W) != rW then (
            isInvariant = false;
            print("NOT invariant under E_(" | toString r | "," | toString c | ")");
        );
    );
if isInvariant then print "W is gl10-invariant." else print "W is NOT gl10-invariant.";