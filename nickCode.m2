k = ZZ/101
R = k[x_0..x_5]
M = genericSkewMatrix(R,4)
I = pfaffians(4,M)
R' = R/I
G = Proj(R')
Q = sheaf coker sub(M, R')
S = sheaf ker sub(M, R')
F = S*Q
for i from 0 to 4 list rank HH^i(F)