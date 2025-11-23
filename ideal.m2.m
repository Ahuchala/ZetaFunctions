R = QQ[x,y,z]
f = x^3+y^3+z^3
S = R[t]/(t*f-1)

R = QQ[x,y,z,t]/(t*(x^3+y^3+z^3)-1,xk = ZZ/101
R = k[x_0..x_5]
M = genericSkewMatrix(R,4)
I = pfaffians(4,M)
R' = R/I
G = Proj(R')
Q = sheaf coker sub(M, R')
F = dual Q ** Q(-2)
for i from 0 to 4 list rank HH^i(F)^2*t-x^4,y^2*t-y^4,z^2*t-z^4)

-- t = 1/f

I = ideal (x^2*t,x^2*t,x^2*t)
for i from 0 to 10 list hilbertFunction(i,S/I)

-- R[t]/(x^2*t,x^2*t,x^2*t,t*f-1)