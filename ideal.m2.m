R = QQ[x,y,z]
f = x^3+y^3+z^3
S = R[t]/(t*f-1)

-- t = 1/f

I = ideal (x^2*t,x^2*t,x^2*t)
for i from 0 to 10 list hilbertFunction(i,S/I)