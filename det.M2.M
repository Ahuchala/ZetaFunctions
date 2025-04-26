k = 2
n = 3

R = QQ[x_{1,1}..x_{n,k}][dx_{1,1}..dx_{n,k}, SkewCommutative => true];

j = 1;
M = transpose matrix {
   flatten(for i from 1 to k list {x_{1,i}..x_{n,i}},
   for i from 1 to n-k list {dx_{1,j}..dx_{n,j}})
};

j = 2;
N = transpose matrix {
   flatten(for i from 1 to k list {x_{1,i}..x_{n,i}},
   for i from 1 to n-k list {dx_{1,j}..dx_{n,j}})
};
// 1/(n-k)!^k normalization
Omega = (det M * det N)/(n-k)!^k


// matrix {
//     flatten(for i from 1 to k list {x_{1,i}..x_{n,i}})
// }
ls =  {x_{1, 1}, x_{1, 2}, x_{2, 1}, x_{2, 2}, x_{3, 1}, x_{3, 2}}
dp_{1,2} = sum (for i from 0 to 5 list diff(ls#i,p_{1,2})*((vars R)_i)_0)
dp_{1,3} = sum (for i from 0 to 5 list diff(ls#i,p_{1,3})*((vars R)_i)_0)
dp_{2,3} = sum (for i from 0 to 5 list diff(ls#i,p_{2,3})*((vars R)_i)_0)

Omega - ( p_{1,2} * dp_{1,3} * dp_{2,3} -  p_{1,3} * dp_{1,2} * dp_{2,3} +  p_{2,3} * dp_{1,2} * dp_{1,3})