

# This file was *autogenerated* from the file debug.sage
from sage.all_cmdline import *   # import sage library

_sage_const_0 = Integer(0); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_34 = Integer(34); _sage_const_6 = Integer(6); _sage_const_3 = Integer(3); _sage_const_23 = Integer(23); _sage_const_15 = Integer(15); _sage_const_5 = Integer(5); _sage_const_4 = Integer(4); _sage_const_20 = Integer(20); _sage_const_61 = Integer(61); _sage_const_29 = Integer(29); _sage_const_126 = Integer(126); _sage_const_136 = Integer(136); _sage_const_571 = Integer(571); _sage_const_201 = Integer(201); _sage_const_90 = Integer(90); _sage_const_7 = Integer(7)
load("mat_mul.sage")
A = matrix([(_sage_const_0 , _sage_const_2 , _sage_const_1 , _sage_const_1 , -_sage_const_1 , _sage_const_0 , _sage_const_34 ), (_sage_const_1 , _sage_const_0 , _sage_const_6 , _sage_const_1 , -_sage_const_1 , -_sage_const_2 , -_sage_const_1 ), (_sage_const_1 , _sage_const_1 , _sage_const_0 , -_sage_const_1 , -_sage_const_3 , _sage_const_23 , _sage_const_15 ), (_sage_const_5 , _sage_const_4 , _sage_const_0 , _sage_const_1 , -_sage_const_6 , _sage_const_4 , _sage_const_2 ), (_sage_const_1 , -_sage_const_20 , -_sage_const_1 , -_sage_const_1 , _sage_const_1 , _sage_const_0 , _sage_const_3 ), (_sage_const_61 , _sage_const_1 , _sage_const_3 , _sage_const_2 , -_sage_const_1 , -_sage_const_1 , -_sage_const_5 ), (-_sage_const_2 , -_sage_const_3 , -_sage_const_1 , -_sage_const_2 , _sage_const_0 , _sage_const_0 , _sage_const_2 )])
v = vector([_sage_const_3 , _sage_const_29 , -_sage_const_1 , _sage_const_1 , _sage_const_0 , _sage_const_0 , _sage_const_2 ])
assert multiply_matrix_vector(A,v) == vector((_sage_const_126 , -_sage_const_4 , _sage_const_61 , _sage_const_136 , -_sage_const_571 , _sage_const_201 , -_sage_const_90 ))

A = random_matrix(ZZ,_sage_const_7 ,_sage_const_7 )
v = vector(random_matrix(ZZ,_sage_const_7 ,_sage_const_1 ))
assert multiply_matrix_vector(A,v) == A*v
