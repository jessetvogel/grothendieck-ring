import sympy as sp
from system import *

# Some tests
x, y, z, w = sp.symbols('x y z w')

print('Class of affine line = {}'.format(sp.factor(System({ x }).compute_class())))
print('Class of {{ x*y = 0 }} = {}'.format(sp.factor(System({ x, y }, [ x * y ]).compute_class())))
print('Class of {{ x*y*z = 0 }} = {}'.format(sp.factor(System({ x, y, z }, [ x * y * z ]).compute_class())))
print('Class of {{ x*y*z*w = 0 }} = {}'.format(sp.factor(System({ x, y, z, w }, [ x * y * z * w ]).compute_class())))

print('Class of {{ x^2 + y^2 = 0 }} = {}'.format(sp.factor(System({ x, y }, [ x**2 + y**2 ]).compute_class())))

# Some linear groups
a, b, c, d, e, f, g, h, i = sp.symbols('a b c d e f g h i')

print('Class of SL_2(C) = {}'.format(sp.factor(System({ a, b, c, d }, [ a*d - b*c - 1 ]).compute_class())))
print('Class of GL_2(C) = {}'.format(sp.factor(System({ a, b, c, d }, [], [ a*d - b*c ]).compute_class())))

print('Class of SL_3(C) = {}'.format(sp.factor(System({ a, b, c, d, e, f, g, h, i }, [ a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g) - 1 ]).compute_class())))
print('Class of GL_3(C) = {}'.format(sp.factor(System({ a, b, c, d, e, f, g, h, i }, [], [ a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g) ]).compute_class())))
