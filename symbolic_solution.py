from sympy import *
#from sympy import Symbol, solve
from sympy.solvers.solveset import linsolve

# Set symbolic parameters
x, a, b, c, d, e, h, ui, ui1, gi, gi1, ro = symbols('x, a, b, c, d, e, h, ui, ui1, gi, gi1, ro')
# Set interpolation polynomials
#f = a + b*x + c*x**2 + d*x**3
f = a + b*x + c*x**2
df = f.diff(x)
int_f = f.integrate(x)
# Set constrains
eq01 = f.subs(x, 0)
#eq02 = df.subs(x, 0)
eq03 = f.subs(x, -h)
#eq04 = df.subs(x, -h)
eq05 = -int_f.subs(x, -h)
# Solve system under constrains
#sol = linsolve([eq01 - ui, eq02 - gi, eq03 - ui1, eq04 - gi1, eq05 - ro], (a, b, c, d, e))
sol = linsolve([eq01 - ui, eq03 - ui1, eq05 - ro], (a, b, c))
#print(sol)
# Get coefficients
sol_get = next(iter(sol))
a_expr = sol_get[0]
b_expr = sol_get[1]
c_expr = sol_get[2]
#d_expr = sol_get[3]
#e_expr = sol_get[4]
print(sol_get)

print("e = ")
print(a_expr)
print("\n")
print("d = ")
print(b_expr)
print("\n")
print("c = ")
print(c_expr)
print("\n")