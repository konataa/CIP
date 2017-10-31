from sympy import *
#from sympy import Symbol, solve
from sympy.solvers.solveset import linsolve

def symsolCIP3():
    # Set symbolic parameters
    x, a, b, c, d, e, h, ui, ui1, gi, gi1, ro = symbols('x, a, b, c, d, e, h, ui, ui1, gi, gi1, ro')
    # Set interpolation polynomials
    f = a + b * x + c * x ** 2 + d * x ** 3
    # f = a + b*x + c*x**2
    df = f.diff(x)
    int_f = f.integrate(x)
    # Set constrains
    eq01 = f.subs(x, 0)
    eq02 = df.subs(x, 0)
    eq03 = f.subs(x, -h)
    eq04 = df.subs(x, -h)
    eq05 = -int_f.subs(x, -h)
    # Solve system under constrains
    sol = linsolve([eq01 - ui, eq02 - gi, eq03 - ui1, eq04 - gi1], (a, b, c, d))
    # sol = linsolve([eq01 - ui, eq03 - ui1, eq05 - ro], (a, b, c))
    # print(sol)
    # Get coefficients
    sol_get = next(iter(sol))
    a_expr = sol_get[0]
    b_expr = sol_get[1]
    c_expr = sol_get[2]
    d_expr = sol_get[3]
    # e_expr = sol_get[4]
    print(sol_get)

    print("d = ")
    print(a_expr)
    print("\n")
    print("c = ")
    print(b_expr)
    print("\n")
    print("b = ")
    print(c_expr)
    print("a = ")
    print(d_expr)
    print("\n")



def symsolCIP5():
    # Set symbolic parameters
    x, a, b, c, d, e, f, h, ui, ui1, gi, gi1, ggi, ggi1 = symbols('x, a, b, c, d, e, f, h, ui, ui1, gi, gi1, ggi, ggi1')
    # Set interpolation polynomial
    func = f + e * x + d * x**2 + c * x**3 + b * x**4 + a * x**5
    df = func.diff(x)
    ddf = df.diff(x)
    eq01 = func.subs(x, 0)
    eq02 = df.subs(x, 0)
    eq03 = func.subs(x, -h)
    eq04 = df.subs(x, -h)
    eq05 = ddf.subs(x, 0)
    eq06 = ddf.subs(x, -h)
    # Solve system under constrains
    sol = linsolve([eq01 - ui, eq02 - gi, eq03 - ui1, eq04 - gi1, eq05 - ggi, eq06 - ggi1], (a, b, c, d, e, f))
    sol_get = next(iter(sol))
    a_expr = sol_get[0]
    b_expr = sol_get[1]
    c_expr = sol_get[2]
    d_expr = sol_get[3]
    e_expr = sol_get[4]
    f_expr = sol_get[5]

    print("a = ")
    print(a_expr)
    print("b = ")
    print(b_expr)
    print("c = ")
    print(c_expr)
    print("d = ")
    print(d_expr)
    print("e = ")
    print(e_expr)
    print("f = ")
    print(f_expr)

def symsolCIP7():
    # Set symbolic parameters
    x, a, b, c, d, e, f, g, l, h, ui, ui1, gi, gi1, ggi, ggi1, gggi, gggi1 = symbols('x, a, b, c, d, e, f, g, l, h, ui, ui1, gi, gi1, ggi, ggi1, gggi, gggi1')
    # Set interpolation polynomial
    func = l + g * x + f * x**2 + e * x**3 + d * x**4 + c * x**5 + b * x**6 + a * x**7
    df = func.diff(x)
    ddf = df.diff(x)
    dddf = ddf.diff(x)
    eq01 = func.subs(x, 0)
    eq02 = df.subs(x, 0)
    eq03 = func.subs(x, -h)
    eq04 = df.subs(x, -h)
    eq05 = ddf.subs(x, 0)
    eq06 = ddf.subs(x, -h)
    eq07 = dddf.subs(x, 0)
    eq08 = dddf.subs(x, -h)
    # Solve system under constrains
    sol = linsolve([eq01 - ui, eq02 - gi, eq03 - ui1, eq04 - gi1, eq05 - ggi, eq06 - ggi1, eq07 - gggi, eq08 - gggi1], (a, b, c, d, e, f, g, l))
    sol_get = next(iter(sol))
    a_expr = sol_get[0]
    b_expr = sol_get[1]
    c_expr = sol_get[2]
    d_expr = sol_get[3]
    e_expr = sol_get[4]
    f_expr = sol_get[5]
    g_expr = sol_get[6]
    l_expr = sol_get[7]

    print("a = ")
    print(a_expr)
    print("b = ")
    print(b_expr)
    print("c = ")
    print(c_expr)
    print("d = ")
    print(d_expr)
    print("e = ")
    print(e_expr)
    print("f = ")
    print(f_expr)
    print("g = ")
    print(g_expr)
    print("l = ")
    print(l_expr)

symsolCIP5()
# symsolCIP7()
