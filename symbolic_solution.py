from sympy import *
#from sympy import Symbol, solve
from sympy.solvers.solveset import linsolve
from sympy.polys.polyfuncs import horner

def symbolCIP3():
    # Set symbolic parameters
    x, a, b, c, d, e, ui, ui1, gi, gi1, ro = symbols('x, a, b, c, d, e, ui, ui1, gi, gi1, ro')
    h = symbols('h')

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

    # Solve system under constrains
    sol = linsolve([eq01 - ui, eq02 - gi, eq03 - ui1, eq04 - gi1], (a, b, c, d))

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



def symbolCIP5():
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
    sol = linsolve([eq01 - ui, eq02 - gi, eq05 - ggi, eq03 - ui1, eq04 - gi1, eq06 - ggi1], (a, b, c, d, e, f))
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

def symbolCIP7():
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
    sol = linsolve([eq01 - ui, eq02 - gi, eq03 - ui1, eq05 - ggi, eq04 - gi1, eq07 - gggi, eq06 - ggi1, eq08 - gggi1], (a, b, c, d, e, f, g, l))
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


def symbolCIP9():
    # Set symbolic parameters
    x, a, b, c, d, e, f, g, l, k, q, h, ui, ui1, gi, gi1, ggi, ggi1, gggi, gggi1, ggggi, ggggi1 \
        = symbols('x, a, b, c, d, e, f, g, l, k, q, h, ui, ui1, gi, gi1, ggi, ggi1, gggi, gggi1, ggggi, ggggi1')

    # Set interpolation polynomial
    func = q + k * x + l * x**2 + g * x**3 + f * x**4 + e * x**5 + d * x**6 + c * x**7 + b * x**8 + a * x**9
    df = func.diff(x)
    ddf = df.diff(x)
    dddf = ddf.diff(x)
    ddddf = dddf.diff(x)
    eq01 = func.subs(x, 0)
    eq02 = df.subs(x, 0)
    eq03 = func.subs(x, -h)
    eq04 = df.subs(x, -h)
    eq05 = ddf.subs(x, 0)
    eq06 = ddf.subs(x, -h)
    eq07 = dddf.subs(x, 0)
    eq08 = dddf.subs(x, -h)
    eq09 = ddddf.subs(x, 0)
    eq10 = ddddf.subs(x, -h)

    # Solve system under constrains

    sol = linsolve([eq01 - ui, eq02 - gi, eq03 - ui1, eq05 - ggi, eq04 - gi1, eq07 - gggi, eq06 - ggi1, eq08 - gggi1, eq09 - ggggi, eq10 - ggggi1], (a, b, c, d, e, f, g, l, k, q))
    sol_get = next(iter(sol))
    a_expr = sol_get[0]
    b_expr = sol_get[1]
    c_expr = sol_get[2]
    d_expr = sol_get[3]
    e_expr = sol_get[4]
    f_expr = sol_get[5]
    g_expr = sol_get[6]
    l_expr = sol_get[7]
    k_expr = sol_get[8]
    q_expr = sol_get[9]
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
    print("k = ")
    print(k_expr)
    print("q = ")
    print(q_expr)

# symsolCIP5()
# symsolCIP7()
# symsolCIP9()


def symboldeltaP3():
    x, a, h, t, nu, c1, c2, c3, c4, ui, ui1, ri, ri1, uxx, u3x, rx, rxx, _uni \
        = symbols('x, a, h, t, nu, c1, c2, c3, c4, ui, ui1, ri, ri1, uxx, u3x, rx, rxx, _uni')

    u = c1 + c2 * x + c3 * x ** 2 + c4 * x ** 3

    du = u.diff(x)
    print(du)
    # Set constrains
    eq01 = u.subs(x, 0)
    eq02 = du.subs(x, 0)
    eq03 = u.subs(x, -h)
    eq04 = du.subs(x, -h)
    # Solve system under constrains
    sol = linsolve([eq01 - ui, eq02 - ri, eq03 - ui1, eq04 - ri1], (c1, c2, c3, c4))
    # Get coefficients
    sol_get = next(iter(sol))
    c1_expr = sol_get[0].subs(ri, ri / h)
    c2_expr = sol_get[1].subs(ri, ri / h)
    c3_expr = sol_get[2].subs(ri, ri / h).subs(ri1, ri1 / h)
    c4_expr = sol_get[3].subs(ri, ri / h).subs(ri1, ri1 / h)
    print(("c1 = {0}".format(c1_expr)))
    print(("c2 = {0}".format(c2_expr)))
    print(("c3 = {0}".format(c3_expr)))
    print(("c4 = {0}".format(c4_expr)))

    # print("u = {0}".format(u))
    # nu = a * t / h
    # Set constrains
    ui = u.subs(x, 0)
    ui1 = u.subs(x, -h)

    ux = u.diff(x)
    uxi = ux.subs(x, 0)

    uxx = ux.diff(x)
    uxxi = uxx.subs(x, 0)

    u3x = uxx.diff(x)
    u3xi = u3x.subs(x, 0)

    r = h * u.diff(x)
    ri = r.subs(x, 0)

    rx = r.diff(x)
    rxi = rx.subs(x, 0)

    rxx = rx.diff(x)
    rxxi = rxx.subs(x, 0)

    # uxx = ux.diff(x)
    # u3x = uxx.diff(x)

    # print("ui = {0}, ui1 = {1}, ri = {2}, ri1 = {3}".format(simplify(ui), simplify(ui1), ri, ri1))
    #
    # c1 = ui
    # c2 = ri / h
    #
    # c3 = (3 * (ui1 - ui) + 2 * ri + ri1) / h**2
    # c4 = (2 * (ui1 - ui) + ri + ri1) / h**3

    # uxx = 2 * c3
    # u3x = 6 * c4
    #
    # rx = h * 2 * c3
    # rxx = h * 6 * c4

    uni = ui - a * t * uxi + a**2 * t**2 / 2 * uxxi - a**3 * t**3 / 6 * u3xi
    # uni = simplify(uni)
    _uni = uni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr).subs(c4, c4_expr).subs(a * t / h, nu)
    _uni = collect(expand(_uni), nu)
    print(_uni)

    rni = simplify(ri - a * t * rxi + a**2 * t**2 / 2 * rxxi)
    _rni = rni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr).subs(c4, c4_expr).subs(a * t / h, nu)
    print(collect(expand(_rni), nu))

def symboldeltaP5():
    x, a, h, tau, nu, ui, ui1, ri, ri1, u3x, u4x, u5x, rxi, rxi1, rxx, r3x, r4x \
        = symbols('x, a, h, tau, nu, ui, ui1, ri, ri1, u3x, u4x, u5x, rxi, rxi1, rxx, r3x, r4x')

    c1, c2, c3, c4, c5, c6 = symbols('c1, c2, c3, c4, c5, c6')

    _uxi, _uxxi, _ux1i, _uxx1i = symbols('_uxi, _uxxi, _ux1i, _uxx1i')

    u = c1 + c2 * x + c3 * x ** 2 + c4 * x ** 3 + c5 * x ** 4 + c6 * x ** 5

    du = u.diff(x)
    ddu = du.diff(x)
    eq01 = u.subs(x, 0)
    eq02 = du.subs(x, 0)
    eq03 = u.subs(x, -h)
    eq04 = du.subs(x, -h)
    eq05 = ddu.subs(x, 0)
    eq06 = ddu.subs(x, -h)
    # Solve system under constrains
    sol = linsolve([eq01 - ui, eq02 - _uxi, eq03 - ui1, eq04 - _ux1i, eq05 - _uxxi, eq06 - _uxx1i],
                   (c1, c2, c3, c4, c5, c6))

    sol_get = next(iter(sol))
    # c1_expr = sol_get[0].subs(ri, ri / h).subs(ri1, ri1 / h).subs(rxi, rxi / h**2).subs(rxi1, rxi1 / h**2)
    # c2_expr = sol_get[1].subs(ri, ri / h).subs(ri1, ri1 / h).subs(rxi, rxi / h**2).subs(rxi1, rxi1 / h**2)
    # c3_expr = sol_get[2].subs(ri, ri / h).subs(ri1, ri1 / h).subs(rxi, rxi / h**2).subs(rxi1, rxi1 / h**2)
    # c4_expr = sol_get[3].subs(ri, ri / h).subs(ri1, ri1 / h).subs(rxi, rxi / h**2).subs(rxi1, rxi1 / h**2)
    # c5_expr = sol_get[4].subs(ri, ri / h).subs(ri1, ri1 / h).subs(rxi, rxi / h**2).subs(rxi1, rxi1 / h**2)
    # c6_expr = sol_get[5].subs(ri, ri / h).subs(ri1, ri1 / h).subs(rxi, rxi / h**2).subs(rxi1, rxi1 / h**2)
    c1_expr = sol_get[0]
    c2_expr = sol_get[1].subs(_uxi, ri / h)
    c3_expr = sol_get[2].subs(_uxxi, rxi / h ** 2)
    c4_expr = sol_get[3].subs(h * _uxi, ri).subs(h * _ux1i, ri1).subs(h ** 2 * _uxxi, rxi).subs(h ** 2 * _uxx1i, rxi1)
    c5_expr = sol_get[4].subs(h * _uxi, ri).subs(h * _ux1i, ri1).subs(h ** 2 * _uxxi, rxi).subs(h ** 2 * _uxx1i, rxi1)
    c6_expr = sol_get[5].subs(h * _uxi, ri).subs(h * _ux1i, ri1).subs(h ** 2 * _uxxi, rxi).subs(h ** 2 * _uxx1i, rxi1)
    print("c1 = {0}".format(c1_expr))
    print("c2 = {0}".format(c2_expr))
    print("c3 = {0}".format(c3_expr))
    print("c4 = {0}".format(c4_expr))
    print("c5 = {0}".format(c5_expr))
    print("c6 = {0}".format(c6_expr))

    # print("u = {0}".format(u))
    # Set constrains
    ui = u.subs(x, 0)

    ux = u.diff(x)
    uxi = ux.subs(x, 0)

    uxx = ux.diff(x)
    uxxi = uxx.subs(x, 0)

    u3x = uxx.diff(x)
    u3xi = u3x.subs(x, 0)

    u4x = u3x.diff(x)
    u4xi = u4x.subs(x, 0)

    u5x = u4x.diff(x)
    u5xi = u5x.subs(x, 0)

    r = h * ux
    ri = r.subs(x, 0)

    rx = r.diff(x)
    rxi = rx.subs(x, 0)

    rxx = rx.diff(x)
    rxxi = rxx.subs(x, 0)

    r3x = rxx.diff(x)
    r3xi = r3x.subs(x, 0)

    r4x = r3x.diff(x)
    r4xi = r4x.subs(x, 0)

    # uxx = ux.diff(x)
    # u3x = uxx.diff(x)

    # print("ui = {0}, ui1 = {1}, ri = {2}, ri1 = {3}".format(simplify(ui), simplify(ui1), ri, ri1))
    #
    # c1 = ui
    # c2 = ri / h
    #
    # c3 = (3 * (ui1 - ui) + 2 * ri + ri1) / h**2
    # c4 = (2 * (ui1 - ui) + ri + ri1) / h**3

    # uxx = 2 * c3
    # u3x = 6 * c4
    #
    # rx = h * 2 * c3
    # rxx = h * 6 * c4

    uni = (ui - a * tau * uxi + a**2 * tau**2 / 2 * uxxi - a**3 * tau**3 / 6
           * u3xi + a**4 * tau**4 / 24 * u4xi - a ** 5 * tau ** 5 / 120 * u5xi)

    _uni = uni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr)\
        .subs(c4, c4_expr).subs(c5, c5_expr).subs(c6, c6_expr).subs(a * tau / h, nu)

    rni = ri - a * tau * rxi + a**2 * tau**2 / 2 * rxxi - a**3 * tau**3 / 6 * r3xi + a**4 * tau**4 / 24 * r4xi

    _rni = rni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr)\
        .subs(c4, c4_expr).subs(c5, c5_expr).subs(c6, c6_expr).subs(a * tau / h, nu)

    rxni = rxi - a * tau * rxxi + a**2 * tau**2 / 2 * r3xi - a**3 * tau**3 / 6 * r4xi

    _rxni = rxni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr)\
        .subs(c4, c4_expr).subs(c5, c5_expr).subs(c6, c6_expr).subs(a * tau / h, nu) * h

    print(simplify(_uni))
    print(simplify(_rni))
    print(simplify(_rxni))

def symboldeltaP7():
    x, a, h, tau, nu, c1, c2, c3, c4, c5, c6, c7, c8 = symbols('x, a, h, tau, nu, c1, c2, c3, c4, c5, c6, c7, c8')

    ui, ui1, uxi, uxi1, uxxi, uxxi1, u3xi, u3xi1 = symbols('ui, ui1, ri, ri1, uxxi, uxx1i, u3xi, u3x1i')

    ri, ri1, rxi, rxi1, rxxi, rxxi1 = symbols('ri, ri1, rxi, rxi1, rxxi, rxxi1')

    u = c1 + c2 * x + c3 * x ** 2 + c4 * x ** 3 + c5 * x ** 4 + c6 * x ** 5 + c7 * x ** 6 + c8 * x ** 7

    du = u.diff(x)
    ddu = du.diff(x)
    dddu = ddu.diff(x)
    eq01 = u.subs(x, 0)
    eq02 = du.subs(x, 0)
    eq03 = u.subs(x, -h)
    eq04 = du.subs(x, -h)
    eq05 = ddu.subs(x, 0)
    eq06 = ddu.subs(x, -h)
    eq07 = dddu.subs(x, 0)
    eq08 = dddu.subs(x, -h)
    # Solve system under constrains
    sol = linsolve([eq01 - ui, eq02 - uxi, eq03 - ui1, eq04 - uxi1, eq05 - uxxi, eq06 - uxxi1, eq07 - u3xi, eq08 - u3xi1],
                   (c1, c2, c3, c4, c5, c6, c7, c8))
    sol_get = next(iter(sol))
    c1_expr = sol_get[0]
    c2_expr = sol_get[1].subs(uxi, ri / h)
    c3_expr = sol_get[2].subs(uxxi, rxi / h ** 2)
    c4_expr = sol_get[3].subs(u3xi, rxxi / h ** 3)
    c5_expr = sol_get[4].subs(h * uxi, ri).subs(h * uxi1, ri1).subs(h ** 2 * uxxi, rxi)\
        .subs(h ** 2 * uxxi1, rxi1).subs(h ** 3 * u3xi, rxxi).subs(h ** 3 * u3xi1, rxxi1)
    c6_expr = sol_get[5].subs(h * uxi, ri).subs(h * uxi1, ri1).subs(h ** 2 * uxxi, rxi)\
        .subs(h ** 2 * uxxi1, rxi1).subs(h ** 3 * u3xi, rxxi).subs(h ** 3 * u3xi1, rxxi1)
    c7_expr = sol_get[6].subs(h * uxi, ri).subs(h * uxi1, ri1).subs(h ** 2 * uxxi, rxi)\
        .subs(h ** 2 * uxxi1, rxi1).subs(h ** 3 * u3xi, rxxi).subs(h ** 3 * u3xi1, rxxi1)
    c8_expr = sol_get[7].subs(h * uxi, ri).subs(h * uxi1, ri1).subs(h ** 2 * uxxi, rxi)\
        .subs(h ** 2 * uxxi1, rxi1).subs(h ** 3 * u3xi, rxxi).subs(h ** 3 * u3xi1, rxxi1)

    print("c1 = {0}".format(c1_expr))
    print("c2 = {0}".format(c2_expr))
    print("c3 = {0}".format(c3_expr))
    print("c4 = {0}".format(c4_expr))
    print("c5 = {0}".format(c5_expr))
    print("c6 = {0}".format(c6_expr))
    print("c7 = {0}".format(c7_expr))
    print("c8 = {0}".format(c8_expr))


    # Set constrains
    ui = u.subs(x, 0)

    ux = u.diff(x)
    uxi = ux.subs(x, 0)

    uxx = ux.diff(x)
    uxxi = uxx.subs(x, 0)

    u3x = uxx.diff(x)
    u3xi = u3x.subs(x, 0)

    u4x = u3x.diff(x)
    u4xi = u4x.subs(x, 0)

    u5x = u4x.diff(x)
    u5xi = u5x.subs(x, 0)

    u6x = u5x.diff(x)
    u6xi = u6x.subs(x, 0)

    u7x = u6x.diff(x)
    u7xi = u7x.subs(x, 0)

    r = h * u.diff(x)
    ri = r.subs(x, 0)

    rx = r.diff(x)
    rxi = rx.subs(x, 0)

    rxx = rx.diff(x)
    rxxi = rxx.subs(x, 0)

    r3x = rxx.diff(x)
    r3xi = r3x.subs(x, 0)

    r4x = r3x.diff(x)
    r4xi = r4x.subs(x, 0)

    r5x = r4x.diff(x)
    r5xi = r5x.subs(x, 0)

    r6x = r5x.diff(x)
    r6xi = r6x.subs(x, 0)

    # uxx = ux.diff(x)
    # u3x = uxx.diff(x)

    # print("ui = {0}, ui1 = {1}, ri = {2}, ri1 = {3}".format(simplify(ui), simplify(ui1), ri, ri1))
    #
    # c1 = ui
    # c2 = ri / h
    #
    # c3 = (3 * (ui1 - ui) + 2 * ri + ri1) / h**2
    # c4 = (2 * (ui1 - ui) + ri + ri1) / h**3

    # uxx = 2 * c3
    # u3x = 6 * c4
    #
    # rx = h * 2 * c3
    # rxx = h * 6 * c4

    uni = ui - a * tau * uxi + a**2 * tau**2 / 2 * uxxi - a**3 * tau**3 / 6 * u3xi + a**4 * tau**4 / 24 * u4xi \
          - a ** 5 * tau ** 5 / 120 * u5xi + a ** 6 * tau ** 6 / 720 * u6xi - a ** 7 * tau ** 7 / 5040 * u7xi

    rni = ri - a * tau * rxi + a**2 * tau**2 / 2 * rxxi - a**3 * tau**3 / 6 * r3xi + a**4 * tau**4 / 24 * r4xi \
          - a ** 5 * tau ** 5 / 120 * r5xi + a ** 6 * tau ** 6 / 720 * r6xi

    rxni = rxi - a * tau * rxxi + a**2 * tau**2 / 2 * r3xi - a**3 * tau**3 / 6 * r4xi \
           + a**4 * tau**4 / 24 * r5xi - a ** 5 * tau ** 5 / 120 * r6xi

    rxxni = rxxi - a * tau * r3xi + a**2 * tau**2 / 2 * r4xi - a**3 * tau**3 / 6 * r5xi \
            + a**4 * tau**4 / 24 * r6xi

    uni = uni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr)\
        .subs(c4, c4_expr).subs(c5, c5_expr).subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(a * tau / h, nu)

    rni = rni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr) \
        .subs(c4, c4_expr).subs(c5, c5_expr).subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(a * tau / h, nu)

    rxni = rxni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr) \
        .subs(c4, c4_expr).subs(c5, c5_expr).subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr)\
        .subs(a * tau / h, nu) * h

    rxxni = rxxni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr) \
        .subs(c4, c4_expr).subs(c5, c5_expr).subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr)\
        .subs(a * tau / h, nu) * h ** 2

    print("uni = {0}".format(simplify(uni)))
    print("rni = {0}".format(simplify(rni)))
    print("rxni = {0}".format(simplify(rxni)))
    print("rxxni = {0}".format(simplify(rxxni)))

    # print("rx = {0}, с2 = {1}, с3 = {2}, с4 = {3}, с5 = {4}".format(simplify(rx), simplify(c2), c3, c4, c5))

def symboldeltaP9():
    x, a, h, tau, nu, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10 \
        = symbols('x, a, h, tau, nu, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10')

    ui, ui1, uxi, uxi1, uxxi, uxxi1, u3xi, u3xi1, u4xi, u4xi1 \
        = symbols('ui, ui1, ri, ri1, uxxi, uxx1i, u3xi, u3x1i, u4xi, u4xi1')

    ri, ri1, rxi, rxi1, rxxi, rxxi1, r3xi, r3xi1 = symbols('ri, ri1, rxi, rxi1, rxxi, rxxi1, r3xi, r3xi1')

    u = c1 + c2 * x + c3 * x ** 2 + c4 * x ** 3 + c5 * x ** 4 + c6 * x ** 5 + c7 * x ** 6 + c8 * x ** 7 \
        + c9 * x ** 8 + c10 * x ** 9

    du = u.diff(x)
    ddu = du.diff(x)
    d3u = ddu.diff(x)
    d4u = d3u.diff(x)

    eq01 = u.subs(x, 0)
    eq02 = du.subs(x, 0)
    eq03 = u.subs(x, -h)
    eq04 = du.subs(x, -h)
    eq05 = ddu.subs(x, 0)
    eq06 = ddu.subs(x, -h)
    eq07 = d3u.subs(x, 0)
    eq08 = d3u.subs(x, -h)
    eq09 = d4u.subs(x, 0)
    eq10 = d4u.subs(x, -h)

    # Solve system under constrains

    sol = linsolve([eq01 - ui, eq02 - uxi, eq03 - ui1, eq04 - uxi1, eq05 - uxxi, eq06 - uxxi1, eq07 - u3xi,
                    eq08 - u3xi1, eq09 - u4xi, eq10 - u4xi1], (c1, c2, c3, c4, c5, c6, c7, c8, c9, c10))

    sol_get = next(iter(sol))

    c1_expr = sol_get[0]
    c2_expr = sol_get[1].subs(uxi, ri / h)
    c3_expr = sol_get[2].subs(uxxi, rxi / h ** 2)
    c4_expr = sol_get[3].subs(u3xi, rxxi / h ** 3)
    c5_expr = sol_get[4].subs(u4xi, r3xi / h ** 4)
    c6_expr = sol_get[5].subs(h * uxi, ri).subs(h * uxi1, ri1).subs(h ** 2 * uxxi, rxi) \
        .subs(h ** 2 * uxxi1, rxi1).subs(h ** 3 * u3xi, rxxi).subs(h ** 3 * u3xi1, rxxi1).subs(h ** 4 * u4xi, r3xi)\
        .subs(h ** 4 * u4xi1, r3xi1)
    c7_expr = sol_get[6].subs(h * uxi, ri).subs(h * uxi1, ri1).subs(h ** 2 * uxxi, rxi) \
        .subs(h ** 2 * uxxi1, rxi1).subs(h ** 3 * u3xi, rxxi).subs(h ** 3 * u3xi1, rxxi1).subs(h ** 4 * u4xi, r3xi)\
        .subs(h ** 4 * u4xi1, r3xi1)
    c8_expr = sol_get[7].subs(h * uxi, ri).subs(h * uxi1, ri1).subs(h ** 2 * uxxi, rxi) \
        .subs(h ** 2 * uxxi1, rxi1).subs(h ** 3 * u3xi, rxxi).subs(h ** 3 * u3xi1, rxxi1).subs(h ** 4 * u4xi, r3xi)\
        .subs(h ** 4 * u4xi1, r3xi1)
    c9_expr = sol_get[8].subs(h * uxi, ri).subs(h * uxi1, ri1).subs(h ** 2 * uxxi, rxi) \
        .subs(h ** 2 * uxxi1, rxi1).subs(h ** 3 * u3xi, rxxi).subs(h ** 3 * u3xi1, rxxi1).subs(h ** 4 * u4xi, r3xi)\
        .subs(h ** 4 * u4xi1, r3xi1)
    c10_expr = sol_get[9].subs(h * uxi, ri).subs(h * uxi1, ri1).subs(h ** 2 * uxxi, rxi) \
        .subs(h ** 2 * uxxi1, rxi1).subs(h ** 3 * u3xi, rxxi).subs(h ** 3 * u3xi1, rxxi1).subs(h ** 4 * u4xi, r3xi)\
        .subs(h ** 4 * u4xi1, r3xi1)
    print("c1 = {0}".format(c1_expr))
    print("c2 = {0}".format(c2_expr))
    print("c3 = {0}".format(c3_expr))
    print("c4 = {0}".format(c4_expr))
    print("c5 = {0}".format(c5_expr))
    print("c6 = {0}".format(c6_expr))
    print("c7 = {0}".format(c7_expr))
    print("c8 = {0}".format(c8_expr))
    print("c9 = {0}".format(c9_expr))
    print("c10 = {0}\n".format(c10_expr))

    # Set constrains
    ui = u.subs(x, 0)

    ux = u.diff(x)
    uxi = ux.subs(x, 0)

    uxx = ux.diff(x)
    uxxi = uxx.subs(x, 0)

    u3x = uxx.diff(x)
    u3xi = u3x.subs(x, 0)

    u4x = u3x.diff(x)
    u4xi = u4x.subs(x, 0)

    u5x = u4x.diff(x)
    u5xi = u5x.subs(x, 0)

    u6x = u5x.diff(x)
    u6xi = u6x.subs(x, 0)

    u7x = u6x.diff(x)
    u7xi = u7x.subs(x, 0)

    u8x = u7x.diff(x)
    u8xi = u8x.subs(x, 0)

    u9x = u8x.diff(x)
    u9xi = u9x.subs(x, 0)

    r = h * u.diff(x)
    ri = r.subs(x, 0)

    rx = r.diff(x)
    rxi = rx.subs(x, 0)

    rxx = rx.diff(x)
    rxxi = rxx.subs(x, 0)

    r3x = rxx.diff(x)
    r3xi = r3x.subs(x, 0)

    r4x = r3x.diff(x)
    r4xi = r4x.subs(x, 0)

    r5x = r4x.diff(x)
    r5xi = r5x.subs(x, 0)

    r6x = r5x.diff(x)
    r6xi = r6x.subs(x, 0)

    r7x = r6x.diff(x)
    r7xi = r7x.subs(x, 0)

    r8x = r7x.diff(x)
    r8xi = r8x.subs(x, 0)

    uni = ui - a * tau * uxi + a ** 2 * tau ** 2 / 2 * uxxi - a ** 3 * tau ** 3 / 6 * u3xi + a ** 4 * tau ** 4 / 24 * u4xi \
          - a ** 5 * tau ** 5 / 120 * u5xi + a ** 6 * tau ** 6 / 720 * u6xi - a ** 7 * tau ** 7 / 5040 * u7xi \
          + a ** 8 * tau ** 8 / 40320 * u8xi - a ** 9 * tau ** 9 / 362880 * u9xi

    rni = ri - a * tau * rxi + a ** 2 * tau ** 2 / 2 * rxxi - a ** 3 * tau ** 3 / 6 * r3xi + a ** 4 * tau ** 4 / 24 * r4xi \
          - a ** 5 * tau ** 5 / 120 * r5xi + a ** 6 * tau ** 6 / 720 * r6xi - a ** 7 * tau ** 7 / 5040 * r7xi \
          + a ** 8 * tau ** 8 / 40320 * r8xi

    rxni = rxi - a * tau * rxxi + a ** 2 * tau ** 2 / 2 * r3xi - a ** 3 * tau ** 3 / 6 * r4xi \
           + a ** 4 * tau ** 4 / 24 * r5xi - a ** 5 * tau ** 5 / 120 * r6xi + a ** 6 * tau ** 6 / 720 * r7xi \
           - a ** 7 * tau ** 7 / 5040 * r8xi

    rxxni = rxxi - a * tau * r3xi + a ** 2 * tau ** 2 / 2 * r4xi - a ** 3 * tau ** 3 / 6 * r5xi \
            + a ** 4 * tau ** 4 / 24 * r6xi - a ** 5 * tau ** 5 / 120 * r7xi + a ** 6 * tau ** 6 / 720 * r8xi

    r3xni = r3xi - a * tau * r4xi + a ** 2 * tau ** 2 / 2 * r5xi - a ** 3 * tau ** 3 / 6 * r6xi \
            + a ** 4 * tau ** 4 / 24 * r7xi - a ** 5 * tau ** 5 / 120 * r8xi

    uni = uni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr) \
        .subs(c4, c4_expr).subs(c5, c5_expr).subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(c9, c9_expr)\
        .subs(c10, c10_expr).subs(a * tau / h, nu)

    rni = rni.subs(c2, c2_expr).subs(c3, c3_expr) \
        .subs(c4, c4_expr).subs(c5, c5_expr).subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(c9, c9_expr)\
        .subs(c10, c10_expr).subs(a * tau / h, nu)

    rxni = rxni.subs(c2, c2_expr).subs(c3, c3_expr) \
        .subs(c4, c4_expr).subs(c5, c5_expr).subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(c9, c9_expr)\
        .subs(c10, c10_expr).subs(a * tau / h, nu) * h

    rxxni = rxxni.subs(c2, c2_expr).subs(c3, c3_expr) \
        .subs(c4, c4_expr).subs(c5, c5_expr).subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(c9, c9_expr)\
        .subs(c10, c10_expr).subs(a * tau / h, nu) * h ** 2

    r3xni = r3xni.subs(c2, c2_expr).subs(c3, c3_expr) \
        .subs(c4, c4_expr).subs(c5, c5_expr).subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(c9, c9_expr)\
        .subs(c10, c10_expr).subs(a * tau / h, nu) * h ** 3

    print("uni = {0}".format(simplify(uni)))
    print("rni = {0}".format(simplify(rni)))
    print("rxni = {0}".format(simplify(rxni)))
    print("rxxni = {0}".format(simplify(rxxni)))
    print("r3xni = {0}".format(simplify(r3xni)))


def get_derivative():
    t, x, y = symbols('t, x, y')

    p = cos(2 * sqrt(2) * pi * t) * sin(2 * pi * x) * sin(2 * pi * y)

    u = -1 / sqrt(2) * sin(2 * sqrt(2) * pi * t) * cos(2 * pi * x) * sin(2 * pi * y)

    v = -1 / sqrt(2) * sin(2 * sqrt(2) * pi * t) * sin(2 * pi * x) * cos(2 * pi * y)

    p_exp = exp(-(10 * (x ** 2 + (y + 0.7) ** 2)))

    p_sin = 1 / pi * sin(pi * x) ** 2 * sin(pi * y) ** 2

    rp = diff(p, x)
    ru = diff(u, x)
    rv = diff(v, x)

    sp = diff(p, y)
    su = diff(u, y)
    sv = diff(v, y)

    rp_exp = diff(p_exp, x)
    sp_exp = diff(p_exp, y)

    rp_sin = diff(p_sin, x)
    sp_sin = diff(p_sin, y)

    print("rp = {0}\n".format(rp))
    print("ru = {0}\n".format(ru))
    print("rv = {0}\n".format(rv))

    print("sp = {0}\n".format(sp))
    print("su = {0}\n".format(su))
    print("sv = {0}\n".format(sv))

    print("rp_exp = {0}\n".format(rp_exp))
    print("sp_exp = {0}\n".format(sp_exp))

    print("rp_sin = {0}\n".format(rp_sin))
    print("sp_sin = {0}\n".format(sp_sin))

def two_dim_deltaP3():
    x, y, a, h, tau, nu, ui, ui1, ri, ri1, uxx, u3x, rx, rxx \
        = symbols('x, y, a, h, tau, nu, ui, ui1, ri, ri1, uxx, u3x, rx, rxx')

    c1, c2, c3, c4, c5, c6, c7, c8, c9, c10 = \
        symbols('c1, c2, c3, c4, c5, c6, c7, c8, c9, c10')

    phi_i_j1, phi_i1_j1, phi_i1_j, phi_i_j, r_i_j1, r_i_j, r_i1_j, s_i_j1, s_i_j, s_i1_j = \
        symbols('phi_i_j1, phi_i1_j1, phi_i1_j, phi_i_j, r_i_j1, r_i_j, r_i1_j, s_i_j1, s_i_j, s_i1_j')

    # p, u, v, rp, ru, rv, sp, su, sv = symbols('p, u, v, rp, ru, rv, sp, su, sv')

    p = c1 + c2 * x + c3 * y + c4 * x * y + c5 * x ** 2 + c6 * y ** 2 + c7 * x ** 2 * y + c8 * x * y ** 2 \
        + c9 * x ** 3 + c10 * y ** 3

    u = c1 + c2 * x + c3 * y + c4 * x * y + c5 * x ** 2 + c6 * y ** 2 + c7 * x ** 2 * y + c8 * x * y ** 2 \
        + c9 * x ** 3 + c10 * y ** 3

    v = c1 + c2 * x + c3 * y + c4 * x * y + c5 * x ** 2 + c6 * y ** 2 + c7 * x ** 2 * y + c8 * x * y ** 2 \
        + c9 * x ** 3 + c10 * y ** 3

    # linear system of 10 equations

    # eq01 = p.subs(x, 0).subs(y, h)
    # eq02 = p.subs(x, h).subs(y, h)
    # eq03 = p.subs(x, h).subs(y, 0)
    # eq04 = p.subs(x, 0).subs(y, 0)
    #
    # eq05 = p.diff(x).subs(x, 0).subs(y, h)
    # eq06 = p.diff(x).subs(x, h).subs(y, h)
    # eq07 = p.diff(x).subs(x, h).subs(y, 0)
    #
    # eq08 = p.diff(y).subs(x, 0).subs(y, h)
    # eq09 = p.diff(y).subs(x, h).subs(y, h)
    # eq10 = p.diff(y).subs(x, h).subs(y, 0)

    eq01 = p.subs(x, -h).subs(y, 0)
    eq02 = p.subs(x, 0).subs(y, 0)
    eq03 = p.subs(x, 0).subs(y, -h)
    eq04 = p.subs(x, -h).subs(y, -h)

    eq05 = p.diff(x).subs(x, -h).subs(y, 0)
    eq06 = p.diff(x).subs(x, 0).subs(y, 0)
    eq07 = p.diff(x).subs(x, 0).subs(y, -h)

    eq08 = p.diff(y).subs(x, -h).subs(y, 0)
    eq09 = p.diff(y).subs(x, 0).subs(y, 0)
    eq10 = p.diff(y).subs(x, 0).subs(y, -h)

    sol = linsolve([eq01 - phi_i1_j, eq02 - phi_i_j, eq03 - phi_i_j1, eq04 - phi_i1_j1, eq05 - r_i1_j, eq06 - r_i_j,
                    eq07 - r_i_j1, eq08 - s_i1_j, eq09 - s_i_j, eq10 - s_i_j1],
                   (c1, c2, c3, c4, c5, c6, c7, c8, c9, c10))

    sol_get = next(iter(sol))

    c1_expr = sol_get[0] #.subs(r_i_j1, r_i_j1 / h).subs(r_i1_j1, r_i1_j1 / h).subs(r_i1_j, r_i1_j / h)\
       # .subs(s_i_j1, s_i_j1 / h).subs(s_i1_j1, s_i1_j1 / h).subs(s_i1_j, s_i1_j / h)

    c2_expr = sol_get[1] #.subs(r_i_j1, r_i_j1 / h).subs(r_i1_j1, r_i1_j1 / h).subs(r_i1_j, r_i1_j / h)\
       # .subs(s_i_j1, s_i_j1 / h).subs(s_i1_j1, s_i1_j1 / h).subs(s_i1_j, s_i1_j / h)

    c3_expr = sol_get[2] #.subs(r_i_j1, r_i_j1 / h).subs(r_i1_j1, r_i1_j1 / h).subs(r_i1_j, r_i1_j / h)\
       # .subs(s_i_j1, s_i_j1 / h).subs(s_i1_j1, s_i1_j1 / h).subs(s_i1_j, s_i1_j / h)

    c4_expr = sol_get[3] #.subs(r_i_j1, r_i_j1 / h).subs(r_i1_j1, r_i1_j1 / h).subs(r_i1_j, r_i1_j / h)\
       # .subs(s_i_j1, s_i_j1 / h).subs(s_i1_j1, s_i1_j1 / h).subs(s_i1_j, s_i1_j / h)

    c5_expr = sol_get[4] #.subs(r_i_j1, r_i_j1 / h).subs(r_i1_j1, r_i1_j1 / h).subs(r_i1_j, r_i1_j / h)\
       # .subs(s_i_j1, s_i_j1 / h).subs(s_i1_j1, s_i1_j1 / h).subs(s_i1_j, s_i1_j / h)

    c6_expr = sol_get[5] #.subs(r_i_j1, r_i_j1 / h).subs(r_i1_j1, r_i1_j1 / h).subs(r_i1_j, r_i1_j / h)\
        #.subs(s_i_j1, s_i_j1 / h).subs(s_i1_j1, s_i1_j1 / h).subs(s_i1_j, s_i1_j / h)

    c7_expr = sol_get[6] #.subs(r_i_j1, r_i_j1 / h).subs(r_i1_j1, r_i1_j1 / h).subs(r_i1_j, r_i1_j / h)\
       # .subs(s_i_j1, s_i_j1 / h).subs(s_i1_j1, s_i1_j1 / h).subs(s_i1_j, s_i1_j / h)

    c8_expr = sol_get[7] #.subs(r_i_j1, r_i_j1 / h).subs(r_i1_j1, r_i1_j1 / h).subs(r_i1_j, r_i1_j / h)\
        #.subs(s_i_j1, s_i_j1 / h).subs(s_i1_j1, s_i1_j1 / h).subs(s_i1_j, s_i1_j / h)

    c9_expr = sol_get[8] #.subs(r_i_j1, r_i_j1 / h).subs(r_i1_j1, r_i1_j1 / h).subs(r_i1_j, r_i1_j / h)\
        #.subs(s_i_j1, s_i_j1 / h).subs(s_i1_j1, s_i1_j1 / h).subs(s_i1_j, s_i1_j / h)

    c10_expr = sol_get[9] #.subs(r_i_j1, r_i_j1 / h).subs(r_i1_j1, r_i1_j1 / h).subs(r_i1_j, r_i1_j / h)\
        #.subs(s_i_j1, s_i_j1 / h).subs(s_i1_j1, s_i1_j1 / h).subs(s_i1_j, s_i1_j / h)

    print("c1 = {0}".format(c1_expr))
    print("c2 = {0}".format(c2_expr))
    print("c3 = {0}".format(c3_expr))
    print("c4 = {0}".format(c4_expr))
    print("c5 = {0}".format(c5_expr))
    print("c6 = {0}".format(c6_expr))
    print("c7 = {0}".format(c7_expr))
    print("c8 = {0}".format(c8_expr))
    print("c9 = {0}".format(c9_expr))
    print("c10 = {0}\n".format(c10_expr))

    # first derivative in the x-direction

    rp = p.diff(x)
    ru = u.diff(x)
    rv = v.diff(x)

    # first derivative in the y-direction

    sp = p.diff(y)
    su = u.diff(y)
    sv = v.diff(y)

    # second and third derivative for p

    p_xx = rp.diff(x)
    p_xy = rp.diff(y)
    p_yy = sp.diff(y)

    p_3x = p_xx.diff(x)
    p_xxy = p_xx.diff(y)
    p_xyy = p_xy.diff(y)
    p_3y = p_yy.diff(y)

    # second and third derivative for u

    u_xx = ru.diff(x)
    u_xy = ru.diff(y)
    u_yy = su.diff(y)

    u_3x = u_xx.diff(x)
    u_xxy = u_xx.diff(y)
    u_xyy = u_xy.diff(y)
    u_3y = u_yy.diff(y)

    # second and third derivative for v

    v_xx = rv.diff(x)
    v_xy = rv.diff(y)
    v_yy = sv.diff(y)

    v_3x = v_xx.diff(x)
    v_xxy = v_xx.diff(y)
    v_xyy = v_xy.diff(y)
    v_3y = v_yy.diff(y)

    # parameters values at the discrete point i

    pi = p.subs(x, 0).subs(y, 0)
    ui = u.subs(x, 0).subs(y, 0)
    vi = v.subs(x, 0).subs(y, 0)

    rpi = rp.subs(x, 0).subs(y, 0)  # first derivatives in the x-direction
    rui = ru.subs(x, 0).subs(y, 0)
    rvi = rv.subs(x, 0).subs(y, 0)

    spi = sp.subs(x, 0).subs(y, 0)  # first derivatives in the y-direction
    sui = su.subs(x, 0).subs(y, 0)
    svi = sv.subs(x, 0).subs(y, 0)

    p_xxi = p_xx.subs(x, 0).subs(y, 0)  # second derivatives
    print("p_xx = {0}\n".format(p_xxi.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr).subs(c4, c4_expr).subs(c5, c5_expr)\
       .subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(c9, c9_expr).subs(c10, c10_expr)))
    p_xyi = p_xy.subs(x, 0).subs(y, 0)
    print("p_xy = {0}\n".format(
        p_xyi.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr).subs(c4, c4_expr).subs(c5, c5_expr) \
        .subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(c9, c9_expr).subs(c10, c10_expr)))
    p_yyi = p_yy.subs(x, 0).subs(y, 0)

    u_xxi = u_xx.subs(x, 0).subs(y, 0)
    u_xyi = u_xy.subs(x, 0).subs(y, 0)
    u_yyi = u_yy.subs(x, 0).subs(y, 0)

    v_xxi = v_xx.subs(x, 0).subs(y, 0)
    v_xyi = v_xy.subs(x, 0).subs(y, 0)
    v_yyi = v_yy.subs(x, 0).subs(y, 0)

    p_3xi = p_3x.subs(x, 0).subs(y, 0)  # third derivatives

    p_xxyi = p_xxy.subs(x, 0).subs(y, 0)
    p_xyyi = p_xyy.subs(x, 0).subs(y, 0)
    p_3yi = p_3y.subs(x, 0).subs(y, 0)

    u_3xi = u_3x.subs(x, 0).subs(y, 0)
    u_xxyi = u_xxy.subs(x, 0).subs(y, 0)
    u_xyyi = u_xyy.subs(x, 0).subs(y, 0)
    u_3yi = u_3y.subs(x, 0).subs(y, 0)

    v_3xi = v_3x.subs(x, 0).subs(y, 0)
    v_xxyi = v_xxy.subs(x, 0).subs(y, 0)
    v_xyyi = v_xyy.subs(x, 0).subs(y, 0)
    v_3yi = v_3y.subs(x, 0).subs(y, 0)

    pni = pi - tau * (rui + svi) + tau**2 / 2 * (p_xxi + p_yyi) - tau**3 / 6 * (u_3xi + u_xyyi + v_xxyi + v_3yi)

    uni = ui - tau * rpi + tau**2 / 2 * (u_xxi + v_xyi) - tau**3 / 6 * (p_3xi + p_xyyi)

    vni = vi - tau * spi + tau**2 / 2 * (u_xyi + v_yyi) - tau**3 / 6 * (p_xxyi + p_3yi)

    # derivatives at the x-direction

    rpni = rpi - tau * (u_xxi + v_xyi) + tau**2 / 2 * (p_3xi + p_xyyi)

    runi = rui - tau * p_xxi + tau ** 2 / 2 * (u_3xi + v_xxyi)

    rvni = rvi - tau * p_xyi + tau**2 / 2 * (u_xxyi + v_xyyi)

    # derivatives at the y-direction

    spni = spi - tau * (u_xyi + v_yyi) + tau**2 / 2 * (p_xxyi + p_3yi)

    suni = sui - tau * p_xyi + tau**2 / 2 * (u_xxyi + v_xyyi)

    svni = svi - tau * p_yyi + tau**2 / 2 * (u_xyyi + v_3yi)

    _pni = pni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr).subs(c4, c4_expr).subs(c5, c5_expr)\
       .subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(c9, c9_expr).subs(c10, c10_expr)

    _uni = uni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr).subs(c4, c4_expr).subs(c5, c5_expr)\
       .subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(c9, c9_expr).subs(c10, c10_expr).subs(tau / h, nu)

    _vni = vni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr).subs(c4, c4_expr).subs(c5, c5_expr)\
        .subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(c9, c9_expr).subs(c10, c10_expr).subs(tau / h, nu)

    _rpni = rpni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr).subs(c4, c4_expr).subs(c5, c5_expr)\
       .subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(c9, c9_expr).subs(c10, c10_expr).subs(tau / h, nu)

    _runi = runi.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr).subs(c4, c4_expr).subs(c5, c5_expr)\
       .subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(c9, c9_expr).subs(c10, c10_expr).subs(tau / h, nu)

    _rvni = rvni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr).subs(c4, c4_expr).subs(c5, c5_expr)\
        .subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(c9, c9_expr).subs(c10, c10_expr).subs(tau / h, nu)

    _spni = spni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr).subs(c4, c4_expr).subs(c5, c5_expr)\
       .subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(c9, c9_expr).subs(c10, c10_expr).subs(tau / h, nu)

    _suni = suni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr).subs(c4, c4_expr).subs(c5, c5_expr)\
       .subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(c9, c9_expr).subs(c10, c10_expr).subs(tau / h, nu)

    _svni = svni.subs(c1, c1_expr).subs(c2, c2_expr).subs(c3, c3_expr).subs(c4, c4_expr).subs(c5, c5_expr)\
       .subs(c6, c6_expr).subs(c7, c7_expr).subs(c8, c8_expr).subs(c9, c9_expr).subs(c10, c10_expr).subs(tau / h, nu)

    print(simplify(_pni))
    # print(simplify(_uni))
    # print(simplify(_vni))
    print(simplify(_rpni))
    # print(simplify(_runi))
    # print(simplify(_rvni))
    print(simplify(_spni))
    # print(simplify(_suni))
    # print(simplify(_svni))

symboldeltaP7()