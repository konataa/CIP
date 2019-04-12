from __future__ import unicode_literals
import numpy as np
import math
import test_functions as test
# import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
# import matplotlib
# matplotlib.rcParams['text.usetex'] = True
# matplotlib.rcParams['text.latex.unicode'] = True
#
# from matplotlib.colors import BoundaryNorm
# from matplotlib.ticker import MaxNLocator
#
# plt.rc('font',**{'family':'serif','serif':['Times']})


def abs_minimum(x, y):
    if (math.fabs(x) >= math.fabs(y)):
        return math.fabs(x)
    else:
        return math.fabs(y)


def minmod(x, y):
    if (x * y > 0):
        return np.sign(x) * abs_minimum(x, y)
    else:
        return 0


def upwind_scheme(u, nt, nx, r, h):
    error = np.zeros(3)
    print(error)
    for n in range(0, nt - 1):
        for i in range(1, nx):
            u[n + 1][i] = u[n][i] - r * (u[n][i] - u[n][i - 1])
        # BC
        u[n + 1][0] = u[n][nx - 1]
    for i in range(0, nx - 1):
        delta = (u[nt - 1][i] - u[0][i]) * h
        error[0] += math.fabs(delta)
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    error[1] = math.sqrt(error[1])
    return error

def fourth_order_polynominal_interpolation(u, ux, nt, nx, h, ksi, del_rho, rho):
    error = np.zeros(3)
    h2 = h ** 2
    h3 = h ** 3
    h4 = h ** 4
    h5 = h ** 5

    xi = ksi
    xi2 = ksi ** 2
    xi3 = ksi ** 3
    xi4 = ksi ** 4
    for n in range(0, nt - 1):
        for i in range(1, nx):
            a = 5 * ((ux[n][i] - ux[n][i - 1]) * h2 - 6 * (u[n][i] + u[n][i - 1]) * h + 12 * rho[n][i - 1]) / (2 * h5)
            b = 2 * (
            (3 * ux[n][i] - 2 * ux[n][i - 1]) * h2 - (16 * u[n][i] + 14 * u[n][i - 1]) * h + 30 * rho[n][i - 1]) / (h4)
            c = 3 * (
            (3 * ux[n][i] - ux[n][i - 1]) * h2 - 4 * (3 * u[n][i] + 2 * u[n][i - 1]) * h + 20 * rho[n][i - 1]) / (
                2 * h3)
            u[n + 1][i] = a * xi4 + b * xi3 + c * xi2 + ux[n][i] * xi + u[n][i]
            ux[n + 1][i] = 4 * a * xi3 + 3 * b * xi2 + 2 * c * xi + ux[n][i]
            del_rho[n + 1][0] = - (a * xi4 / 5 + b * xi3 / 4 + c * xi2 / 3 + ux[n][0] * xi / 2 + u[n][0]) * xi
            del_rho[n + 1][i] = - (a * xi4 / 5 + b * xi3 / 4 + c * xi2 / 3 + ux[n][i] * xi / 2 + u[n][i]) * xi
            rho[n + 1][i - 1] = rho[n][i - 1] + del_rho[n + 1][i - 1] - del_rho[n + 1][i]

        # BC
        u[n + 1][0] = u[n][nx - 1]
        ux[n + 1][0] = ux[n][nx - 1]
        #        rho[n + 1][0]   = rho[n][nx - 1]
    for i in range(0, nx - 1):
        delta = (u[nt - 1][i] - u[0][i]) * h
        error[0] += math.fabs(delta)
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    error[1] = math.sqrt(error[1])
    #print(error)
    return error

#        rho[n + 1][0]   = rho[n][nx - 1]


def second_order_polynominal_interpolation(u, ux, nt, nx, h, ksi, del_rho, rho):
    error = np.zeros(3)
    print(error)
    h2 = h ** 2
    h3 = h ** 3

    xi = ksi
    xi2 = ksi ** 2
    xi3 = ksi ** 3
    for n in range(0, nt - 1):
        for i in range(1, nx):
            a = (h * u[n][i] + h * u[n][i - 1] - 2 * rho[n][i - 1]) / h3
            b = (2 * h * u[n][i] + h * u[n][i - 1] - 3 * rho[n][i - 1]) / h2
            u[n + 1][i] = a * xi2  + b * xi + u[n][i]
            del_rho[n + 1][0] = -( a * xi3 + b * xi2 + u[n][0] * xi)
            del_rho[n + 1][i] = -( a * xi3 + b * xi2 + u[n][i] * xi)
            rho[n + 1][i - 1] = rho[n][i - 1] + del_rho[n+1][i - 1] - del_rho[n+1][i]

        # BC
        u[n + 1][0] = u[n][nx - 1]
        ux[n + 1][0]    = ux[n][nx - 1]
    for i in range(0, nx - 1):
        delta = (u[nt - 1][i] - u[0][i]) * h
        error[0] += math.fabs(delta)
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    error[1] = math.sqrt(error[1])
    #print(error)
    return error


def rcip(x, nt, nx, h, xi, error):
    ux = np.zeros((nt, nx))
    u = np.zeros((nt, nx))
    u[0] = [np.sin(x[i]) for i in range(0, nx)]
    ux[0] = [np.cos(x[i]) for i in range(0, nx)]
    for n in range(0, nt - 1):
        print('N = {0}, h = {1}'.format(n, h))
        for i in range(0, nx - 1):
            s_i = (u[n][i + 1] - u[n][i]) / h
            b_i = ((s_i - ux[n][i]) / (ux[n][i + 1] - s_i) - 1) / h

            if (ux[n][i] * ux[n][i + 1] < 0):
                r1 = ux[n][i] + u[n][i] * b_i
                r2 = s_i * b_i + (s_i - ux[n][i]) / h

                u[n + 1][i] = (u[n][i] + r1*xi + r2*xi**2) / (1 + b_i*xi)
                ux[n + 1][i] = ((r1 + 2*r2*xi)*(1 + b_i*xi) - b_i * (u[n][i] + r1*xi + r2*xi**2)) / (1 + b_i*xi)**2
            else:
                c1 = ux[n][i]
                c2 = -(2*ux[n][i] + ux[n][i + 1] - 3*s_i) / h
                c3 = (ux[n][i] + ux[n][i + 1] - 2*s_i) / h**2

                u[n + 1][i] = u[n][i] + c1*xi + c2*xi**2 + c3*xi**3
                ux[n + 1][i] = c1 + 2*c2*xi + 3*c3*xi**2
            print('s_i = {0}, b_i = {1}, f = {2}, d = {3}'.format(s_i, b_i, u[n][i], ux[n][i]))
        u[n + 1][0] = u[n][nx - 1]
        ux[n + 1][0] = ux[n][nx - 1]

    for i in range(0, nx):
        delta = u[nt - 1][i] - u[0][i]
        error[0] += math.fabs(delta) * h
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)

    error[0] = error[0] / nt
    error[1] = math.sqrt(error[1])
    return u

def cip3(x, tau, h, xi, tfin, error):

    t = 0
    nx = len(x)
    ux = np.zeros(nx)
    u = np.zeros(nx)
    uxn = np.zeros(nx)
    un = np.zeros(nx)

    ua = np.zeros(nx)

    u = [test.test_func_sin(x[i]) for i in range(0, nx)]
    ua = [test.test_func_sin(x[i]) for i in range(0, nx)]

    ux = [test.der_test_func_sin(x[i]) for i in range(0, nx)]

    # u.T = [0 for i in range(0, nt)]
    # u.T[nx - 1] = [0 for i in range(0, nt)]

    while t < tfin:
        # if t + tau > tfin:
        #     tau = tfin - t
        #     xi = -a*tau
    # for n in range(0, nt - 1):

        a = (ux[0] + ux[nx - 2]) / h ** 2 - 2 * (u[0] - u[nx - 2]) / h ** 3

        b = (2 * ux[0] + ux[nx - 2]) / h - 3 * (u[0] - u[nx - 2]) / h ** 2

        un[0] = u[0] + ux[0] * xi + b * xi ** 2 + a * xi ** 3

        uxn[0] = ux[0] + 2 * b * xi + 3 * a * xi ** 2

        for i in range(1, nx - 1):

            a = (ux[i] + ux[i - 1]) / h**2 - 2 * (u[i] - u[i - 1]) / h**3

            b = (2 * ux[i] + ux[i - 1]) / h - 3 * (u[i] - u[i - 1]) / h**2

            un[i] = u[i] + ux[i] * xi + b * xi**2 + a * xi**3

            uxn[i] = ux[i] + 2 * b * xi + 3 * a * xi**2

        u[:] = un[:]
        ux[:] = uxn[:]

        t = t + tau

    for i in range(0, nx):
        delta = u[i] - ua[i]
        # error[0] += math.fabs(delta) * h
        error[1] += delta ** 2
        if(math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    error[0] = h * np.sum(np.abs(np.subtract(ua, u)))
    error[1] = math.sqrt(error[1])
    return u, ua


def cip5(x, tau, h, xi, tfin, nt, error):
    # ux = np.zeros((nt, nx))
    # u = np.zeros((nt, nx))
    # u2x = np.zeros((nt, nx))
    # u[0] = [np.sin(x[i]) for i in range(0, nx)]
    # u.T[nx - 1] = [0 for i in range(0, nt)]

    # for n in range(0, nt - 1):
    #
    #     a = (u2x[n][0] * h ** 2 - u2x[n][nx - 2] * h ** 2 - 6 * ux[n][0] * h - 6 * ux[n][nx - 2] * h + 12 * u[n][
    #         0] - 12 * u[n][nx - 2]) / (2 * h ** 5)
    #     b = (3 * u2x[n][0] * h ** 2 / 2 - u2x[n][nx - 2] * h ** 2 - 8 * ux[n][0] * h - 7 * ux[n][nx - 2] * h + 15 *
    #          u[n][0] - 15 * u[n][nx - 2]) / h ** 4
    #     c = (3 * u2x[n][0] * h ** 2 - u2x[n][nx - 2] * h ** 2 - 12 * ux[n][0] * h - 8 * ux[n][nx - 2] * h + 20 * u[n][
    #         0] - 20 * u[n][nx - 2]) / (2 * h ** 3)
    #
    #     d = u2x[n][0] / 2
    #
    #     e = ux[n][0]
    #
    #     f = u[n][0]
    #
    #     u[n + 1][0] = f + e * xi + d * xi ** 2 + c * xi ** 3 + b * xi ** 4 + a * xi ** 5
    #
    #     ux[n + 1][0] = e + 2 * d * xi + 3 * c * xi ** 2 + 4 * b * xi ** 3 + 5 * a * xi ** 4
    #
    #     u2x[n + 1][0] = 2 * d + 6 * c * xi + 12 * b * xi ** 2 + 20 * a * xi ** 3
    #
    #     for i in range(1, nx - 1):
    #         # Main #
    #         a = (u2x[n][i]*h**2 - u2x[n][i - 1]*h**2 - 6*ux[n][i]*h - 6*ux[n][i - 1]*h + 12*u[n][i] - 12*u[n][i - 1])/(2*h**5)
    #         b = (3*u2x[n][i]*h**2/2 - u2x[n][i - 1]*h**2 - 8*ux[n][i]*h - 7*ux[n][i - 1]*h + 15*u[n][i] - 15*u[n][i - 1])/h**4
    #         c = (3*u2x[n][i]*h**2 - u2x[n][i - 1]*h**2 - 12*ux[n][i]*h - 8*ux[n][i - 1]*h + 20*u[n][i] - 20*u[n][i - 1])/(2*h**3)
    #         d = u2x[n][i] / 2
    #
    #         e = ux[n][i]
    #
    #         f = u[n][i]
    #
    #         # print("a = {0:3f}, b = {1:3f}, c = {2:3f}, d = {3:3f}, e = {4:3f}, f = {5:3f}".format(a,b,c,d,e,f))
    #
    #         u[n + 1][i] = f + e * xi + d * xi**2 + c * xi**3 + b * xi**4 + a * xi**5
    #
    #         ux[n + 1][i] = e + 2 * d * xi + 3 * c * xi**2 + 4 * b * xi**3 + 5 * a * xi**4
    #
    #         u2x[n + 1][i] = 2 * d + 6 * c * xi + 12 * b * xi**2 + 20 * a * xi**3

    # for i in range(0, nx):
    #     delta = u[nt - 1][i] - u[0][i]
    #     error[0] += math.fabs(delta) * h
    #     error[1] += delta ** 2
    #     if (math.fabs(delta) > error[2]):
    #         error[2] = math.fabs(delta)
    # error[0] = error[0] / nt
    # error[1] = math.sqrt(error[1])

    t = 0
    nx = len(x)
    u = np.zeros(nx)
    ux = np.zeros(nx)
    u2x = np.zeros(nx)

    un = np.zeros(nx)
    uxn = np.zeros(nx)
    u2xn = np.zeros(nx)

    ua = np.zeros(nx)

    u = [test.test_func_sin(x[i]) for i in range(0, nx)]
    ua = [test.test_func_sin(x[i]) for i in range(0, nx)]

    ux = [test.der_test_func_sin(x[i]) for i in range(0, nx)]
    u2x = [test.der2_test_func_sin(x[i]) for i in range(0, nx)]

    # u.T = [0 for i in range(0, nt)]
    # u.T[nx - 1] = [0 for i in range(0, nt)]

    while t < tfin:
    #     if t + tau > tfin:
    #         tau = tfin - t
    #         xi = -a * tau
    # for n in range(0, nt - 1):

        a = (u2x[0] * h ** 2 - u2x[nx - 2] * h ** 2 - 6 * ux[0] * h - 6 * ux[nx - 2] * h + 12 * u[0] - 12 * u[nx - 2]) / (2 * h ** 5)

        b = (3 * u2x[0] * h ** 2 / 2 - u2x[nx - 2] * h ** 2 - 8 * ux[0] * h - 7 * ux[nx - 2] * h + 15 * u[0] - 15 * u[nx - 2]) / h ** 4

        c = (3 * u2x[0] * h ** 2 - u2x[nx - 2] * h ** 2 - 12 * ux[0] * h - 8 * ux[nx - 2] * h + 20 * u[0] - 20 * u[nx - 2]) / (2 * h ** 3)

        d = u2x[0] / 2

        e = ux[0]

        f = u[0]

        un[0] = f + e * xi + d * xi ** 2 + c * xi ** 3 + b * xi ** 4 + a * xi ** 5

        uxn[0] = e + 2 * d * xi + 3 * c * xi ** 2 + 4 * b * xi ** 3 + 5 * a * xi ** 4

        u2xn[0] = 2 * d + 6 * c * xi + 12 * b * xi ** 2 + 20 * a * xi ** 3

        for i in range(1, nx - 1):
            a = (u2x[i] * h ** 2 - u2x[i - 1] * h ** 2 - 6 * ux[i] * h - 6 * ux[i - 1] * h + 12 * u[i] - 12 * u[i - 1]) / (2 * h ** 5)

            b = (3*u2x[i]*h**2/2 - u2x[i - 1]*h**2 - 8*ux[i]*h - 7*ux[i - 1]*h + 15*u[i] - 15*u[i - 1])/h**4

            c = (3*u2x[i]*h**2 - u2x[i - 1]*h**2 - 12*ux[i]*h - 8*ux[i - 1]*h + 20*u[i] - 20*u[i - 1])/(2*h**3)

            d = u2x[i] / 2

            e = ux[i]

            f = u[i]

            #         # print("a = {0:3f}, b = {1:3f}, c = {2:3f}, d = {3:3f}, e = {4:3f}, f = {5:3f}".format(a,b,c,d,e,f))

            un[i] = f + e * xi + d * xi**2 + c * xi**3 + b * xi**4 + a * xi**5

            uxn[i] = e + 2 * d * xi + 3 * c * xi**2 + 4 * b * xi**3 + 5 * a * xi**4

            u2xn[i] = 2 * d + 6 * c * xi + 12 * b * xi**2 + 20 * a * xi**3

        u[:] = un[:]
        ux[:] = uxn[:]
        u2x[:] = u2xn[:]


        t = t + tau

        # u[n + 1][0] = u[n][nx - 1]
        # ux[n + 1][0] = ux[n][nx - 1]
    # print(u)
    for i in range(0, nx):
        delta = u[i] - ua[i]
        # error[0] += math.fabs(delta) * h
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    error[0] = h * np.sum(np.abs(np.subtract(ua, u)))
    error[1] = math.sqrt(error[1])
    return u, ua

def cip7(x, tau, nt, h, xi, tfin, error):
    nx = len(x)
    u = np.zeros(nx)
    ux = np.zeros(nx)
    u2x = np.zeros(nx)
    u3x = np.zeros(nx)

    un = np.zeros(nx)
    uxn = np.zeros(nx)
    u2xn = np.zeros(nx)
    u3xn = np.zeros(nx)

    ua = np.zeros(nx)

    u = [test.test_func_sin(x[i]) for i in range(0, nx)]
    ua = [test.test_func_sin(x[i]) for i in range(0, nx)]

    ux = [test.der_test_func_sin(x[i]) for i in range(0, nx)]
    u2x = [test.der2_test_func_sin(x[i]) for i in range(0, nx)]
    u3x = [test.der3_test_func_sin(x[i]) for i in range(0, nx)]

    t, a = 0, 1

    # u[0] = [np.sin(x[i]) for i in range(0, nx)]
    # u.T[nx - 1] = [0 for i in range(0, nt)]
    while t < tfin:
        # if t + tau > tfin:
        #     tau = tfin - t
        #     xi = -a * tau
    # for n in range(0, nt - 1):

        a = (u3x[0] * h ** 3 + u3x[nx - 2] * h ** 3 - 12 * u2x[0] * h ** 2 + 12 * u2x[nx - 2] * h ** 2 + 60 *
             ux[0] * h + 60 * ux[nx - 2] * h - 120 * u[0] + 120 * u[nx - 2]) / (
                    6 * h ** 7)

        b = (4 * u3x[0] * h ** 3 + 3 * u3x[nx - 2] * h ** 3 - 45 * u2x[0] * h ** 2 + 39 * u2x[nx - 2] * h ** 2
              + 216 * ux[0] * h + 204 * ux[nx - 2] * h - 420 * u[0] + 420 * u[nx - 2]) / (6 * h ** 6)

        c = (u3x[0] * h ** 3 + u3x[nx - 2] * h ** 3 / 2 - 10 * u2x[0] * h ** 2 + 7 * u2x[nx - 2] * h ** 2
             + 45 * ux[0] * h + 39 * ux[nx - 2] * h - 84 * u[0] + 84 * u[nx - 2]) / h ** 5

        d = (4 * u3x[0] * h ** 3 + u3x[nx - 2] * h ** 3 - 30 * u2x[0] * h ** 2 + 15 * u2x[nx - 2] * h ** 2
             + 120 * ux[0] * h + 90 * ux[nx - 2] * h - 210 * u[0] + 210 * u[nx - 2]) / (6 * h ** 4)

        e = u3x[0] / 6

        f = u2x[0] / 2

        g = ux[0]

        l = u[0]

        un[0] = l + g * xi + f * xi ** 2 + e * xi ** 3 + d * xi ** 4 + c * xi ** 5 + b * xi ** 6 + a * xi ** 7
        uxn[0] = g + 2 * f * xi + 3 * e * xi ** 2 + 4 * d * xi ** 3 + 5 * c * xi ** 4 + 6 * b * xi ** 5 + 7 * a * xi ** 6
        u2xn[0] = 2 * f + 6 * e * xi + 12 * d * xi ** 2 + 20 * c * xi ** 3 + 30 * b * xi ** 4 + 42 * a * xi ** 5
        u3xn[0] = 6 * e + 24 * d * xi + 60 * c * xi ** 2 + 120 * b * xi ** 3 + 210 * a * xi ** 4

        for i in range(1, nx - 1):
            a = (u3x[i] * h ** 3 + u3x[i - 1] * h ** 3 - 12 * u2x[i] * h ** 2 + 12 * u2x[i - 1] * h ** 2
                 + 60 * ux[i] * h + 60 * ux[i - 1] * h - 120 * u[i] + 120 * u[i - 1]) / (6 * h ** 7)

            b = (4 * u3x[i] * h ** 3 + 3 * u3x[i - 1] * h ** 3 - 45 * u2x[i] * h ** 2 + 39 * u2x[i - 1] * h ** 2
                 + 216 * ux[i] * h + 204 * ux[i - 1] * h - 420 * u[i] + 420 * u[i - 1]) / (6 * h ** 6)

            c = (u3x[i] * h ** 3 + u3x[i - 1] * h ** 3 / 2 - 10 * u2x[i] * h ** 2 + 7 * u2x[i - 1] * h ** 2
                 + 45 * ux[i] * h + 39 * ux[i - 1] * h - 84 * u[i] + 84 * u[i - 1]) / h ** 5

            d = (4 * u3x[i] * h ** 3 + u3x[i - 1] * h ** 3 - 30 * u2x[i] * h ** 2
                 + 15 * u2x[i - 1] * h ** 2 + 120 * ux[i] * h + 90 * ux[i - 1] * h - 210 * u[i] + 210 * u[i - 1]) / (
                    6 * h ** 4)

            e = u3x[i] / 6

            f = u2x[i] / 2

            g = ux[i]

            l = u[i]

            un[i] = l + g * xi + f * xi**2 + e * xi**3 + d * xi**4 + c * xi**5 + b * xi**6 + a * xi**7
            uxn[i] = g + 2 * f * xi + 3 * e * xi**2 + 4 * d * xi**3 + 5 * c * xi**4 + 6 * b * xi**5 + 7 * a * xi**6
            u2xn[i] = 2 * f + 6 * e * xi + 12 * d * xi**2 + 20 * c * xi**3 + 30 * b * xi**4 + 42 * a * xi**5
            u3xn[i] = 6 * e + 24 * d * xi + 60 * c * xi**2 + 120 * b * xi**3 + 210 * a * xi**4
        u[:] = un[:]
        ux[:] = uxn[:]
        u2x[:] = u2xn[:]
        u3x[:] = u3xn[:]

        t = t + tau

    for i in range(0, nx):
        delta = u[i] - ua[i]
        # error[0] += math.fabs(delta)*h
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    error[0] = h * np.sum(np.abs(np.subtract(ua, u)))
    error[1] = math.sqrt(error[1])
    return u


def cip9(x, tau, nt, h, xi, tfin, error):
    nx = len(x)

    u = np.zeros(nx)
    ux = np.zeros(nx)
    u2x = np.zeros(nx)
    u3x = np.zeros(nx)
    u4x = np.zeros(nx)

    un = np.zeros(nx)
    uxn = np.zeros(nx)
    u2xn = np.zeros(nx)
    u3xn = np.zeros(nx)
    u4xn = np.zeros(nx)

    ua = np.zeros(nx)

    u = [test.test_func_sin(x[i]) for i in range(0, nx)]
    ua = [test.test_func_sin(x[i]) for i in range(0, nx)]

    ux = [test.der_test_func_sin(x[i]) for i in range(0, nx)]
    u2x = [test.der2_test_func_sin(x[i]) for i in range(0, nx)]
    u3x = [test.der3_test_func_sin(x[i]) for i in range(0, nx)]
    u4x = [test.der4_test_func_sin(x[i]) for i in range(0, nx)]

    t, a = 0, 1

    while t < tfin:
        # if t + tau > tfin:
        #     tau = tfin - t
        #     xi = -a * tau
        a = (u4x[0] * h ** 4 - u4x[nx - 2] * h ** 4 - 20 * u3x[0] * h ** 3 - 20 * u3x[nx - 2] * h ** 3
             + 180 * u2x[0] * h ** 2 - 180 * u2x[nx - 2] * h ** 2 - 840 * ux[0] * h - 840 * ux[nx - 2] * h
             + 1680 * u[0] - 1680 * u[nx - 2]) / (24 * h ** 9)

        b = (5 * u4x[0] * h ** 4 - 4 * u4x[nx - 2] * h ** 4 - 96 * u3x[0] * h ** 3 - 84 * u3x[nx - 2] * h ** 3
             + 840 * u2x[0] * h ** 2 - 780 * u2x[nx - 2] * h ** 2 - 3840 * ux[0] * h - 3720 * ux[nx - 2] * h
             + 7560 * u[0] - 7560 * u[nx - 2]) / (24 * h ** 8)

        c =(5 * u4x[0] * h ** 4 - 3 * u4x[nx - 2] * h ** 4 - 90 * u3x[0] * h ** 3 - 66 * u3x[nx - 2] * h ** 3
            + 756 * u2x[0] * h ** 2 - 636 * u2x[nx - 2] * h ** 2 - 3360 * ux[0] * h - 3120 * ux[nx - 2] * h
            + 6480 * u[0] - 6480 * u[nx - 2]) / (12 * h ** 7)

        d = (5 * u4x[0] * h ** 4 - 2 * u4x[nx - 2] * h ** 4 - 80 * u3x[0] * h ** 3 - 46 * u3x[nx - 2] * h ** 3
             + 630 * u2x[0] * h ** 2 - 462 * u2x[nx - 2] * h ** 2 - 2688 * ux[0] * h - 2352 * ux[nx - 2] * h
             + 5040 * u[0] - 5040 * u[nx - 2]) / (12 * h ** 6)

        e = (5 * u4x[0] * h ** 4 - u4x[nx - 2] * h ** 4 - 60 * u3x[0] * h ** 3 - 24 * u3x[nx - 2] * h ** 3
             + 420 * u2x[0] * h ** 2 - 252 * u2x[nx - 2] * h ** 2 - 1680 * ux[0] * h - 1344 * ux[nx - 2] * h
             + 3024 * u[0] - 3024 * u[nx - 2]) / (24 * h ** 5)

        f = u4x[0] / 24

        g = u3x[0] / 6

        l = u2x[0] / 2

        k = ux[0]

        q = u[0]

        un[0] = q + k * xi + l * xi**2 + g * xi**3 + f * xi**4 + e * xi**5 + d * xi**6 + c * xi**7 + b * xi**8 + a * xi**9
        uxn[0] = k + 2 * l * xi + 3 * g * xi**2 + 4 * f * xi**3 + 5 * e * xi**4 + 6 * d * xi**5 + 7 * c * xi**6 + 8 * b * xi**7 + 9 * a * xi**8
        u2xn[0] = 2 * l + 6 * g * xi + 12 * f * xi**2 + 20 * e * xi**3 + 30 * d * xi**4 + 42 * c * xi**5 + 56 * b * xi**6 + 72 * a * xi**7
        u3xn[0] = 6 * g + 24 * f * xi + 60 * e * xi**2 + 120 * d * xi**3 + 210 * c * xi**4 + 336 * b * xi**5 + 504 * a * xi**6
        u4xn[0] = 24 * f + 120 * e * xi + 360 * d * xi**2 + 840 * c * xi**3 + 1680 * b * xi**4 + 3024 * a * xi**5

        for i in range(1, nx - 1):
            a = (u4x[i] * h ** 4 - u4x[i - 1] * h ** 4 - 20 * u3x[i] * h ** 3 - 20 * u3x[i - 1] * h ** 3
                 + 180 * u2x[i] * h ** 2 - 180 * u2x[i - 1] * h ** 2 - 840 * ux[i] * h - 840 * ux[i - 1] * h
                 + 1680 * u[i] - 1680 * u[i - 1]) / (24 * h ** 9)

            b = (5 * u4x[i] * h ** 4 - 4 * u4x[i - 1] * h ** 4 - 96 * u3x[i] * h ** 3 - 84 * u3x[i - 1] * h ** 3
                 + 840 * u2x[i] * h ** 2 - 780 * u2x[i - 1] * h ** 2 - 3840 * ux[i] * h - 3720 * ux[i - 1] * h
                 + 7560 * u[i] - 7560 * u[i - 1]) / (24 * h ** 8)

            c = (5 * u4x[i] * h ** 4 - 3 * u4x[i - 1] * h ** 4 - 90 * u3x[i] * h ** 3 - 66 * u3x[i - 1] * h ** 3
                 + 756 * u2x[i] * h ** 2 - 636 * u2x[i - 1] * h ** 2 - 3360 * ux[i] * h - 3120 * ux[i - 1] * h
                 + 6480 * u[i] - 6480 * u[i - 1]) / (12 * h ** 7)

            d = (5 * u4x[i] * h ** 4 - 2 * u4x[i - 1] * h ** 4 - 80 * u3x[i] * h ** 3 - 46 * u3x[i - 1] * h ** 3
                 + 630 * u2x[i] * h ** 2 - 462 * u2x[i - 1] * h ** 2 - 2688 * ux[i] * h - 2352 * ux[i - 1] * h
                 + 5040 * u[i] - 5040 * u[i - 1]) / (12 * h ** 6)

            e = (5 * u4x[i] * h ** 4 - u4x[i - 1] * h ** 4 - 60 * u3x[i] * h ** 3 - 24 * u3x[i - 1] * h ** 3
                 + 420 * u2x[i] * h ** 2 - 252 * u2x[i - 1] * h ** 2 - 1680 * ux[i] * h - 1344 * ux[i - 1] * h
                 + 3024 * u[i] - 3024 * u[i - 1]) / (24 * h ** 5)

            f = u4x[i] / 24

            g = u3x[i] / 6

            l = u2x[i] / 2

            k = ux[i]

            q = u[i]

            un[i] = q + k * xi + l * xi ** 2 + g * xi ** 3 + f * xi ** 4 + e * xi ** 5 + d * xi ** 6 + c * xi ** 7 + b * xi ** 8 + a * xi ** 9
            uxn[i] = k + 2 * l * xi + 3 * g * xi ** 2 + 4 * f * xi ** 3 + 5 * e * xi ** 4 + 6 * d * xi ** 5 + 7 * c * xi ** 6 + 8 * b * xi ** 7 + 9 * a * xi ** 8
            u2xn[i] = 2 * l + 6 * g * xi + 12 * f * xi ** 2 + 20 * e * xi ** 3 + 30 * d * xi ** 4 + 42 * c * xi ** 5 + 56 * b * xi ** 6 + 72 * a * xi ** 7
            u3xn[i] = 6 * g + 24 * f * xi + 60 * e * xi ** 2 + 120 * d * xi ** 3 + 210 * c * xi ** 4 + 336 * b * xi ** 5 + 504 * a * xi ** 6
            u4xn[i] = 24 * f + 120 * e * xi + 360 * d * xi ** 2 + 840 * c * xi ** 3 + 1680 * b * xi ** 4 + 3024 * a * xi ** 5

        u[:] = un[:]
        ux[:] = uxn[:]
        u2x[:] = u2xn[:]
        u3x[:] = u3xn[:]
        u4x[:] = u4xn[:]

        t = t + tau

    for i in range(0, nx):
        delta = u[i] - ua[i]
        # error[0] += math.fabs(delta)*h
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    error[0] = h * np.sum(np.abs(np.subtract(ua, u)))
    error[1] = math.sqrt(error[1])

    return u


def deltaP3(x, nu, tau, h, tfin, error):
    nx = len(x)
    u = np.zeros(nx)
    r = np.zeros(nx)

    un = np.zeros(nx)
    rn = np.zeros(nx)

    u = test.test_func_sin(x)
    ua = test.test_func_sin(x)

    r = test.der_test_func_sin(x) * h

    t = 0

    while(t < tfin):
        un[0] = (nu - 1) ** 2 * (2 * nu + 1) * u[0] + nu ** 2 * (3 - 2 * nu) * u[nx - 2] + nu ** 2 * (1 - nu) * r[nx - 2] - \
                nu * (nu - 1) ** 2 * r[0]

        rn[0] = 6 * nu * (1 - nu) * (u[0] - u[nx - 2]) + (nu - 1) * (3 * nu - 1) * r[0] + nu * (3 * nu - 2) * r[nx - 2]

        for i in range(1, nx - 1):
            un[i] = (nu - 1)**2 * (2*nu + 1) * u[i] + nu**2 * (3 - 2 * nu) * u[i - 1] + nu**2 * (1 - nu) * r[i - 1] - \
                    nu * (nu - 1)**2 * r[i]

            rn[i] = 6 * nu * (1 - nu) * (u[i] - u[i - 1]) + (nu - 1) * (3 * nu - 1) * r[i] + nu * (3 * nu - 2) * r[i - 1]

        u[:] = un[:]
        r[:] = rn[:]

        t = t + tau

    for i in range(0, nx):
        delta = u[i] - ua[i]
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    error[0] = h * np.sum(np.abs(np.subtract(ua, u)))
    error[1] = math.sqrt(error[1])

    return u, ua


def deltaP3sym(x, nu, tau, h, tfin, error):
    nx = len(x)
    u = np.zeros(nx)
    r = np.zeros(nx)
    a = 1
    un = np.zeros(nx)
    rn = np.zeros(nx)
    nu = tau / h
    u = test.test_func_sin(x)
    ua = test.test_func_sin(x)

    r = test.der_test_func_sin(x) * h

    t = 0
    while(t < tfin):

        # c1 = u[0]
        # c2 = r[0] / h
        #
        # c3 = (3 * (u[nx - 2] - u[0]) + 2 * r[0] + r[nx - 2]) / h ** 2
        # c4 = (2 * (u[nx - 2] - u[0]) + r[0] + r[nx - 2]) / h ** 3

        # un[0] = -a ** 3 * c4 * tau ** 3 + a ** 2 * c3 * tau ** 2 - a * c2 * tau + c1
        # rn[0] = h * (3 * a ** 2 * c4 * tau ** 2 - 2 * a * c3 * tau + c2)
        un[0] = nu ** 3 * (-r[0] - r[nx - 2] + 2 * u[0] - 2 * u[nx - 2]) + nu ** 2 \
                * (2 * r[0] + r[nx - 2] - 3 * u[0] + 3 * u[nx - 2]) - nu * r[0] + u[0]

        rn[0] = 3*nu**2*(r[0] + r[nx - 2] - 2*u[0] + 2*u[nx - 2]) \
                - 2*nu*(2*r[0] + r[nx - 2] - 3*u[0] + 3*u[nx - 2]) + r[0]

        for i in range(1, nx - 1):
            # c1 = u[i]
            # c2 = r[i] / h
            #
            # c3 = (3 * (u[i - 1] - u[i]) + 2 * r[i] + r[i - 1]) / h ** 2
            # c4 = (2 * (u[i - 1] - u[i]) + r[i] + r[i - 1]) / h ** 3

            # un[i] = -a ** 3 * c4 * tau**3 + a**2 * c3 * tau ** 2 - a * c2 * tau + c1
            # rn[i] = h * (3 * a**2 * c4 * tau**2 - 2 * a * c3 * tau + c2)

            un[i] = nu**3*(-r[i] - r[i - 1] + 2*u[i] - 2*u[i - 1]) \
                    + nu**2*(2*r[i] + r[i - 1] - 3*u[i] + 3*u[i - 1]) - nu*r[i] + u[i]

            rn[i] = 3*nu**2*(r[i] + r[i - 1] - 2*u[i] + 2*u[i - 1]) \
                    - 2*nu*(2*r[i] + r[i - 1] - 3*u[i] + 3*u[i - 1]) + r[i]


        u[:] = un[:]
        r[:] = rn[:]

        t = t + tau
    for i in range(0, nx):
        delta = u[i] - ua[i]
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    error[0] = h * np.sum(np.abs(np.subtract(ua, u)))
    error[1] = math.sqrt(error[1])

    return u, ua

def deltaP5sym(x, nu, tau, h, tfin, error):
    nx = len(x)
    u = np.zeros(nx)
    r = np.zeros(nx)
    a = 1
    un = np.zeros(nx)
    rn = np.zeros(nx)
    rxn = np.zeros(nx)
    nu = tau / h

    u = test.test_func_sin(x)
    ua = test.test_func_sin(x)

    r = test.der_test_func_sin(x) * h
    rx = test.der2_test_func_sin(x) * h**2

    t = 0
    while(t < tfin):
        ui = u[0]
        ui1 = u[nx - 2]
        ri = r[0]
        ri1 = r[nx - 2]
        rxi = rx[0]
        rxi1 = rx[nx - 2]

        # c6 = (rx[0] * h - rx[nx - 2] * h - 6 * r[0] - 6 * r[nx - 2] + 12 * u[0] - 12 * u[
        #     nx - 2]) / (2 * h ** 5)
        #
        # c5 = (3 * rx[0] * h / 2 - rx[nx - 2] * h - 8 * r[0] - 7 * r[nx - 2] + 15 * u[0] - 15 * u[
        #     nx - 2]) / h ** 4
        #
        # c4 = (3 * rx[0] * h - rx[nx - 2] * h - 12 * r[0] - 8 * r[nx - 2] + 20 * u[0] - 20 * u[
        #     nx - 2]) / (2 * h ** 3)
        #
        # c3 = rx[0] / (2 * h)
        #
        # c2 = r[0] / h
        #
        # c1 = u[0]

        un[0] = nu**5*(6*ri + 6*ri1 - rxi + rxi1 - 12*ui + 12*ui1)/2 + nu**4*(-16*ri - 14*ri1 + 3*rxi - 2*rxi1 + 30*ui - 30*ui1)/2 + nu**3*(12*ri + 8*ri1 - 3*rxi + rxi1 - 20*ui + 20*ui1)/2 + nu**2*rxi/2 - nu*ri + ui

        rn[0] = 5*nu**4*(-6*ri - 6*ri1 + rxi - rxi1 + 12*ui - 12*ui1)/2 + 2*nu**3*(16*ri + 14*ri1 - 3*rxi + 2*rxi1 - 30*ui + 30*ui1) + 3*nu**2*(-12*ri - 8*ri1 + 3*rxi - rxi1 + 20*ui - 20*ui1)/2 - nu*rxi + ri

        rxn[0] = 10*nu**3*(6*ri + 6*ri1 - rxi + rxi1 - 12*ui + 12*ui1) + 6*nu**2*(-16*ri - 14*ri1 + 3*rxi - 2*rxi1 + 30*ui - 30*ui1) + 3*nu*(12*ri + 8*ri1 - 3*rxi + rxi1 - 20*ui + 20*ui1) + rxi


        for i in range(1, nx - 1):
            ui = u[i]
            ui1 = u[i - 1]
            ri = r[i]
            ri1 = r[i - 1]
            rxi = rx[i]
            rxi1 = rx[i - 1]
            # c6 = (rx[i] * h - rx[i - 1] * h - 6 * r[i] - 6 * r[i - 1] + 12 * u[i] - 12 * u[
            #     i - 1]) / (2 * h ** 5)
            #
            # c5 = (3 * rx[i] * h / 2 - rx[i - 1] * h - 8 * r[i] - 7 * r[i - 1] + 15 * u[i] - 15 * u[
            #     i - 1]) / h ** 4
            #
            # c4 = (3 * rx[i] * h - rx[i - 1] * h - 12 * r[i] - 8 * r[i - 1] + 20 * u[i] - 20 * u[
            #     i - 1]) / (2 * h ** 3)
            #
            # c3 = rx[i] / (2 * h)
            #
            # c2 = r[i] / h
            #
            # c1 = u[i]

            un[i] = nu**5*(6*ri + 6*ri1 - rxi + rxi1 - 12*ui + 12*ui1)/2 + nu**4*(-16*ri - 14*ri1 + 3*rxi - 2*rxi1 + 30*ui - 30*ui1)/2 + nu**3*(12*ri + 8*ri1 - 3*rxi + rxi1 - 20*ui + 20*ui1)/2 + nu**2*rxi/2 - nu*ri + ui


            rn[i] = 5*nu**4*(-6*ri - 6*ri1 + rxi - rxi1 + 12*ui - 12*ui1)/2 + 2*nu**3*(16*ri + 14*ri1 - 3*rxi + 2*rxi1 - 30*ui + 30*ui1) + 3*nu**2*(-12*ri - 8*ri1 + 3*rxi - rxi1 + 20*ui - 20*ui1)/2 - nu*rxi + ri


            rxn[i] = 10*nu**3*(6*ri + 6*ri1 - rxi + rxi1 - 12*ui + 12*ui1) + 6*nu**2*(-16*ri - 14*ri1 + 3*rxi - 2*rxi1 + 30*ui - 30*ui1) + 3*nu*(12*ri + 8*ri1 - 3*rxi + rxi1 - 20*ui + 20*ui1) + rxi



        u[:] = un[:]
        r[:] = rn[:]
        rx[:] = rxn[:]

        t = t + tau
    for i in range(0, nx):
        delta = u[i] - ua[i]
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    error[0] = h * np.sum(np.abs(np.subtract(ua, u)))
    error[1] = math.sqrt(error[1])

    return u, ua


def deltaP7sym(x, nu, tau, h, tfin, error):
    nx = len(x)
    u = np.zeros(nx)
    r = np.zeros(nx)
    a = 1
    un = np.zeros(nx)
    rn = np.zeros(nx)
    rxn = np.zeros(nx)
    rxxn = np.zeros(nx)

    u = test.test_func_sin(x)
    ua = test.test_func_sin(x)

    r = test.der_test_func_sin(x) * h
    rx = test.der2_test_func_sin(x) * h ** 2
    rxx = test.der3_test_func_sin(x) * h ** 3

    t = 0
    while(t < tfin):
        ui = u[0]
        ui1 = u[nx - 2]
        ri = r[0]
        ri1 = r[nx - 2]
        rxi = rx[0]
        rxi1 = rx[nx - 2]
        rxxi = rxx[0]
        rxxi1 = rxx[nx - 2]

        # c8 = (rxx[0] * h ** 2 + rxx[nx - 2] * h ** 2 - 12 * rx[0] * h + 12 * rx[nx - 2] * h
        #          + 60 * r[0] + 60 * r[nx - 2] - 120 * u[0] + 120 * u[nx - 2]) / (6 * h ** 7)
        #
        # c7 = (4 * rxx[0] * h ** 2 + 3 * rxx[nx - 2] * h ** 2 - 45 * rx[0] * h + 39 * rx[nx - 2] * h
        #       + 216 * r[0] + 204 * r[nx - 2] - 420 * u[0] + 420 * u[nx - 2]) / (6 * h ** 6)
        #
        # c6 = (rxx[0] * h ** 2 + rxx[nx - 2] * h ** 2 / 2 - 10 * rx[0] * h + 7 * rx[nx - 2] * h
        #      + 45 * r[0] + 39 * r[nx - 2] - 84 * u[0] + 84 * u[nx - 2]) / h ** 5
        #
        # c5 = (4 * rxx[0] * h ** 2 + rxx[nx - 2] * h ** 2 - 30 * rx[0] * h
        #          + 15 * rx[nx - 2] * h + 120 * r[0] + 90 * r[nx - 2] - 210 * u[0] + 210 * u[nx - 2]) / (6 * h ** 4)
        #
        # c4 = rxx[0] / (6 * h)
        #
        # c3 = rx[0] / (2 * h)
        #
        # c2 = r[0] / h
        #
        # c1 = u[0]

        un[0] = nu**7*(-60*ri - 60*ri1 + 12*rxi - 12*rxi1 - rxxi - rxxi1 + 120*ui - 120*ui1)/6 + nu**6*(216*ri + 204*ri1 - 45*rxi + 39*rxi1 + 4*rxxi + 3*rxxi1 - 420*ui + 420*ui1)/6 + nu**5*(-90*ri - 78*ri1 + 20*rxi - 14*rxi1 - 2*rxxi - rxxi1 + 168*ui - 168*ui1)/2 + nu**4*(120*ri + 90*ri1 - 30*rxi + 15*rxi1 + 4*rxxi + rxxi1 - 210*ui + 210*ui1)/6 - nu**3*rxxi/6 + nu**2*rxi/2 - nu*ri + ui

        rn[0] = 7*nu**6*(60*ri + 60*ri1 - 12*rxi + 12*rxi1 + rxxi + rxxi1 - 120*ui + 120*ui1)/6 + nu**5*(-216*ri - 204*ri1 + 45*rxi - 39*rxi1 - 4*rxxi - 3*rxxi1 + 420*ui - 420*ui1) + 5*nu**4*(90*ri + 78*ri1 - 20*rxi + 14*rxi1 + 2*rxxi + rxxi1 - 168*ui + 168*ui1)/2 + 2*nu**3*(-120*ri - 90*ri1 + 30*rxi - 15*rxi1 - 4*rxxi - rxxi1 + 210*ui - 210*ui1)/3 + nu**2*rxxi/2 - nu*rxi + ri

        rxn[0] = 7*nu**5*(-60*ri - 60*ri1 + 12*rxi - 12*rxi1 - rxxi - rxxi1 + 120*ui - 120*ui1) + 5*nu**4*(216*ri + 204*ri1 - 45*rxi + 39*rxi1 + 4*rxxi + 3*rxxi1 - 420*ui + 420*ui1) + 10*nu**3*(-90*ri - 78*ri1 + 20*rxi - 14*rxi1 - 2*rxxi - rxxi1 + 168*ui - 168*ui1) + 2*nu**2*(120*ri + 90*ri1 - 30*rxi + 15*rxi1 + 4*rxxi + rxxi1 - 210*ui + 210*ui1) - nu*rxxi + rxi

        rxxn[0] = 35*nu**4*(60*ri + 60*ri1 - 12*rxi + 12*rxi1 + rxxi + rxxi1 - 120*ui + 120*ui1) + 20*nu**3*(-216*ri - 204*ri1 + 45*rxi - 39*rxi1 - 4*rxxi - 3*rxxi1 + 420*ui - 420*ui1) + 30*nu**2*(90*ri + 78*ri1 - 20*rxi + 14*rxi1 + 2*rxxi + rxxi1 - 168*ui + 168*ui1) - 4*nu*(120*ri + 90*ri1 - 30*rxi + 15*rxi1 + 4*rxxi + rxxi1 - 210*ui + 210*ui1) + rxxi

        for i in range(1, nx - 1):
            ui = u[i]
            ui1 = u[i - 1]
            ri = r[i]
            ri1 = r[i - 1]
            rxi = rx[i]
            rxi1 = rx[i - 1]
            rxxi = rxx[i]
            rxxi1 = rxx[i - 1]

            # c8 = (rxx[i] * h ** 2 + rxx[i - 1] * h ** 2 - 12 * rx[i] * h + 12 * rx[i - 1] * h
            #      + 60 * r[i] + 60 * r[i - 1] - 120 * u[i] + 120 * u[i - 1]) / (6 * h ** 7)
            #
            # c7 = (4 * rxx[i] * h ** 2 + 3 * rxx[i - 1] * h ** 2 - 45 * rx[i] * h + 39 * rx[i - 1] * h
            #   + 216 * r[i] + 204 * r[i - 1] - 420 * u[i] + 420 * u[i - 1]) / (6 * h ** 6)
            #
            # c6 = (rxx[i] * h ** 2 + rxx[i - 1] * h ** 2 / 2 - 10 * rx[i] * h + 7 * rx[i - 1] * h
            #  + 45 * r[i] + 39 * r[i - 1] - 84 * u[i] + 84 * u[i - 1]) / h ** 5
            #
            # c5 = (4 * rxx[i] * h ** 2 + rxx[i - 1] * h ** 2 - 30 * rx[i] * h
            #      + 15 * rx[i - 1] * h + 120 * r[i] + 90 * r[i - 1] - 210 * u[i] + 210 * u[i - 1]) / (6 * h ** 4)
            #
            # c4 = rxx[i] / (6 * h)
            #
            # c3 = rx[i] / (2 * h)
            #
            # c2 = r[i] / h
            #
            # c1 = u[i]

            un[i] = nu ** 7 * (
                    -60 * ri - 60 * ri1 + 12 * rxi - 12 * rxi1 - rxxi - rxxi1 + 120 * ui - 120 * ui1) / 6 + nu ** 6 * (
                            216 * ri + 204 * ri1 - 45 * rxi + 39 * rxi1 + 4 * rxxi + 3 * rxxi1 - 420 * ui + 420 * ui1) / 6 + nu ** 5 * (
                            -90 * ri - 78 * ri1 + 20 * rxi - 14 * rxi1 - 2 * rxxi - rxxi1 + 168 * ui - 168 * ui1) / 2 + nu ** 4 * (
                            120 * ri + 90 * ri1 - 30 * rxi + 15 * rxi1 + 4 * rxxi + rxxi1 - 210 * ui + 210 * ui1) / 6 - nu ** 3 * rxxi / 6 + nu ** 2 * rxi / 2 - nu * ri + ui

            rn[i] = 7 * nu ** 6 * (
                    60 * ri + 60 * ri1 - 12 * rxi + 12 * rxi1 + rxxi + rxxi1 - 120 * ui + 120 * ui1) / 6 + nu ** 5 * (
                            -216 * ri - 204 * ri1 + 45 * rxi - 39 * rxi1 - 4 * rxxi - 3 * rxxi1 + 420 * ui - 420 * ui1) + 5 * nu ** 4 * (
                            90 * ri + 78 * ri1 - 20 * rxi + 14 * rxi1 + 2 * rxxi + rxxi1 - 168 * ui + 168 * ui1) / 2 + 2 * nu ** 3 * (
                            -120 * ri - 90 * ri1 + 30 * rxi - 15 * rxi1 - 4 * rxxi - rxxi1 + 210 * ui - 210 * ui1) / 3 + nu ** 2 * rxxi / 2 - nu * rxi + ri

            rxn[i] = 7 * nu ** 5 * (
                    -60 * ri - 60 * ri1 + 12 * rxi - 12 * rxi1 - rxxi - rxxi1 + 120 * ui - 120 * ui1) + 5 * nu ** 4 * (
                             216 * ri + 204 * ri1 - 45 * rxi + 39 * rxi1 + 4 * rxxi + 3 * rxxi1 - 420 * ui + 420 * ui1) + 10 * nu ** 3 * (
                             -90 * ri - 78 * ri1 + 20 * rxi - 14 * rxi1 - 2 * rxxi - rxxi1 + 168 * ui - 168 * ui1) + 2 * nu ** 2 * (
                             120 * ri + 90 * ri1 - 30 * rxi + 15 * rxi1 + 4 * rxxi + rxxi1 - 210 * ui + 210 * ui1) - nu * rxxi + rxi

            rxxn[i] = 35 * nu ** 4 * (
                    60 * ri + 60 * ri1 - 12 * rxi + 12 * rxi1 + rxxi + rxxi1 - 120 * ui + 120 * ui1) + 20 * nu ** 3 * (
                              -216 * ri - 204 * ri1 + 45 * rxi - 39 * rxi1 - 4 * rxxi - 3 * rxxi1 + 420 * ui - 420 * ui1) + 30 * nu ** 2 * (
                              90 * ri + 78 * ri1 - 20 * rxi + 14 * rxi1 + 2 * rxxi + rxxi1 - 168 * ui + 168 * ui1) - 4 * nu * (
                              120 * ri + 90 * ri1 - 30 * rxi + 15 * rxi1 + 4 * rxxi + rxxi1 - 210 * ui + 210 * ui1) + rxxi

        u[:] = un[:]
        r[:] = rn[:]
        rx[:] = rxn[:]
        rxx[:] = rxxn[:]

        t = t + tau

    for i in range(0, nx):
        delta = u[i] - ua[i]
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    error[0] = h * np.sum(np.abs(np.subtract(ua, u)))
    error[1] = math.sqrt(error[1])

    return u, ua


def deltaP9sym(x, nu, tau, h, tfin, error):
    nx = len(x)
    u = np.zeros(nx)
    r = np.zeros(nx)
    un = np.zeros(nx)
    rn = np.zeros(nx)
    rxn = np.zeros(nx)
    rxxn = np.zeros(nx)
    r3xn = np.zeros(nx)

    u = test.test_func_sin(x)
    ua = test.test_func_sin(x)

    r = test.der_test_func_sin(x) * h
    rx = test.der2_test_func_sin(x) * h ** 2
    rxx = test.der3_test_func_sin(x) * h ** 3
    r3x = test.der4_test_func_sin(x) * h ** 4

    t, a = 0, 1

    while (t < tfin):
        ui = u[0]
        ui1 = u[nx - 2]
        ri = r[0]
        ri1 = r[nx - 2]
        rxi = rx[0]
        rxi1 = rx[nx - 2]
        rxxi = rxx[0]
        rxxi1 = rxx[nx - 2]
        r3xi = r3x[0]
        r3xi1 = r3x[nx - 2]

        un[0] = nu**9*(-r3xi + r3xi1 + 840*ri + 840*ri1 - 180*rxi + 180*rxi1 + 20*rxxi + 20*rxxi1 - 1680*ui + 1680*ui1)\
                /24 + nu**8*(5*r3xi - 4*r3xi1 - 3840*ri - 3720*ri1 + 840*rxi - 780*rxi1 - 96*rxxi - 84*rxxi1 + 7560*ui
                             - 7560*ui1)/24 + nu**7*(-5*r3xi + 3*r3xi1 + 3360*ri + 3120*ri1 - 756*rxi + 636*rxi1
                                                     + 90*rxxi + 66*rxxi1 - 6480*ui + 6480*ui1)/12 \
                + nu**6*(5*r3xi - 2*r3xi1 - 2688*ri - 2352*ri1 + 630*rxi - 462*rxi1 - 80*rxxi - 46*rxxi1 + 5040*ui
                         - 5040*ui1)/12 + nu**5*(-5*r3xi + r3xi1 + 1680*ri + 1344*ri1 - 420*rxi + 252*rxi1 + 60*rxxi
                                                 + 24*rxxi1 - 3024*ui + 3024*ui1)/24 + nu**4*r3xi/24 - nu**3*rxxi/6 \
                + nu**2*rxi/2 - nu*ri + ui

        rn[0] = 3*nu**8*(r3xi - r3xi1 - 840*ri - 840*ri1 + 180*rxi - 180*rxi1 - 20*rxxi - 20*rxxi1 + 1680*ui - 1680*ui1)\
                /8 + nu**7*(-5*r3xi + 4*r3xi1 + 3840*ri + 3720*ri1 - 840*rxi + 780*rxi1 + 96*rxxi + 84*rxxi1 - 7560*ui
                            + 7560*ui1)/3 + 7*nu**6*(5*r3xi - 3*r3xi1 - 3360*ri - 3120*ri1 + 756*rxi - 636*rxi1
                                                     - 90*rxxi - 66*rxxi1 + 6480*ui - 6480*ui1)/12 \
                + nu**5*(-5*r3xi + 2*r3xi1 + 2688*ri + 2352*ri1 - 630*rxi + 462*rxi1 + 80*rxxi + 46*rxxi1 - 5040*ui
                         + 5040*ui1)/2 + 5*nu**4*(5*r3xi - r3xi1 - 1680*ri - 1344*ri1 + 420*rxi - 252*rxi1 - 60*rxxi
                                                  - 24*rxxi1 + 3024*ui - 3024*ui1)/24 - nu**3*r3xi/6 + nu**2*rxxi/2 \
                - nu*rxi + ri

        rxn[0] = 3*nu**7*(-r3xi + r3xi1 + 840*ri + 840*ri1 - 180*rxi + 180*rxi1 + 20*rxxi + 20*rxxi1 - 1680*ui
                          + 1680*ui1) + 7*nu**6*(5*r3xi - 4*r3xi1 - 3840*ri - 3720*ri1 + 840*rxi - 780*rxi1 - 96*rxxi
                                                 - 84*rxxi1 + 7560*ui - 7560*ui1)/3 + 7*nu**5*(-5*r3xi + 3*r3xi1
                                                                                               + 3360*ri + 3120*ri1
                                                                                               - 756*rxi + 636*rxi1
                                                                                               + 90*rxxi + 66*rxxi1
                                                                                               - 6480*ui + 6480*ui1)/2 \
                 + 5*nu**4*(5*r3xi - 2*r3xi1 - 2688*ri - 2352*ri1 + 630*rxi - 462*rxi1 - 80*rxxi - 46*rxxi1 + 5040*ui
                            - 5040*ui1)/2 + 5*nu**3*(-5*r3xi + r3xi1 + 1680*ri + 1344*ri1 - 420*rxi + 252*rxi1 + 60*rxxi
                                                     + 24*rxxi1 - 3024*ui + 3024*ui1)/6 + nu**2*r3xi/2 - nu*rxxi + rxi


        rxxn[0] = 21*nu**6*(r3xi - r3xi1 - 840*ri - 840*ri1 + 180*rxi - 180*rxi1 - 20*rxxi - 20*rxxi1 + 1680*ui
                            - 1680*ui1) + 14*nu**5*(-5*r3xi + 4*r3xi1 + 3840*ri + 3720*ri1 - 840*rxi + 780*rxi1
                                                    + 96*rxxi + 84*rxxi1 - 7560*ui + 7560*ui1) \
                  + 35*nu**4*(5*r3xi - 3*r3xi1 - 3360*ri - 3120*ri1 + 756*rxi - 636*rxi1 - 90*rxxi - 66*rxxi1 + 6480*ui
                              - 6480*ui1)/2 + 10*nu**3*(-5*r3xi + 2*r3xi1 + 2688*ri + 2352*ri1 - 630*rxi + 462*rxi1
                                                        + 80*rxxi + 46*rxxi1 - 5040*ui + 5040*ui1) \
                  + 5*nu**2*(5*r3xi - r3xi1 - 1680*ri - 1344*ri1 + 420*rxi - 252*rxi1 - 60*rxxi - 24*rxxi1 + 3024*ui
                             - 3024*ui1)/2 - nu*r3xi + rxxi

        r3xn[0] = nu**5*(-126*r3xi + 126*r3xi1 + 105840*ri + 105840*ri1 - 22680*rxi + 22680*rxi1 + 2520*rxxi
                         + 2520*rxxi1 - 211680*ui + 211680*ui1) \
                  + nu**4*(350*r3xi - 280*r3xi1 - 268800*ri - 260400*ri1 + 58800*rxi - 54600*rxi1 - 6720*rxxi
                           - 5880*rxxi1 + 529200*ui - 529200*ui1) \
                  + nu**3*(-350*r3xi + 210*r3xi1 + 235200*ri + 218400*ri1 - 52920*rxi + 44520*rxi1 + 6300*rxxi
                           + 4620*rxxi1 - 453600*ui + 453600*ui1) \
                  + nu**2*(150*r3xi - 60*r3xi1 - 80640*ri - 70560*ri1 + 18900*rxi - 13860*rxi1 - 2400*rxxi - 1380*rxxi1
                           + 151200*ui - 151200*ui1) + 5*nu*(-5*r3xi + r3xi1 + 1680*ri + 1344*ri1 - 420*rxi + 252*rxi1
                                                             + 60*rxxi + 24*rxxi1 - 3024*ui + 3024*ui1) + r3xi

        for i in range(1, nx - 1):
            ui = u[i]
            ui1 = u[i - 1]
            ri = r[i]
            ri1 = r[i - 1]
            rxi = rx[i]
            rxi1 = rx[i - 1]
            rxxi = rxx[i]
            rxxi1 = rxx[i - 1]
            r3xi = r3x[i]
            r3xi1 = r3x[i - 1]

            un[i] = nu ** 9 * (
                    -r3xi + r3xi1 + 840 * ri + 840 * ri1 - 180 * rxi + 180 * rxi1 + 20 * rxxi + 20 * rxxi1 - 1680 * ui + 1680 * ui1) \
                    / 24 + nu ** 8 * (
                            5 * r3xi - 4 * r3xi1 - 3840 * ri - 3720 * ri1 + 840 * rxi - 780 * rxi1 - 96 * rxxi - 84 * rxxi1 + 7560 * ui
                            - 7560 * ui1) / 24 + nu ** 7 * (
                            -5 * r3xi + 3 * r3xi1 + 3360 * ri + 3120 * ri1 - 756 * rxi + 636 * rxi1
                            + 90 * rxxi + 66 * rxxi1 - 6480 * ui + 6480 * ui1) / 12 \
                    + nu ** 6 * (
                            5 * r3xi - 2 * r3xi1 - 2688 * ri - 2352 * ri1 + 630 * rxi - 462 * rxi1 - 80 * rxxi - 46 * rxxi1 + 5040 * ui
                            - 5040 * ui1) / 12 + nu ** 5 * (
                            -5 * r3xi + r3xi1 + 1680 * ri + 1344 * ri1 - 420 * rxi + 252 * rxi1 + 60 * rxxi
                            + 24 * rxxi1 - 3024 * ui + 3024 * ui1) / 24 + nu ** 4 * r3xi / 24 - nu ** 3 * rxxi / 6 \
                    + nu ** 2 * rxi / 2 - nu * ri + ui

            rn[i] = 3 * nu ** 8 * (
                    r3xi - r3xi1 - 840 * ri - 840 * ri1 + 180 * rxi - 180 * rxi1 - 20 * rxxi - 20 * rxxi1 + 1680 * ui - 1680 * ui1) \
                    / 8 + nu ** 7 * (
                            -5 * r3xi + 4 * r3xi1 + 3840 * ri + 3720 * ri1 - 840 * rxi + 780 * rxi1 + 96 * rxxi + 84 * rxxi1 - 7560 * ui
                            + 7560 * ui1) / 3 + 7 * nu ** 6 * (
                            5 * r3xi - 3 * r3xi1 - 3360 * ri - 3120 * ri1 + 756 * rxi - 636 * rxi1
                            - 90 * rxxi - 66 * rxxi1 + 6480 * ui - 6480 * ui1) / 12 \
                    + nu ** 5 * (
                            -5 * r3xi + 2 * r3xi1 + 2688 * ri + 2352 * ri1 - 630 * rxi + 462 * rxi1 + 80 * rxxi + 46 * rxxi1 - 5040 * ui
                            + 5040 * ui1) / 2 + 5 * nu ** 4 * (
                            5 * r3xi - r3xi1 - 1680 * ri - 1344 * ri1 + 420 * rxi - 252 * rxi1 - 60 * rxxi
                            - 24 * rxxi1 + 3024 * ui - 3024 * ui1) / 24 - nu ** 3 * r3xi / 6 + nu ** 2 * rxxi / 2 \
                    - nu * rxi + ri

            rxn[i] = 3 * nu ** 7 * (
                    -r3xi + r3xi1 + 840 * ri + 840 * ri1 - 180 * rxi + 180 * rxi1 + 20 * rxxi + 20 * rxxi1 - 1680 * ui
                    + 1680 * ui1) + 7 * nu ** 6 * (
                             5 * r3xi - 4 * r3xi1 - 3840 * ri - 3720 * ri1 + 840 * rxi - 780 * rxi1 - 96 * rxxi
                             - 84 * rxxi1 + 7560 * ui - 7560 * ui1) / 3 + 7 * nu ** 5 * (-5 * r3xi + 3 * r3xi1
                                                                                         + 3360 * ri + 3120 * ri1
                                                                                         - 756 * rxi + 636 * rxi1
                                                                                         + 90 * rxxi + 66 * rxxi1
                                                                                         - 6480 * ui + 6480 * ui1) / 2 \
                     + 5 * nu ** 4 * (
                             5 * r3xi - 2 * r3xi1 - 2688 * ri - 2352 * ri1 + 630 * rxi - 462 * rxi1 - 80 * rxxi - 46 * rxxi1 + 5040 * ui
                             - 5040 * ui1) / 2 + 5 * nu ** 3 * (
                             -5 * r3xi + r3xi1 + 1680 * ri + 1344 * ri1 - 420 * rxi + 252 * rxi1 + 60 * rxxi
                             + 24 * rxxi1 - 3024 * ui + 3024 * ui1) / 6 + nu ** 2 * r3xi / 2 - nu * rxxi + rxi

            rxxn[i] = 21 * nu ** 6 * (
                    r3xi - r3xi1 - 840 * ri - 840 * ri1 + 180 * rxi - 180 * rxi1 - 20 * rxxi - 20 * rxxi1 + 1680 * ui
                    - 1680 * ui1) + 14 * nu ** 5 * (
                              -5 * r3xi + 4 * r3xi1 + 3840 * ri + 3720 * ri1 - 840 * rxi + 780 * rxi1
                              + 96 * rxxi + 84 * rxxi1 - 7560 * ui + 7560 * ui1) \
                      + 35 * nu ** 4 * (
                              5 * r3xi - 3 * r3xi1 - 3360 * ri - 3120 * ri1 + 756 * rxi - 636 * rxi1 - 90 * rxxi - 66 * rxxi1 + 6480 * ui
                              - 6480 * ui1) / 2 + 10 * nu ** 3 * (
                              -5 * r3xi + 2 * r3xi1 + 2688 * ri + 2352 * ri1 - 630 * rxi + 462 * rxi1
                              + 80 * rxxi + 46 * rxxi1 - 5040 * ui + 5040 * ui1) \
                      + 5 * nu ** 2 * (
                              5 * r3xi - r3xi1 - 1680 * ri - 1344 * ri1 + 420 * rxi - 252 * rxi1 - 60 * rxxi - 24 * rxxi1 + 3024 * ui
                              - 3024 * ui1) / 2 - nu * r3xi + rxxi

            r3xn[i] = nu ** 5 * (
                    -126 * r3xi + 126 * r3xi1 + 105840 * ri + 105840 * ri1 - 22680 * rxi + 22680 * rxi1 + 2520 * rxxi
                    + 2520 * rxxi1 - 211680 * ui + 211680 * ui1) \
                      + nu ** 4 * (
                              350 * r3xi - 280 * r3xi1 - 268800 * ri - 260400 * ri1 + 58800 * rxi - 54600 * rxi1 - 6720 * rxxi
                              - 5880 * rxxi1 + 529200 * ui - 529200 * ui1) \
                      + nu ** 3 * (
                              -350 * r3xi + 210 * r3xi1 + 235200 * ri + 218400 * ri1 - 52920 * rxi + 44520 * rxi1 + 6300 * rxxi
                              + 4620 * rxxi1 - 453600 * ui + 453600 * ui1) \
                      + nu ** 2 * (
                              150 * r3xi - 60 * r3xi1 - 80640 * ri - 70560 * ri1 + 18900 * rxi - 13860 * rxi1 - 2400 * rxxi - 1380 * rxxi1
                              + 151200 * ui - 151200 * ui1) + 5 * nu * (
                              -5 * r3xi + r3xi1 + 1680 * ri + 1344 * ri1 - 420 * rxi + 252 * rxi1
                              + 60 * rxxi + 24 * rxxi1 - 3024 * ui + 3024 * ui1) + r3xi

        u[:] = un[:]
        r[:] = rn[:]
        rx[:] = rxn[:]
        rxx[:] = rxxn[:]
        r3x[:] = r3xn[:]

        t = t + tau

    for i in range(0, nx):
        delta = u[i] - ua[i]
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    error[0] = h * np.sum(np.abs(np.subtract(ua, u)))
    error[1] = math.sqrt(error[1])

    return u, ua



def two_dim_deltaP3(err, nx, ny):
    x = np.linspace(0, 1, nx + 1)
    y = np.linspace(0, 1, ny + 1)
    # nx = 50
    # ny = 50
    r = 0.1
    tfin = 1
    h = 1 / nx
    tau = r * h
    count = 0
    p = np.zeros((nx, ny))
    pa = np.zeros((nx, ny))

    u = np.zeros((nx, ny))
    ua = np.zeros((nx, ny))

    v = np.zeros((nx, ny))
    va = np.zeros((nx, ny))

    rp = np.zeros((nx, ny))
    ru = np.zeros((nx, ny))
    rv = np.zeros((nx, ny))

    sp = np.zeros((nx, ny))
    su = np.zeros((nx, ny))
    sv = np.zeros((nx, ny))

    pn = np.zeros((nx, ny))
    un = np.zeros((nx, ny))
    vn = np.zeros((nx, ny))

    rpn = np.zeros((nx, ny))
    run = np.zeros((nx, ny))
    rvn = np.zeros((nx, ny))

    spn = np.zeros((nx, ny))
    sun = np.zeros((nx, ny))
    svn = np.zeros((nx, ny))

    for i in range(0, nx):
        for j in range(0, ny):
            p[i][j] = test.two_dim_test_sin(x[i], y[j])
            rp[i][j] = test.two_dim_test_sin_r(x[i], y[j])
            sp[i][j] = test.two_dim_test_sin_s(x[i], y[j])
            pa[i][j] = test.two_dim_test_sin(x[i], y[j])

            pn[i][j] = test.two_dim_test_sin(x[i], y[j])
            rpn[i][j] = test.two_dim_test_sin_r(x[i], y[j])
            spn[i][j] = test.two_dim_test_sin_s(x[i], y[j])
    #     pn[i][0] = test.two_dim_test_exp(x[i], y[0])
    #     rpn[i][0] = test.two_dim_test_exp(x[i], y[0])
    #     spn[i][0] = test.two_dim_test_exp(x[i], y[0])
    #     pn[i][ny - 1] = test.two_dim_test_exp(x[i], y[ny - 1])
    #     rpn[i][ny - 1] = test.two_dim_test_exp(x[i], y[ny - 1])
    #     spn[i][ny - 1] = test.two_dim_test_exp(x[i], y[ny - 1])
    # #
    # for j in range(0, ny):
    #     pn[0][j] = test.two_dim_test_exp(x[0], y[j])
    #     rpn[0][j] = test.two_dim_test_exp(x[0], y[j])
    #     spn[0][j] = test.two_dim_test_exp(x[0], y[j])
    #     pn[nx - 1][j] = test.two_dim_test_exp(x[nx - 1], y[j])
    #     rpn[nx - 1][j] = test.two_dim_test_exp(x[nx - 1], y[j])
    #     spn[nx - 1][j] = test.two_dim_test_exp(x[nx - 1], y[j])

        # p[i][j] = test.two_dim_test_func_p(x[i], y[j], 0)
            # u[i][j] = test.two_dim_test_func_u(x[i], y[j], 0)
            # v[i][j] = test.two_dim_test_func_v(x[i], y[j], 0)

            # rp[i][j] = test.two_dim_test_der_rp(x[i], y[j], 0) * h
            # ru[i][j] = test.two_dim_test_der_ru(x[i], y[j], 0) * h
            # rv[i][j] = test.two_dim_test_der_rv(x[i], y[j], 0) * h
            #
            # sp[i][j] = test.two_dim_test_der_sp(x[i], y[j], 0) * h
            # su[i][j] = test.two_dim_test_der_su(x[i], y[j], 0) * h
            # sv[i][j] = test.two_dim_test_der_sv(x[i], y[j], 0) * h
    t = 0
    # print("p = {0}".format(p[0][:]))

    dx, dy = 0.005, 0.005
    # fig, ax = plt.subplots(nrows=1)
    # levels0 = MaxNLocator(nbins=15).tick_values(p.min(), p.max())
    # cmap = plt.get_cmap('PiYG')
    # #
    # ax.set_xlabel(r'\textbf{x}')
    # ax.set_ylabel(r'\textbf{y}')
    while (t < tfin):
        print(t)
        phi_i_j = p[0][0]
        phi_i1_j = p[nx - 2][0]
        phi_i_j1 = p[0][ny - 2]
        phi_i1_j1 = p[nx - 2][ny - 2]
        # #
        r_i1_j = rp[nx - 2][0]
        r_i_j1 = rp[0][ny - 2]
        r_i_j = rp[0][0]
        # #
        s_i1_j = sp[nx - 2][0]
        s_i_j1 = sp[0][ny - 2]
        s_i_j = sp[0][0]
        # #
        pn[0][0] = (3 * h ** 3 * (phi_i_j - tau * (r_i_j + s_i_j))
                    + 3 * h * tau ** 2 * (h * r_i1_j + 2 * h * r_i_j + 2 * h * s_i_j + h * s_i_j1 + 3 * phi_i1_j
                                          - 6 * phi_i_j + 3 * phi_i_j1) + tau ** 3 * (
                            -3 * h * r_i1_j - 4 * h * r_i_j
                            + h * r_i_j1 + h * s_i1_j - 4 * h * s_i_j
                            - 3 * h * s_i_j1 - 8 * phi_i1_j
                            + 2 * phi_i1_j1 + 14 * phi_i_j
                            - 8 * phi_i_j1)) / (3 * h ** 3)
        # #
        rpn[0][0] = (h ** 3 * r_i_j - h * tau * (
                2 * h * r_i1_j + 5 * h * r_i_j - h * r_i_j1 - h * s_i1_j + h * s_i_j
                + 7 * phi_i1_j - phi_i1_j1 - 7 * phi_i_j + phi_i_j1)
                     + tau ** 2 * (
                           3 * h * r_i1_j + 3 * h * r_i_j - h * s_i1_j + h * s_i_j + 7 * phi_i1_j - phi_i1_j1
                             - 7 * phi_i_j + phi_i_j1)) / h ** 3
        #
        spn[0][0] = (h ** 3 * s_i_j - h * tau * (
                h * r_i_j - h * r_i_j1 - h * s_i1_j + 5 * h * s_i_j + 2 * h * s_i_j1
                + phi_i1_j - phi_i1_j1 - 7 * phi_i_j + 7 * phi_i_j1)
                     + tau ** 2 * (h * r_i_j - h * r_i_j1 + 3 * h * s_i_j + 3 * h * s_i_j1 + phi_i1_j - phi_i1_j1
                                   - 7 * phi_i_j + 7 * phi_i_j1)) / h ** 3

        for i in range(1, nx - 1):
            phi_i_j = p[i][0]
            phi_i1_j = p[i - 1][0]
            phi_i_j1 = p[i][ny - 2]
            phi_i1_j1 = p[i - 1][ny - 2]
            # #
            r_i1_j = rp[i - 1][0]
            r_i_j1 = rp[i][ny - 2]
            r_i_j = rp[i][0]
            # #
            s_i1_j = sp[i - 1][0]
            s_i_j1 = sp[i][ny - 2]
            s_i_j = sp[i][0]
            # #
            pn[i][0] = (3 * h ** 3 * (phi_i_j - tau * (r_i_j + s_i_j))
                        + 3 * h * tau ** 2 * (h * r_i1_j + 2 * h * r_i_j + 2 * h * s_i_j + h * s_i_j1 + 3 * phi_i1_j
                                              - 6 * phi_i_j + 3 * phi_i_j1) + tau ** 3 * (
                                -3 * h * r_i1_j - 4 * h * r_i_j
                                + h * r_i_j1 + h * s_i1_j - 4 * h * s_i_j
                                - 3 * h * s_i_j1 - 8 * phi_i1_j
                                + 2 * phi_i1_j1 + 14 * phi_i_j
                                - 8 * phi_i_j1)) / (3 * h ** 3)
            # #
            rpn[i][0] = (h ** 3 * r_i_j - h * tau * (
                    2 * h * r_i1_j + 5 * h * r_i_j - h * r_i_j1 - h * s_i1_j + h * s_i_j
                    + 7 * phi_i1_j - phi_i1_j1 - 7 * phi_i_j + phi_i_j1)
                         + tau ** 2 * (
                                 3 * h * r_i1_j + 3 * h * r_i_j - h * s_i1_j + h * s_i_j + 7 * phi_i1_j - phi_i1_j1
                                 - 7 * phi_i_j + phi_i_j1)) / h ** 3

            spn[i][0] = (h ** 3 * s_i_j - h * tau * (
                    h * r_i_j - h * r_i_j1 - h * s_i1_j + 5 * h * s_i_j + 2 * h * s_i_j1
                    + phi_i1_j - phi_i1_j1 - 7 * phi_i_j + 7 * phi_i_j1)
                         + tau ** 2 * (
                                 h * r_i_j - h * r_i_j1 + 3 * h * s_i_j + 3 * h * s_i_j1 + phi_i1_j - phi_i1_j1
                                 - 7 * phi_i_j + 7 * phi_i_j1)) / h ** 3

            for j in range(1, ny - 1):
                phi_i_j = p[0][j]
                phi_i1_j = p[nx - 2][j]
                phi_i_j1 = p[0][j - 1]
                phi_i1_j1 = p[nx - 2][j - 1]
                # #
                r_i1_j = rp[nx - 2][j]
                r_i_j1 = rp[0][j - 1]
                r_i_j = rp[0][j]
                # #
                s_i1_j = sp[nx - 2][j]
                s_i_j1 = sp[0][j - 1]
                s_i_j = sp[0][j]
                # #
                pn[0][j] = (3 * h ** 3 * (phi_i_j - tau * (r_i_j + s_i_j))
                            + 3 * h * tau ** 2 * (h * r_i1_j + 2 * h * r_i_j + 2 * h * s_i_j + h * s_i_j1 + 3 * phi_i1_j
                                                  - 6 * phi_i_j + 3 * phi_i_j1) + tau ** 3 * (
                                    -3 * h * r_i1_j - 4 * h * r_i_j
                                    + h * r_i_j1 + h * s_i1_j - 4 * h * s_i_j
                                    - 3 * h * s_i_j1 - 8 * phi_i1_j
                                    + 2 * phi_i1_j1 + 14 * phi_i_j
                                    - 8 * phi_i_j1)) / (3 * h ** 3)
                # #
                rpn[0][j] = (h ** 3 * r_i_j - h * tau * (
                        2 * h * r_i1_j + 5 * h * r_i_j - h * r_i_j1 - h * s_i1_j + h * s_i_j
                        + 7 * phi_i1_j - phi_i1_j1 - 7 * phi_i_j + phi_i_j1)
                             + tau ** 2 * (
                                     3 * h * r_i1_j + 3 * h * r_i_j - h * s_i1_j + h * s_i_j + 7 * phi_i1_j - phi_i1_j1
                                     - 7 * phi_i_j + phi_i_j1)) / h ** 3
                # #
                spn[0][j] = (h ** 3 * s_i_j - h * tau * (
                        h * r_i_j - h * r_i_j1 - h * s_i1_j + 5 * h * s_i_j + 2 * h * s_i_j1
                        + phi_i1_j - phi_i1_j1 - 7 * phi_i_j + 7 * phi_i_j1)
                             + tau ** 2 * (
                                     h * r_i_j - h * r_i_j1 + 3 * h * s_i_j + 3 * h * s_i_j1 + phi_i1_j - phi_i1_j1
                                     - 7 * phi_i_j + 7 * phi_i_j1)) / h ** 3

            for j in range(1, ny - 1):
                phi_i_j = p[i][j]
                phi_i1_j = p[i - 1][j]
                phi_i_j1 = p[i][j - 1]
                phi_i1_j1 = p[i - 1][j - 1]

                r_i1_j = rp[i - 1][j]
                r_i_j1 = rp[i][j - 1]
                r_i_j = rp[i][j]

                s_i1_j = sp[i - 1][j]
                s_i_j1 = sp[i][j - 1]
                s_i_j = sp[i][j]

                pn[i][j] = (3 * h ** 3 * (phi_i_j - tau * (r_i_j + s_i_j))
                            + 3 * h * tau ** 2 * (
                                    h * r_i1_j + 2 * h * r_i_j + 2 * h * s_i_j + h * s_i_j1 + 3 * phi_i1_j
                                    - 6 * phi_i_j + 3 * phi_i_j1) + tau ** 3 * (-3 * h * r_i1_j - 4 * h * r_i_j
                                                                                + h * r_i_j1 + h * s_i1_j - 4 * h * s_i_j
                                                                                - 3 * h * s_i_j1 - 8 * phi_i1_j
                                                                                + 2 * phi_i1_j1 + 14 * phi_i_j
                                                                                - 8 * phi_i_j1)) / (3 * h ** 3)
                rpn[i][j] = (h ** 3 * r_i_j - h * tau * (
                        2 * h * r_i1_j + 5 * h * r_i_j - h * r_i_j1 - h * s_i1_j + h * s_i_j
                        + 7 * phi_i1_j - phi_i1_j1 - 7 * phi_i_j + phi_i_j1)
                             + tau ** 2 * (
                                     3 * h * r_i1_j + 3 * h * r_i_j - h * s_i1_j + h * s_i_j + 7 * phi_i1_j - phi_i1_j1
                                     - 7 * phi_i_j + phi_i_j1)) / h ** 3

                spn[i][j] = (h ** 3 * s_i_j - h * tau * (
                        h * r_i_j - h * r_i_j1 - h * s_i1_j + 5 * h * s_i_j + 2 * h * s_i_j1
                        + phi_i1_j - phi_i1_j1 - 7 * phi_i_j + 7 * phi_i_j1)
                             + tau ** 2 * (
                                     h * r_i_j - h * r_i_j1 + 3 * h * s_i_j + 3 * h * s_i_j1 + phi_i1_j - phi_i1_j1
                                     - 7 * phi_i_j + 7 * phi_i_j1)) / h ** 3

        # print("pn = {0}".format(pn[i][:]))
        # print("rpn = {0}".format(rpn[i][:]))
        # print("spn = {0}\n".format(spn[i][:]))

                    # if(n == 1):
                    #     un[i][j] = c1 + c10*h**3 + c2*h + c3*h + c4*h**2 + c5*h**2 + c6*h**2 + c7*h**3 + c8*h**3 + c9*h**3 - h*tau**3*(c8 + 3*c9)/3 + h*tau**2*(c4 + 2*c5 + 4*c7*h + 2*c8*h + 6*c9*h)/2 - h*tau*(c2 + c4*h + 2*c5*h + 2*c7*h**2 + c8*h**2 + 3*c9*h**2)
                    #
                    #     run[i][j] = h*(c2 + c4*h + 2*c5*h + 2*c7*h**2 + c8*h**2 + 3*c9*h**2 + tau**2*(c7 + 3*c9) - 2*tau*(c5 + c7*h + 3*c9*h))
                    #
                    #     sun[i][j] = h*(3*c10*h**2 + c3 + c4*h + 2*c6*h + c7*h**2 + 2*c8*h**2 + tau**2*(c7 + c8) - tau*(c4 + 2*c7*h + 2*c8*h))
                    #
                    #
                    # if(n == 2):
                    #
                    #     vn[i][j] = c1 + c10*h**3 + c2*h + c3*h + c4*h**2 + c5*h**2 + c6*h**2 + c7*h**3 + c8*h**3 + c9*h**3 - h*tau**3*(3*c10 + c7)/3 + h*tau**2*(6*c10*h + c4 + 2*c6 + 2*c7*h + 4*c8*h)/2 - h*tau*(3*c10*h**2 + c3 + c4*h + 2*c6*h + c7*h**2 + 2*c8*h**2)
                    #
                    #     rvn[i][j] = h*(c2 + c4*h + 2*c5*h + 2*c7*h**2 + c8*h**2 + 3*c9*h**2 + tau**2*(c7 + c8) - tau*(c4 + 2*c7*h + 2*c8*h))
                    #
                    #     svn[i][j] = h*(3*c10*h**2 + c3 + c4*h + 2*c6*h + c7*h**2 + 2*c8*h**2 + tau**2*(3*c10 + c8) - 2*tau*(3*c10*h + c6 + c8*h))
        for i in range(0, nx):
            for j in range(0, ny):
                p[i][j] = pn[i][j]
                rp[i][j] = rpn[i][j]
                sp[i][j] = spn[i][j]

        # p[:][:] = pn[:][:]
        u[:][:] = un[:][:]
        v[:][:] = vn[:][:]
        # print("p = {0}\n".format(p))

        # rp[:][:] = rpn[:][:]
        ru[:][:] = run[:][:]
        rv[:][:] = rvn[:][:]

        # sp[:][:] = spn[:][:]
        su[:][:] = sun[:][:]
        sv[:][:] = svn[:][:]
        # ax.clear()

        # CS = ax.contourf(x[:-1] + dx / 2., y[:-1] + dy / 2., p, levels=levels0)
        # plt.contour(CS, colors='k')
        # plt.grid(c='k', ls='-', alpha=0.3)
        #
        # # fig.colorbar(cf0, ax=ax)
        # plt.savefig("2D_pics/pic{0}.png".format(count), dpi=250)
        count = count + 1
        t = t + tau

    # print("pa = {0}".format(pa))

    for i in range(0, ny):
        delta = p[1][i] - pa[1][i]
        err[1] += delta ** 2
        if (math.fabs(delta) > err[2]):
            err[2] = math.fabs(delta)

    err[0] = h * np.sum(np.abs(np.subtract(pa[1][:], p[1][:])))
    err[1] = math.sqrt(err[1])
    print("error = {0}".format(err[0]))

    # plt.show()

    return p, pa



def convergence(error):
    result = np.zeros(error.size)
    for i in range(1, error.size):
        result[i] = (np.log(error[i - 1] / error[i])) / (np.log(2))
    print(result)
    return result