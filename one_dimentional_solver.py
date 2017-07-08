import numpy as np
import math


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


def upwind_scheme(u, nt, nx, r, x):
    error = np.zeros(3)
    print(error)
    for n in range(0, nt - 1):
        for i in range(1, nx):
            u[n + 1][i] = u[n][i] - r * (u[n][i] - u[n][i - 1])
        # BC
        u[n + 1][0] = u[n][nx - 1]
    for i in range(0, nx - 1):
        delta = u[nt - 1][i] - u[0][i]
        #print("delta = ",delta)
        error[0] += math.fabs(delta)
        #print(delta)
        #print(error[0])
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    error[1] = math.sqrt(error[1])
   # print(error)
    return error


def cubic_polynominal_interpolation(u, ux, nt, nx, h, ksi):
    error = np.zeros(3)
    for n in range(0, nt - 1):
        for i in range(1, nx):
            a = (ux[n][i] + ux[n][i - 1]) / h ** 2 - 2 * (u[n][i] - u[n][i - 1]) / h ** 3

            b = (2 * ux[n][i] + ux[n][i - 1]) / h - 3 * (u[n][i] - u[n][i - 1]) / h ** 2

            u[n + 1][i] = a * ksi ** 3 + b * ksi ** 2 + ux[n][i] * ksi + u[n][i]

            ux[n + 1][i] = 3 * a * ksi ** 2 + 2 * b * ksi + ux[n][i]

            #s[n + 1][i - 1] = (ux[n + 1][i] - ux[n + 1][i - 1]) / h

        #for i in range(1, nx - 1):
        #    s_i = minmod(s[n + 1][i - 1], s[n + 1][i])
        #    ux[n + 1][i] = minmod(ux[n + 1][i - 1], 3 * s_i)

        # BC
        u[n + 1][0] = u[n][nx - 1]
        ux[n + 1][0] = ux[n][nx - 1]

    for i in range(0, nx - 1):
        delta = u[nt - 1][i] - u[0][i]
        error[0] += math.fabs(delta)
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    #error[1] = math.sqrt(error[1])
    print(error)
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
        delta = u[nt - 1][i] - u[0][i]
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
        delta = u[nt - 1][i] - u[0][i]
        error[0] += math.fabs(delta)
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    error[1] = math.sqrt(error[1])
    #print(error)
    return error


def general_formula(u, ux, nt, nx, h, ksi):
    error = np.zeros(3)
    h2 = h ** 2
    xi = ksi
    xi2 = ksi ** 2
    xi3 = ksi ** 3
    for n in range(0, nt - 1):
       # print("======= n = ",n," =========")
        for i in range(0, nx - 1):
        #    print("========= x = ",x[i], "===========")
        #    print("h = ",h)
        #    print("fi+1 = ", u[n][i + 1])
        #    print("fi = ", u[n][i])
        #    print("di+1 = ", ux[n][i + 1])
        #    print("di = ", ux[n][i])
            S_i = (u[n][i + 1] - u[n][i]) / h
        #    print("S_i = ", S_i)
            if (ux[n][i] * ux[n][i + 1] < 0):
        #        print("<0")
                B_i = (((S_i-ux[n][i])/(ux[n][i+1]-S_i))-1)/h
            else:
        #        print(">=0")
                B_i = 0.0
        #    print("B_i = ", B_i)
            A1 = ux[n][i] + u[n][i] * B_i
        #    print("A1 = ",A1)
            A3 = (ux[n][i] - S_i + (ux[n][i + 1] - S_i) * (1 + B_i * h))/ h2
        #    print("A3 = ", A3)
            A2 = S_i * B_i + ((S_i - ux[n][i]) / h) - A3 * h
        #    print("A2 = ", A2)

            u[n + 1][i] = (u[n][i] + A1 * xi + A2 * xi2 + A3 * xi3) / (1 + B_i * xi)
        #    print((A1 + 2 * A2 * xi + 3 * A3 * xi2) / (1 + B_i * xi))
            ux[n + 1][i] = (A1 + 2 * A2 * xi + 3 * A3* xi2) / (1 + B_i * xi) - ((u[n][i] + A1 * xi + A2 * xi2 + A3 * xi3) * B_i / (1 + B_i * xi) ** 2)
        #    print("d_i = ", ux[n+1][i])
        #    print("f_i = ", u[n+1][i])
        # BC
        u[n + 1][0] = u[n][nx - 1]
        ux[n + 1][0]    = ux[n][nx - 1]
    for i in range(0, nx - 1):
        delta = u[nt - 1][i] - u[0][i]
        error[0] += math.fabs(delta)
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    error[1] = math.sqrt(error[1])
    return error


#print(u)


