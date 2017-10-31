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
        error[0] += delta
        error[1] += delta ** 2
        if (delta > error[2]):
            error[2] = delta
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


def general_formula(u, ux, nt, nx, h, ksi):
    error = np.zeros(3)
    h2 = h ** 2
    xi = ksi
    xi2 = ksi ** 2
    xi3 = ksi ** 3
    for n in range(0, nt - 1):
        for i in range(0, nx - 1):
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

############# Задание №1 todo_2017 ##############


def cip3(u, ux, nt, nx, h, xi):
    error = np.zeros(3)
    for n in range(0, nt - 1):
        for i in range(1, nx - 1):
            # BC #
            a = (ux[n][0] + ux[n][nx - 1]) / h**2 - 2 * (u[n][0] - u[n][nx - 1]) / h**3

            b = (2 * ux[n][0] + ux[n][nx - 1]) / h - 3 * (u[n][0] - u[n][nx - 1]) / h**2

            c = ux[n][0]

            d = u[n][0]

            u[n + 1][0] = d + c * xi + b * xi**2 + a * xi**3

            ux[n + 1][0] = c + 2 * b * xi + 3 * a * xi**2

            # Main Loop #

            a = (ux[n][i] + ux[n][i - 1]) / h**2 - 2 * (u[n][i] - u[n][i - 1]) / h**3

            b = (2 * ux[n][i] + ux[n][i - 1]) / h - 3 * (u[n][i] - u[n][i - 1]) / h**2

            u[n + 1][i] = u[n][i] + ux[n][i] * xi + b * xi**2 + a * xi**3

            ux[n + 1][i] = ux[n][i] + 2 * b * xi + 3 * a * xi**2

    for i in range(0, nx):
        delta = (u[nt - 1][i] - u[0][i]) * h
        error[0] += math.fabs(delta)
        error[1] += delta ** 2
        if(math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)

    error[1] = math.sqrt(error[1])
    return error[0]


def cip5(u, ux, uxx, nt, nx, h, xi):
    error = np.zeros(3)
    for n in range(0, nt - 1):
        for i in range(1, nx):
            # BC #

            a = (uxx[n][0] - uxx[n][nx - 1]) / (2 * h**3) - 6 * (ux[n][0] + ux[n][nx - 1]) / (2 * h**4) + 12 * (u[n][0] - u[n][nx - 1]) / (2 * h**5)

            b = ((3/2) * uxx[n][0] - uxx[n][nx - 1]) / h**2 - (8 * ux[n][0] + 7 * ux[n][nx - 1]) / h**3 + 15 * (u[n][0] - u[n][nx - 1]) / h**4

            c = (3 * uxx[n][0] - uxx[n][nx - 1]) / (2 * h) - (12 * ux[n][0] + 8 * ux[n][nx - 1]) / (2 * h**2) + 20 * (u[n][0] - u[n][nx - 1]) / (2 * h**3)

            d = uxx[n][0] / 2

            e = ux[n][0]

            f = u[n][0]
            # print("a = {0:3f}, b = {1:3f}, c = {2:3f}, d = {3:3f}, e = {4:3f}, f = {5:3f}".format(a,b,c,d,e,f))

            u[n + 1][0] = f + e * xi + d * xi**2 + c * xi**3 + b * xi**4 + a * xi**5

            ux[n + 1][0] = e + 2 * d * xi + 3 * c * xi**2 + 4 * b * xi**3 + 5 * a * xi**4

            uxx[n + 1][0] = 2 * d + 6 * c * xi + 12 * b * xi**2 + 20 * a * xi**3



            # Main #
            a = (uxx[n][i] - uxx[n][i - 1]) / (2 * h**3) - 6 * (ux[n][i] + ux[n][i - 1]) / (2 * h**4) + 12 * (u[n][i] - u[n][i - 1]) / (2 * h**5)

            b = ((3/2) * uxx[n][i] - uxx[n][i - 1]) / h**2 - (8 * ux[n][i] + 7 * ux[n][i - 1]) / h**3 + 15 * (u[n][i] - u[n][i - 1]) / h**4

            c = (3 * uxx[n][i] - uxx[n][i - 1]) / (2 * h) - (12 * ux[n][i] + 8 * ux[n][i - 1]) / (2 * h**2) + 20 * (u[n][i] - u[n][i - 1]) / (2 * h**3)

            d = uxx[n][i] / 2

            e = ux[n][i]

            f = u[n][i]
            # print("a = {0:3f}, b = {1:3f}, c = {2:3f}, d = {3:3f}, e = {4:3f}, f = {5:3f}".format(a,b,c,d,e,f))

            u[n + 1][i] = f + e * xi + d * xi**2 + c * xi**3 + b * xi**4 + a * xi**5

            ux[n + 1][i] = e + 2 * d * xi + 3 * c * xi**2 + 4 * b * xi**3 + 5 * a * xi**4

            uxx[n + 1][i] = 2 * d + 6 * c * xi + 12 * b * xi**2 + 20 * a * xi**3

            # print("u[n + 1][i] = {0:3f}, ux[n + 1][i] = {1:3f}, uxx[n + 1][i] = {2:3f}".format(u[n + 1][i],ux[n + 1][i],uxx[n + 1][i]))

    for i in range(0, nx - 1):
        delta = (u[nt - 1][i] - u[0][i]) * h
        error[0] += math.fabs(delta)
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    error[1] = math.sqrt(error[1])
    return error

def cip7(u, ux, u2x, u3x, nt, nx, h, xi):
    error = np.zeros(3)
    for n in range(0, nt - 1):
        for i in range(1, nx):
            # BC #
            a = (u3x[n][0] + u3x[n][nx - 1]) / (6 * h**4) - 2 * (u2x[n][0] - u2x[n][nx - 1]) / h**5 + 10 * (ux[n][0] + ux[n][nx - 1]) / h**6 - 20 * (u[n][0] - u[n][nx - 1]) / h**7

            b = (4 * u3x[n][0] + 3 * u3x[n][nx - 1]) / (6 * h**3) - (45 * u2x[n][0] - 39 * u2x[n][nx - 1]) / (6 * h**4) + (216 * ux[n][0] + 204 * ux[n][nx - 1]) / (6 * h**5) - 420 * (u[n][0] - u[n][nx - 1]) / (6 * h**6)

            c = (u3x[n][0] + 0.5 * u3x[n][nx - 1]) / h**2 - (10 * u2x[n][0] - 7 * u2x[n][nx - 1]) / h**3 + (45 * ux[n][0] + 39 * ux[n][nx - 1]) / h**4 - 84 * (u[n][0] - u[n][nx - 1]) / h**5

            d = (4 * u3x[n][0] + u3x[n][nx - 1]) / (6 * h) - 15 * (2 * u2x[n][0] - u2x[n][nx - 1]) / (6 * h**2) + 30 * (4 * ux[n][0] + 3 * ux[n][nx - 1]) / (6 * h**3) - 210 * (u[n][0] - u[n][nx - 1]) / (6 * h**4)

            e = u3x[n][0] / 6

            f = u2x[n][0] / 2

            g = ux[n][0]

            l = u[n][0]

            u[n + 1][0] = l + g * xi + f * xi ** 2 + e * xi ** 3 + d * xi ** 4 + c * xi ** 5 + b * xi ** 6 + a * xi ** 7

            ux[n + 1][0] = g + 2 * f * xi + 3 * e * xi ** 2 + 4 * d * xi ** 3 + 5 * c * xi ** 4 + 6 * b * xi ** 5 + 7 * a * xi ** 6

            u2x[n + 1][0] = 2 * f + 6 * e * xi + 12 * d * xi ** 2 + 20 * c * xi ** 3 + 30 * b * xi ** 4 + 42 * a * xi ** 5

            u3x[n + 1][0] = 6 * e + 24 * d * xi + 60 * c * xi ** 2 + 120 * b * xi ** 3 + 210 * a * xi ** 4

            # Main #
            a = (u3x[n][i] + u3x[n][i - 1]) / (6 * h**4) - 2 * (u2x[n][i] - u2x[n][i - 1]) / h**5 + 10 * (ux[n][i] + ux[n][i - 1]) / h**6 - 20 * (u[n][i] - u[n][i - 1]) / h**7

            b = (4 * u3x[n][i] + 3 * u3x[n][i - 1]) / (6 * h**3) - (45 * u2x[n][i] - 39 * u2x[n][i - 1]) / (6 * h**4) + (216 * ux[n][i] + 204 * ux[n][i - 1]) / (6 * h**5) - 420 * (u[n][i] - u[n][i - 1]) / (6 * h**6)

            c = (u3x[n][i] + 0.5 * u3x[n][i - 1]) / h**2 - (10 * u2x[n][i] - 7 * u2x[n][i - 1]) / h**3 + (45 * ux[n][i] + 39 * ux[n][i - 1]) / h**4 - 84 * (u[n][i] - u[n][i - 1]) / h**5

            d = (4 * u3x[n][i] + u3x[n][i - 1]) / (6 * h) - 15 * (2 * u2x[n][i] - u2x[n][i - 1]) / (6 * h**2) + 30 * (4 * ux[n][i] + 3 * ux[n][i - 1]) / (6 * h**3) - 210 * (u[n][i] - u[n][i - 1]) / (6 * h**4)

            e = u3x[n][i] / 6

            f = u2x[n][i] / 2

            g = ux[n][i]

            l = u[n][i]

            u[n + 1][i] = l + g * xi + f * xi**2 + e * xi**3 + d * xi**4 + c * xi**5 + b * xi**6 + a * xi**7

            ux[n + 1][i] = g + 2 * f * xi + 3 * e * xi**2 + 4 * d * xi**3 + 5 * c * xi**4 + 6 * b * xi**5 + 7 * a * xi**6

            u2x[n + 1][i] = 2 * f + 6 * e * xi + 12 * d * xi**2 + 20 * c * xi**3 + 30 * b * xi**4 + 42 * a * xi**5

            u3x[n + 1][i] = 6 * e + 24 * d * xi + 60 * c * xi**2 + 120 * b * xi**3 + 210 * a * xi**4

    for i in range(0, nx - 1):
        delta = (u[nt - 1][i] - u[0][i]) * h
        error[0] += math.fabs(delta)
        error[1] += delta ** 2
        if (math.fabs(delta) > error[2]):
            error[2] = math.fabs(delta)
    error[1] = math.sqrt(error[1])
    return error

    return error