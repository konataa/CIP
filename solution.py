import one_dimentional_solver as ods
import test_functions as test
import numpy as np

a = 1
l = 1
T = 20

xL = -1
xR = 1

# CFL
r = 0.4

# output files
out_file = open('results.dat', 'w')


# functions calling

grid_steps = np.array([8, 16, 32, 64, 128, 256])
grid_steps = np.array([100, 200, 400, 800, 1600])

# grid_steps = np.array([8])

# nx = 200
for nx in grid_steps:
    h = (xR - xL) / nx
    tau = r * h / a
    xi = -a * tau
    nt = 10
    x = np.linspace(xL, xR, nx)
    x_dence = np.linspace(xL, xR, nx)
    ux = np.zeros((nt, nx))
    u2x = np.zeros((nt, nx))
    u3x = np.zeros((nt, nx))
    u = np.zeros((nt, nx))

    for i in range(0, nx):
        u[0][i] = test.test_func(x[i])

    # for i in range(0, nx):
    #     u_dence[i] = test.test_func(x_dence[i])

    # for i in range(0, nx - 1):
    #     rho[0][i] = 0.5 * (u[0][i] + u[0][i + 1]) * h

    # l = ods.upwind_scheme(u, nt, nx, r, h)
    #
    # print("upwind")
    # print("N = ", nx)
    # print(l)
    #
    # l = ods.cubic_polynominal_interpolation(u, ux, nt, nx, h, xi)
    #
    # print("CIP")
    # print("N = ", nx)
    # print(l)
    #
    # l = ods.fourth_order_polynominal_interpolation(u, ux, nt, nx, h, xi, del_rho, rho)
    #
    # print("CIP CSL4")
    # print("N = ", nx)
    # print(l)
    #
    # l = ods.second_order_polynominal_interpolation(u, ux, nt, nx, h, xi, del_rho, rho)
    #
    # print("CIP CSL2")
    # print("N = ", nx)
    # print(l)

    l = ods.cip3(u, ux, nt, nx, h, xi)
    print("cip3")
    print("N = ", nx)
    print(l)
    # print(u)

    #
    # l = ods.cip5(u, ux, u2x, nt, nx, h, xi)
    # print("cip5")
    # print("N = ", nx)
    # print(l)
    #
    # l = ods.cip7(u, ux, u2x, u3x, nt, nx, h, xi)
    # print("cip7")
    # print("N = ", nx)
    # print(l)

out_file.close()

