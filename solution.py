import one_dimentional_solver as ods
import test_functions as test
from matplotlib import pyplot as plt
import csv

import numpy as np

a = 1
l = 1
T = 1

xL = -1
xR = 1

# CFL
r = 0.5

# output files
out_file = open('results.dat', 'w')
outfile = open('convergence.csv', 'w')
wrtr = csv.writer(outfile)

vtest = np.vectorize(test.test_func)


# plot

fig = plt.figure()
fig2 = plt.figure()
fig3 = plt.figure()
fig4 = plt.figure()
x1 = fig.add_subplot(312)
x2 = fig.add_subplot(312)
x3 = fig2.add_subplot(312)
x4 = fig2.add_subplot(312)
x5 = fig3.add_subplot(312)
x6 = fig3.add_subplot(312)
x7 = fig4.add_subplot(312)
x8 = fig4.add_subplot(312)

# functions calling


grid_steps = np.array([100, 200, 400, 600, 800, 1600])
i = 0

for nx in grid_steps:
    nt = nx*2                   
    h = (xR - xL) / nx

    tau = r * h / a
    ksi = -a * tau

    x = np.linspace(xL, xR, nx)
    x_dence = np.linspace(xL, xR, nx)
    ux = np.zeros((nt, nx))

    u = np.zeros((nt, nx))
    u_dence = np.zeros(2 ** 10)
    s = np.zeros((nt, nx - 1))

    rho = np.zeros((nt, nx))
    del_rho = np.zeros((nt, nx))

    for i in range(0, nx):
        u[0][i] = test.test_func(x[i])
        # u[0][i]     = unit_box(x[i])
    for i in range(0, nx):
        u_dence[i] = test.test_func(x_dence[i])
        # u_dence[i]  =  unit_box(x_dence[i])
        # u[0][i] = fi(x[i])

        # РґР»СЏ РёРЅС‚РµРіСЂР°Р»Р° (РЅСѓР»РµРІРѕРіРѕ РјРѕРјРµРЅС‚Р°)
    for i in range(0, nx - 1):
        rho[0][i] = 0.5 * (u[0][i] + u[0][i + 1]) * h

    l = ods.upwind_scheme(u, nt, nx, r, x)
    print("upwind")
    print("N = ", nx)
    print(l)

    #x7.plot(x, u[0])
    #x8.plot(x, u[nt - 1])
    #l = ods.cubic_polynominal_interpolation(u, ux, nt, nx, h, ksi)
    print("CIP")
    print("N = ", nx)
    print(l)

    out_file.write(str(u))
    out_file.writelines(str(ux))
    #x1.plot(x, u[0])
    #x2.plot(x, u[nt - 1])

   # l = ods.fourth_order_polynominal_interpolation(u, ux, nt, nx, h, ksi, del_rho, rho)
    print("CIP CSL4")
    print("N = ", nx)
    print(l)


    l = ods.second_order_polynominal_interpolation(u, ux, nt, nx, h, ksi, del_rho, rho)
    print("CIP CSL2")
    print("N = ", nx)
    print(l)    # i = i + 1


outfile.close()
out_file.close()

plt.show()
