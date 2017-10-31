import one_dimentional_solver as ods
import test_functions as test
from matplotlib import pyplot as plt

import numpy as np

a = 1
l = 1
T = 0.5

xL = -1
xR = 1

# CFL
r = 0.4

# output files
out_file = open('results.dat', 'w')

# plot

# fig = plt.figure()
# #
# x1 = fig.add_subplot(111)
# x2 = fig.add_subplot(111)
# #
# fig2 = plt.figure()
# #
# x3 = fig2.add_subplot(111)
# x4 = fig2.add_subplot(111)

fig3 = plt.figure()
#
x5 = fig3.add_subplot(111)
x6 = fig3.add_subplot(111)

# fig4 = plt.figure()
#
# x7 = fig4.add_subplot(111)
# x8 = fig4.add_subplot(111)
#
# fig5 = plt.figure()
#
# x9 = fig5.add_subplot(111)
# x10 = fig5.add_subplot(111)

# functions calling


nx = 200

nt = 20
h = (xR - xL) / nx

tau = r * h / a
xi = -a * tau

x = np.linspace(xL, xR, nx)

u = np.zeros((nt, nx))
ux = np.zeros((nt, nx))
u2x = np.zeros((nt, nx))
u3x = np.zeros((nt, nx))

for i in range(0, nx):
    u[0][i] = test.test_func_sin(x[i])

# l = ods.upwind_scheme(u, nt, nx, r, h)
#
# print("upwind")
# print("N = ", nx)
# print(l)
#
# x1.plot(x, u[0])
# x2.plot(x, u[nt - 1])
#
# l = ods.cubic_polynominal_interpolation(u, ux, nt, nx, h, xi)
#
# print("CIP")
# print("N = ", nx)
# print(l)
#
# x3.plot(x, u[0])
# x4.plot(x, u[nt - 1])
l = ods.cip3(u, ux, nt, nx, h, xi)
x5.plot(x, u[0])
x6.plot(x, u[nt - 1])
print(l)
# print(u)
# l = ods.cip5(u, ux, u2x, nt, nx, h, xi)
# x7.plot(x, u[0])
# x8.plot(x, u[nt - 1])
#
# l = ods.cip7(u, ux, u2x, u3x, nt, nx, h, xi)
#
# x9.plot(x, u[0])
# x8.plot(x, u[nt - 1])

out_file.close()

plt.show()
