from __future__ import unicode_literals
import one_dimentional_solver as ods
import matplotlib
pgf_with_lualatex = {
    "pgf.texsystem": "lualatex",
    "pgf.preamble": [
         r"\usepackage[utf8]{inputenc}",
         r"\usepackage[russian]{babel}"
         ]
}
from matplotlib import rc
from matplotlib import pyplot as plt
matplotlib.rcParams.update(pgf_with_lualatex)

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True

import numpy as np



rc('font',**{'family':'serif'})
plt.rc('text', usetex=True)
plt.rc('text.latex', unicode=True)
# rc('text.latex',preamble='\\usepackage[utf8]{inputenc}')
# rc('text.latex',preamble='\')

a, xL, xR, tfin = 1, -1, 1, 1.0

# CFL
r = 0.5

# functions calling

grid_steps = np.array([8, 16, 32, 64])
# grid_steps = np.array([100, 200, 400, 800, 1600])

# fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True)


# def init():  # only required for blitting to give a clean slate.
#     line.set_ydata([np.nan] * len(x))
#     return line,
#
#
# def animate(i):
#     line.set_ydata(u[i])  # update the data.
#     return line,
c_cip3 = np.zeros(grid_steps.size)
c_cip5 = np.zeros(grid_steps.size)
c_cip7 = np.zeros(grid_steps.size)
c_cip9 = np.zeros(grid_steps.size)

c_deltaP3 = np.zeros(grid_steps.size)
c_deltaP5 = np.zeros(grid_steps.size)
c_deltaP7 = np.zeros(grid_steps.size)
c_deltaP9 = np.zeros(grid_steps.size)


# cL_cip3 = np.zeros(grid_steps.size)
# cL_cip5 = np.zeros(grid_steps.size)
# cL_cip7 = np.zeros(grid_steps.size)
# cL_cip9 = np.zeros(grid_steps.size)

# cL_deltaP3 = np.zeros(grid_steps.size)
# cL_deltaP5 = np.zeros(grid_steps.size)
# cL_deltaP7 = np.zeros(grid_steps.size)
# cL_deltaP9 = np.zeros(grid_steps.size)


# c_deltaP3_2D = np.zeros(grid_steps.size)

for i in range(grid_steps.size):
    nx = grid_steps[i]
    h = (xR - xL) / nx
    tau = r * h / a
    xi = -a * tau
    x = np.linspace(xL, xR, nx + 1)
    error = np.zeros((8, 3))
    err = np.zeros(3)
    nt = int(nx / r)
    print('h = {0}, nx = {1}, tau = {2}, nt = {3}'.format(h, nx, tau, nt))
    res_cip3, ua3 = ods.cip3(x, tau, h, xi, tfin, error[0])
    res_cip5, ua5 = ods.cip5(x, tau, h, xi, tfin, nt, error[1])
    res_cip7 = ods.cip7(x, tau, nt, h, xi, tfin, error[2])
    res_cip9 = ods.cip9(x, tau, nt, h, xi, tfin, error[3])
    res_deltaP3, ua = ods.deltaP3sym(x, r, tau, h, tfin, error[4])
    res_deltaP5, ua = ods.deltaP5sym(x, r, tau, h, tfin, error[5])
    res_deltaP7, ua = ods.deltaP7sym(x, r, tau, h, tfin, error[6])
    res_deltaP9, ua = ods.deltaP9sym(x, r, tau, h, tfin, error[7])
#     res_deltaP3_2D, pa = ods.two_dim_deltaP3(err, grid_steps[i], grid_steps[i])
#     print(err[0])



   # c_deltaP3_2D[i] = err[0]
#     # print(res_rcip)
    c_cip3[i] = error[0][0]
    c_cip5[i] = error[1][0]
    c_cip7[i] = error[2][0]
    c_cip9[i] = error[3][0]
    c_deltaP3[i] = error[4][0]
    c_deltaP5[i] = error[5][0]
    c_deltaP7[i] = error[6][0]
    c_deltaP9[i] = error[7][0]
    #
    # cL_cip3[i] = error[0][2]
    # cL_cip5[i] = error[1][2]
    # cL_cip7[i] = error[2][2]
    # cL_cip9[i] = error[3][2]

    # cL_deltaP3[i] = error[4][2]
    # cL_deltaP5[i] = error[5][2]
    # cL_deltaP7[i] = error[6][2]
    # cL_deltaP9[i] = error[7][2]

convergence_plot_grid = np.zeros(grid_steps.size - 1)
# conv2D = ods.convergence(c_deltaP3_2D)
#
conv3 = ods.convergence(c_cip3)
conv5 = ods.convergence(c_cip5)
conv7 = ods.convergence(c_cip7)
conv9 = ods.convergence(c_cip9)

print("L DeltaP")

convP3 = ods.convergence(c_deltaP3)
convP5 = ods.convergence(c_deltaP5)
convP7 = ods.convergence(c_deltaP7)
convP9 = ods.convergence(c_deltaP9)

# print("Linf")
#
# Lconv3 = ods.convergence(cL_cip3)
# Lconv5 = ods.convergence(cL_cip5)
# Lconv7 = ods.convergence(cL_cip7)
# Lconv9 = ods.convergence(cL_cip9)


# print("Linf DeltaP")
#
#
# LconvdeltaP3 = ods.convergence(cL_deltaP3)
# LconvdeltaP3 = ods.convergence(cL_deltaP5)
# LconvdeltaP3 = ods.convergence(cL_deltaP7)
# LconvdeltaP3 = ods.convergence(cL_deltaP9)


# print('CIP3 {0}'.format(c_cip3))
# print('CIP5 {0}'.format(c_cip5))
# print('CIP7 {0}'.format(c_cip7))
# print('CIP9 {0}'.format(c_cip9))

print('deltaP3 {0}'.format(c_deltaP3))

print('deltaP5 {0}'.format(c_deltaP5))
print('deltaP7 {0}'.format(c_deltaP7))
print('deltaP9 {0}'.format(c_deltaP9))

# print('deltaP3 {0}'.format(c_deltaP3_2D))

print("Linf")

# print('CIP3 {0}'.format(cL_cip3))
# print('CIP5 {0}'.format(cL_cip5))
# print('CIP7 {0}'.format(cL_cip7))
# print('CIP9 {0}'.format(cL_cip9))
# print('deltaP3 {0}'.format(cL_deltaP3))
# print('deltaP5 {0}'.format(cL_deltaP5))
# print('deltaP7 {0}'.format(cL_deltaP7))
# print('deltaP9 {0}'.format(cL_deltaP9))
#
convergence_plot3 = np.zeros(grid_steps.size - 1)

convergence_plot5 = np.zeros(grid_steps.size - 1)

convergence_plot7 = np.zeros(grid_steps.size - 1)

convergence_plot9 = np.zeros(grid_steps.size - 1)

convergence_plotP3 = np.zeros(grid_steps.size - 1)

convergence_plotP5 = np.zeros(grid_steps.size - 1)

convergence_plotP7 = np.zeros(grid_steps.size - 1)

convergence_plotP9 = np.zeros(grid_steps.size - 1)

# convergence_plotP3_2D = np.zeros(grid_steps.size - 1)


for j in range(1, grid_steps.size):
    # convergence_plotP3_2D[j - 1] = conv2D[j]
    convergence_plot3[j - 1] = conv3[j]
    convergence_plot5[j - 1] = conv5[j]
    convergence_plot7[j - 1] = conv7[j]
    convergence_plot9[j - 1] = conv9[j]
    convergence_plotP3[j - 1] = convP3[j]
    convergence_plotP5[j - 1] = convP5[j]
    convergence_plotP7[j - 1] = convP7[j]
    convergence_plotP9[j - 1] = convP9[j]

    convergence_plot_grid[j - 1] = j

# plt.subplot(121)
plt.plot(x, res_deltaP9, 'o-', label='$\Delta$-P9')
plt.plot(x, ua, '-')
plt.ylabel(r'u(x, t)')
# plt.xlabel(r'x')
plt.xlim(x.min(), x.max())
plt.ylim(-1, 1)
# plt.title(u'а)')
# plt.legend()
# plt.subplot(122)
# plt.plot(convergence_plot_grid, convergence_plotP9, '-o', label=u'порядок сходимости')
# plt.xlim(convergence_plot_grid.min() - 0.01, convergence_plot_grid.max() + 0.01)
# plt.ylim(convergence_plotP9.min() - 0.01, convergence_plotP9.max() + 0.01)
# plt.title(u'б)')
plt.legend()
plt.grid()
# plt.savefig("/Users/ramazanova/Documents/master_diss/pic/delta_P9.pdf", dpi=180)
plt.show()
