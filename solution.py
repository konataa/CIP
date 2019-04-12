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

c_cip3 = np.zeros(grid_steps.size)


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
    c_cip3[i] = error[0][0]

convergence_plot_grid = np.zeros(grid_steps.size - 1)
conv3 = ods.convergence(c_cip3)

convergence_plot3 = np.zeros(grid_steps.size - 1)

for j in range(1, grid_steps.size):
    # convergence_plotP3_2D[j - 1] = conv2D[j]
    convergence_plot3[j - 1] = conv3[j]
    convergence_plot_grid[j - 1] = j

# plt.subplot(121)
plt.plot(x, res_cip3, 'o-', label='CIP3')
plt.plot(x, ua3, '-')
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
plt.show()
