import numpy as np

def test_func(x):
    if (x >= -0.8 and x <= -0.6):
        return (np.exp(-np.log1p(2) * (x + 0.7) ** 2 / 0.0009))
    else:
        if (x >= -0.4 and x <= -0.2):
            return 1
        else:
            if (x >= 0.0 and x <= 0.2):
                return (1 - np.fabs(10 * x - 1))
            else:
                if (x >= 0.4 and x <= 0.6):
                    return (1 - 100 * (x - 0.5) ** 2) ** 0.5
                else:
                    return 0



def fi(x):
    if (np.fabs(x) >= 0.1 and np.fabs(x) <= 0.3):
        return 1
    else:
        if ((np.fabs(x) >= 0 and np.fabs(x) < 0) or (np.fabs(x) > 0.3 and np.fabs(x) <= 1)):
            return 0

def unit_box(x):
    if(x>=-0.25 and x<=0.25):
        return 1
    else:
        return 0

def test_func_r(i):
    if(i>=13 and i<=21):
        return -1
    if (i >= 40 and i <= 48):
        return 1
    else:
        return 0



def test_func_sin(x):
    return np.sin(2*np.pi*x)

def der_test_func_sin(x):
    return 2 * np.pi * np.cos(2*np.pi*x)

def der2_test_func_sin(x):
    return -4 * np.pi**2 * np.sin(2*np.pi*x)

def der3_test_func_sin(x):
    return -8 * np.pi ** 3 * np.cos(2*np.pi*x)

def der4_test_func_sin(x):
    return 16 * np.pi ** 4 * np.sin(2*np.pi*x)

def two_dim_test_func_p(x, y, t):

    p = np.cos(2 * np.sqrt(2) * np.pi * t) * np.sin(2 * np.pi * x) * np.sin(2 * np.pi * y)

    return p


def two_dim_test_func_u(x, y, t):

    u = -1 / np.sqrt(2) * np.sin(2 * np.sqrt(2) * np.pi * t) * np.cos(2 * np.pi * x) * np.sin(2 * np.pi * y)

    return u


def two_dim_test_func_v(x, y, t):

    v = -1 / np.sqrt(2) * np.sin(2 * np.sqrt(2) * np.pi * t) * np.sin(2 * np.pi * x) * np.cos(2 * np.pi * y)

    return v


def two_dim_test_der_rp(x, y, t):

    rp = 2 * np.pi * np.sin(2 * np.pi * y) * np.cos(2 * np.pi * x) * np.cos(2 * np.sqrt(2) * np.pi * t)

    return rp


def two_dim_test_der_ru(x, y, t):

    ru = np.sqrt(2) * np.pi * np.sin(2 * np.pi * x) * np.sin(2 * np.pi * y) * np.sin(2 * np.sqrt(2) * np.pi * t)

    return ru


def two_dim_test_der_rv(x, y, t):

    rv = -np.sqrt(2) * np.pi * np.sin(2 * np.sqrt(2) * np.pi * t) * np.cos(2 * np.pi * x) * np.cos(2 * np.pi * y)

    return rv


def two_dim_test_der_sp(x, y, t):

    sp = 2 * np.pi * np.sin(2 * np.pi * x) * np.cos(2 * np.pi * y) * np.cos(2 * np.sqrt(2) * np.pi * t)

    return sp


def two_dim_test_der_su(x, y, t):

    su = -np.sqrt(2) * np.pi * np.sin(2 * np.sqrt(2) * np.pi * t) * np.cos(2 * np.pi * x) * np.cos(2 * np.pi * y)

    return su


def two_dim_test_der_sv(x, y, t):

    sv = np.sqrt(2) * np.pi * np.sin(2 * np.pi * x) * np.sin(2 * np.pi * y) * np.sin(2 * np.sqrt(2) * np.pi * t)

    return sv


def two_dim_test_exp(x, y):
    return np.exp(-(10 * (x ** 2 + (y + 0.7) ** 2)))

def two_dim_test_exp_r(x, y):
    return -20*x*np.exp(-10*x**2 - 10*(y + 0.7)**2)

def two_dim_test_exp_s(x, y):
    return (-20*y - 14.0)*np.exp(-10*x**2 - 10*(y + 0.7)**2)


def two_dim_test_sin(x, y):
    return 1 / np.pi * np.sin(np.pi * x) ** 2 * np.sin(np.pi * y) ** 2

def two_dim_test_sin_r(x, y):
    return 2 * np.sin(np.pi * x) * np.sin(np.pi * y) ** 2 * np.cos(np.pi * x)

def two_dim_test_sin_s(x, y):
    return 2 * np.sin(np.pi * x) ** 2 * np.sin(np.pi * y) * np.cos(np.pi * y)
