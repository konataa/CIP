import math

def test_func(x):
    if (x >= -0.8 and x <= -0.6):
        return (math.exp(-math.log1p(2) * (x + 0.7) ** 2 / 0.0009))
    else:
        if (x >= -0.4 and x <= -0.2):
            return 1
        else:
            if (x >= 0.0 and x <= 0.2):
                return (1 - math.fabs(10 * x - 1))
            else:
                if (x >= 0.4 and x <= 0.6):
                    return (1 - 100 * (x - 0.5) ** 2) ** 0.5
                else:
                    return 0



def fi(x):
    if (math.fabs(x) >= 0.1 and math.fabs(x) <= 0.3):
        return 1
    else:
        if ((math.fabs(x) >= 0 and math.fabs(x) < 0) or (math.fabs(x) > 0.3 and math.fabs(x) <= 1)):
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
