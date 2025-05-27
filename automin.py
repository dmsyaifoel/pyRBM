from math import sin as _sin, cos as _cos, acos as _acos, pi as _pi, sqrt as _sqrt, ceil as _ceil, floor as _floor

'''
A rudimentary Interval Arithmetic implementation. Intended to determine the minimum and maximum of a function upon
symbolic buildup, and to be integrated with sym.py
'''

class Automin:
    def __init__(self, min_val, max_val):
        self.min_val = float(min_val)
        self.max_val = float(max_val)

    def __add__(self, other):
        if not isinstance(other, Automin):
            other = Automin(other)

        new_min = self.min_val + other.min_val
        new_max = self.max_val + other.max_val
        return Automin(new_min, new_max)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if not isinstance(other, Automin):
            other = Automin(other)

        new_min = self.min_val - other.max_val
        new_max = self.max_val - other.min_val
        return Automin(new_min, new_max)

    def __rsub__(self, other):
        return Automin(other) - self

    def __mul__(self, other):
        if not isinstance(other, Automin):
            other = Automin(other)

        products = [
            self.min_val*other.min_val,
            self.min_val*other.max_val,
            self.max_val*other.min_val,
            self.max_val*other.max_val
        ]
        new_min = min(products)
        new_max = max(products)
        return Automin(new_min, new_max)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if not isinstance(other, Automin):
            other = Automin(other)

        if other.min_val <= 0 <= other.max_val:
            raise ValueError(f'Cannot divide by an interval which includes zero.')

        reciprocal_min = 1/other.max_val
        reciprocal_max = 1/other.min_val

        products = [
            self.min_val*reciprocal_min,
            self.min_val*reciprocal_max,
            self.max_val*reciprocal_min,
            self.max_val*reciprocal_max
        ]
        new_min = min(products)
        new_max = max(products)
        return Automin(new_min, new_max)

    def __rtruediv__(self, other):
        return Automin(other)/self

    def __pow__(self, power):
        if not isinstance(power, (int, float)):
            raise TypeError('Power must be a scalar (int or float).')

        if power == 0:
            return Automin(1.0, 1.0)

        if power%1 == 0:
            if power > 0:
                if self.min_val >= 0:
                    new_min = self.min_val**power
                    new_max = self.max_val**power
                elif self.max_val <= 0:
                    if power%2 == 1:
                        new_min = self.min_val**power
                        new_max = self.max_val**power
                    else:
                        new_min = self.max_val**power
                        new_max = self.min_val**power
                else:
                    if power%2 == 1:
                        new_min = self.min_val**power
                        new_max = self.max_val**power
                    else:
                        new_min = 0
                        new_max = max(self.min_val**power, self.max_val**power)
            else:
                if self.min_val <= 0 <= self.max_val:
                    raise ValueError("Cannot raise to a negative power if the interval includes zero.")

                reciprocal_min = 1/self.max_val
                reciprocal_max = 1/self.min_val
                temp_automin = Automin(reciprocal_min, reciprocal_max)

                return temp_automin**abs(power)
        else:
            if self.min_val < 0:
                raise ValueError("Cannot raise negative numbers to fractional powers for real results.")

            new_min = self.min_val ** power
            new_max = self.max_val ** power

        return Automin(new_min, new_max)

    def __neg__(self):
        new_min = -self.max_val
        new_max = -self.min_val
        return Automin(new_min, new_max)

    def __repr__(self):
        return f'{[self.min_val, self.max_val]}'

def automin_sqrt(x_automin):
    if not isinstance(x_automin, Automin):
        x_automin = Automin(x_automin)
    return x_automin**.5

def automin_sin(x_automin):
    if not isinstance(x_automin, Automin):
        x_automin = Automin(x_automin)

    if x_automin.max_val - x_automin.min_val >= 2*_pi:
        new_min = -1
        new_max = 1
    else:
        val1 = _sin(x_automin.min_val)
        val2 = _sin(x_automin.max_val)

        temp_min = min(val1, val2)
        temp_max = max(val1, val2)

        k_start = int(_floor(x_automin.min_val/(_pi/2)))
        k_end = int(_ceil(x_automin.max_val/(_pi/2)))

        for k in range(k_start, k_end + 1):
            crit_point_val = k*(_pi/2)
            if x_automin.min_val <= crit_point_val <= x_automin.max_val:
                if _sin(crit_point_val) == 1:
                    temp_max = 1
                elif _sin(crit_point_val) == -1:
                    temp_min = -1

        new_min = temp_min
        new_max = temp_max

    return Automin(new_min, new_max)

def automin_cos(x_automin):
    if not isinstance(x_automin, Automin):
        x_automin = Automin(x_automin)

    if x_automin.max_val - x_automin.min_val >= 2*_pi:
        new_min = -1
        new_max = 1
    else:
        val1 = _cos(x_automin.min_val)
        val2 = _cos(x_automin.max_val)

        temp_min = min(val1, val2)
        temp_max = max(val1, val2)

        k_start = int(_floor(x_automin.min_val/(_pi/2)))
        k_end = int(_ceil(x_automin.max_val/(_pi/2)))

        for k in range(k_start, k_end + 1):
            crit_point_val = k*(_pi/2)
            if x_automin.min_val <= crit_point_val <= x_automin.max_val:
                if _cos(crit_point_val) == 1:
                    temp_max = 1
                elif _cos(crit_point_val) == -1:
                    temp_min = -1

        new_min = temp_min
        new_max = temp_max

    return Automin(new_min, new_max)

def automin_acos(x_automin):
    if not isinstance(x_automin, Automin):
        x_automin = Automin(x_automin)

    if x_automin.max_val < -1.0 or x_automin.min_val > 1.0:
        raise ValueError(f"Domain error for acos: interval [{x_automin.min_val}, {x_automin.max_val}] is outside [-1, 1].")

    clamped_min = max(-1.0, x_automin.min_val)
    clamped_max = min(1.0, x_automin.max_val)

    new_value = _acos(x_automin.value)

    new_min = _acos(clamped_max)
    new_max = _acos(clamped_min)

    return Automin(new_value, new_min, new_max)
