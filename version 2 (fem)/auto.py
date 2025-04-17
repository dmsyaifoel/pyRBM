from math import sin as msin, cos as mcos, acos as macos

def syms(s):
  return [Sym(name) for name in s.replace(' ', '').split(',')]

# def symvec(name, n):
#   return [Sym(name, j, n) for j in range(n)]

def sin(x):
    return Sin(x)

def cos(x):
    return Cos(x)

def acos(x):
    return Arccos(x)

class Sym:
    def __init__(self, name, j=None, n=None):
        self.j = j
        self.n = n
        if j is None:
            self.name = name
        else:
            self.name = name + str(j)

    def __repr__(self):
        return str(self.name)

    def __add__(self, other):
        if other == 0:
            return self
        return Add(self, other)

    def __radd__(self, other):
        if other == 0:
            return self
        return Add(other, self)

    def __mul__(self, other):
        if other == 1:
            return self
        return Mul(self, other)

    def __rmul__(self, other):
        if other == 1:
            return self
        return Mul(other, self)

    def __neg__(self):
        return Mul(-1, self)

    def __sub__(self, other):
        if isinstance(other, (int, float)):
            if other == 0:
                return self
            return Add(self, -other)
        return Add(self, Mul(-1, other))

    def __rsub__(self, other):
        return Add(-self, other)

    def __pow__(self, other):
        assert isinstance(other, (int, float))
        return Pow(self, other)

    def __truediv__(self, other):
        return Mul(self, Pow(other, -1))

    def val(self, dic):
        if self.n is None:
            val = dic[self.name]
            grad = [1 if self.name == key else 0 for key in dic]
        else:
            val = dic[self.j]
            grad = self.n * [0]
            grad[self.j] = 1
        return val, grad


class Add(Sym):
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def __repr__(self):
        return '(' + str(self.a) + ' + ' + str(self.b) + ')'

    def val(self, dic):
        if isinstance(self.a, (int, float)):
            va, ga = self.a, len(dic) * [0]
        else:
            va, ga = self.a.val(dic)
        if isinstance(self.b, (int, float)):
            vb, gb = self.b, len(dic) * [0]
        else:
            vb, gb = self.b.val(dic)
        val = va + vb
        grad = [ga[i] + gb[i] for i in range(len(dic))]
        return val, grad

class Mul(Sym):
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def __repr__(self):
        return '(' + str(self.a) + '*' + str(self.b) + ')'

    def val(self, dic):
        if isinstance(self.a, (int, float)):
            va, ga = self.a, len(dic) * [0]
        else:
            va, ga = self.a.val(dic)
        if isinstance(self.b, (int, float)):
            vb, gb = self.b, len(dic) * [0]
        else:
            vb, gb = self.b.val(dic)
        val = va * vb
        grad = [va * gb[i] + vb * ga[i] for i in range(len(dic))]
        return val, grad


class Pow(Sym):
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def __repr__(self):
        return '(' + str(self.a) + '**' + str(self.b) + ')'

    def val(self, dic):
        va, ga = self.a.val(dic)
        val = va ** self.b
        grad = [self.b * va ** (self.b - 1) * ga[i] for i in range(len(dic))]
        return val, grad


class Sin(Sym):
    def __init__(self, a):
        self.a = a

    def __repr__(self):
        return 'sin(' + str(self.a) + ')'

    def val(self, dic):
        va, ga = self.a.val(dic)
        val = msin(va)
        grad = [mcos(va) * ga[i] for i in range(len(dic))]
        return val, grad


class Cos(Sym):
    def __init__(self, a):
        self.a = a

    def __repr__(self):
        return 'cos(' + str(self.a) + ')'

    def val(self, dic):
        va, ga = self.a.val(dic)
        val = mcos(va)
        grad = [-msin(va) * ga[i] for i in range(len(dic))]
        return val, grad

class Arccos(Sym):
    def __init__(self, a):
        self.a = a

    def __repr__(self):
        return 'acos(' + str(self.a) + ')'

    def val(self, dic):
        va, ga = self.a.val(dic)
        val = macos(va)
        grad = [-ga[i] / (1 - va ** 2) ** .5 for i in range(len(dic))]
        return val, grad