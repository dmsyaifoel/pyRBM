from math import sin as _sin, cos as _cos, acos as _acos, pi as _pi

def isconst(x):
  return isinstance(x, (int, float))

def zeros(n):
  return n*[0]

class Sym:
  '''
  Base class for symbolic autodifferention
  '''
  def __init__(self, name):
    self.name = name
    self.key = hash(self)

  def __repr__(self):
    return str(self.name)

  def __add__(self, other):
    if other == 0:
      return self
    if isinstance(self, Add):
      if isinstance(other, Add):
        return Add(self.l + other.l, self.c + other.c)
      if isconst(other):
        return Add(self.l, self.c + other)
      return Add(self.l + [other], self.c)
    if isconst(other):
      return Add([self], other)
    return Add([self, other], 0)

  def __radd__(self, other):
    return self + (other)

  def __mul__(self, other):
    if other == 1:
      return self
    if other == 0:
      return 0
    if isinstance(self, Mul):
      if isinstance(other, Mul):
        return Mul(self.l + other.l, self.c*other.c)
      if isconst(other):
        return Mul(self.l, self.c * other)
      return Mul(self.l + [other], self.c)
    if isconst(other):
      return Mul([self], other)
    return Mul([self, other], 1)

  def __rmul__(self, other):
    return self*(other)

  def __neg__(self):
    return self*(-1)

  def __sub__(self, other):
    if other == 0:
      return self
    return self + -other

  def __rsub__(self, other):
    if other == 0:
      return self
    return -self + (other)

  def __pow__(self, other):
    if other == 0:
      return 1
    if other == 1:
      return self
    if isinstance(self, Pow):
      return self.a**(self.b*other)
    return Pow(self, other)

  def __truediv__(self, other):
    if other == 1:
      return self
    if isconst(other):
      return self*(1/other)
    return self*(other**-1)

  def val(self, dic):
    return dic[self.name]

  def grad(self, dic, n, cache):
    if self.key in cache:
      return cache[self.key]
    grad = zeros(n)
    if self.name in dic['map']:
      grad[dic['map'][self.name]] = 1
    val = dic[self.name]
    cache[self.key] = (val, grad)
    return val, grad

class Add(Sym):
  def __init__(self, l, c):
    self.l = l
    self.c = c
    self.key = hash(self)

  def __repr__(self):
    s = '('
    for i in self.l:
      s += str(i) + ' + '
    if self.c != 0:
      return s + str(self.c) + ')'
    return s[:-3] + ')'

  def val(self, dic):
    val = self.c
    for i in self.l:
      val += i.val(dic)
    return val

  def grad(self, dic, n, cache):
    if self.key in cache:
      return cache[self.key]
    val = self.c
    grad = zeros(n)
    for i in self.l:
      v_, g_ = i.grad(dic, n, cache)
      val += v_
      for i in range(n):
        grad[i] += g_[i]
    cache[self.key] = (val, grad)
    return val, grad

class Mul(Sym):
  def __init__(self, l, c):
    self.l = l
    self.c = c
    self.key = hash(self)

  def __repr__(self):
    s = '('
    for i in self.l:
      s += str(i) + '*'
    if self.c != 1:
      return s + str(self.c) + ')'
    return s[:-1] + ')'

  def val(self, dic):
    val = self.c
    for i in self.l:
      val *= i.val(dic)
    return val

  def grad(self, dic, n, cache):
    if self.key in cache:
      return cache[self.key]
    grad = zeros(n)
    vals = []
    grads = []
    val = self.c
    for i in self.l:
      v_, g_ = i.grad(dic, n, cache)
      val *= v_
      vals.append(v_)
      grads.append(g_)
    for k in range(len(self.l)):
      prod = self.c
      for j in range(len(self.l)):
        if j != k:
          prod *= vals[j]
      for i in range(n):
        grad[i] += grads[k][i]*prod
    cache[self.key] = (val, grad)
    return val, grad

class Pow(Sym):
  def __init__(self, a, b):
    self.a = a
    self.b = b
    self.key = hash(self)

  def __repr__(self):
    return '(' + str(self.a) + '**' + str(self.b) + ')'

  def val(self, dic):
    return self.a.val(dic)**self.b

  def grad(self, dic, n, cache):
    if self.key in cache:
      return cache[self.key]
    va, ga = self.a.grad(dic, n, cache)
    val = va**self.b
    grad = [self.b*va**(self.b - 1)*ga[i] for i in range(n)]
    cache[self.key] = (val, grad)
    return val, grad

class Sin(Sym):
  def __init__(self, a):
    self.a = a
    self.key = hash(self)

  def __repr__(self):
    return 'sin(' + str(self.a) + ')'

  def val(self, dic):
    return _sin(self.a.val(dic))

  def grad(self, dic, n, cache):
    if self.key in cache:
      return cache[self.key]
    va, ga = self.a.grad(dic, n, cache)
    val = _sin(va)
    grad = [_cos(va)*ga[i] for i in range(n)]
    cache[self.key] = (val, grad)
    return val, grad

class Cos(Sym):
  def __init__(self, a):
    self.a = a
    self.key = hash(self)

  def __repr__(self):
    return 'cos(' + str(self.a) + ')'

  def val(self, dic):
    return _cos(self.a.val(dic))

  def grad(self, dic, n, cache):
    if self.key in cache:
      return cache[self.key]
    va, ga = self.a.grad(dic, n, cache)
    val = _cos(va)
    grad = [-_sin(va)*ga[i] for i in range(n)]
    cache[self.key] = (val, grad)
    return val, grad

class Acos(Sym):
  def __init__(self, a):
    self.a = a
    self.key = hash(self)

  def __repr__(self):
    return 'acos(' + str(self.a) + ')'

  def val(self, dic):
    v = self.a.val(dic)
    if v >= 1:
      return 0
    elif v <= -1:
      return _pi
    return _acos(v)

  def grad(self, dic, n, cache):
    key = self.key
    if key in cache:
      return cache[key]
    va, ga = self.a.grad(dic, n, cache)
    if va >= 1:
      val, grad = 0, zeros(n)
    elif va <= -1:
      val, grad = _pi, zeros(n)
    else:
      val = _acos(va)
      grad = [-ga[i]/(1 - va**2)**.5 for i in range(n)]
    cache[key] = (val, grad)
    return val, grad

def sin(x):
  if isconst(x):
    return _sin(x)
  return Sin(x)

def cos(x):
  if isconst(x):
    return _cos(x)
  return Cos(x)

def acos(x):
  if isconst(x):
    return _acos(x)
  return Acos(x)

def sqrt(x):
  return x**.5

def sym(name):
  return Sym(name)

def syms(s):
  return [sym(name) for name in s.replace(' ', '').split(',')]

def symlist(name, n):
  return [sym(name + str(j)) for j in range(n)]

def parsedict(dic):
  dic2 = dict()
  dic3 = dict()
  for key, item in dic.items():
    if isconst(item):
      dic2[key] = item
    else:
      item = [float(i) for i in item]
      for i in range(len(item)):
        dic2[key + str(i)] = item[i]
        dic3[key + str(i)] = i
  if len(dic3) == 0:
    dic3 = {name: idx for idx, name in enumerate(dic.keys())}
  dic2['map'] = dic3
  return dic2