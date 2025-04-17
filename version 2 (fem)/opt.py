def dist(x1, x2):
  return sum([(x1[i] - x2[i])**2 for i in range(len(x1))])**.5

def gradientdescent(f, x0, atol=1e-6, maxloop=1000, frac=.01):
  x0 = list(x0)
  dim = len(x0)
  for i in range(maxloop):
    y, grad = f.val(x0)
    print(y, grad)
    xn = [x0[i] - frac*grad[i] for i in range(dim)]
    if dist(xn, x0) < atol: return xn, i
    x0 = xn
  return xn, 'maxloop'

def deriv(f, x, dx=1e-6): return (f(x + dx) - f(x - dx))/2/dx

def newtonroot(f, x0=0, atol=1e-6, maxloop=1000, divtol=1e-50):
  for i in range(maxloop):
    d = deriv(f, x0)
    if abs(d) < divtol: return x0, 'divtol'
    xn = x0 - f(x0)/d
    if abs(xn - x0) < atol: return xn, i
    x0 = xn
  return xn, 'maxloop'

def deriv2(f, x, dx=1e-6): return (f(x + dx) - 2*f(x) + f(x - dx))/dx**2

def newtonmin(f, x0=0, atol=1e-6, maxloop=1000, divtol=1e-50):
  for i in range(maxloop):
    d2 = deriv2(f, x0)
    if abs(d2) < divtol: return x0, 'divtol'
    xn = x0 - deriv(f, x0)/d2
    if abs(xn - x0) < atol: return xn, i
    x0 = xn
  return xn, 'maxloop'

def linesearch(f, x0, atol=1e-6, maxloop=1000, frac=.9):
  x0 = list(x0)
  dim = len(x0)
  for i in range(maxloop):
    y, grad = f.val(x0)
    xt = lambda t: [x0[j] + grad[j]*t for j in range(dim)]
    func = lambda t: f(xt(t))
    tn, _ = newtonmin(func, 0)
    xn = xt(tn*frac)
    if dist(xn, x0) < atol: return xn, i
    x0 = xn
  return xn, 'maxloop'