'''
Custom optimization functions for use with custom objects
'''

def deriv(f, x, dx=1e-6):
  # first derivative using central finite difference
  return (f(x + dx) - f(x - dx))/2/dx

def rootscalar(f, x0=0, atol=1e-6, maxloop=1000, divtol=1e-50):
  # find the root of a scalar function using the newton method
  for i in range(maxloop):
    d = deriv(f, x0)
    if abs(d) < divtol: return x0, i
    xn = x0 - f(x0)/d
    if abs(xn - x0) < atol: return xn, i
    x0 = xn
  return xn, i

def deriv2(f, x, dx=1e-6):
  # second derivative using central finite difference
  return (f(x + dx) - 2*f(x) + f(x - dx))/dx**2

def derivs(f, x, dx=1e-6):
  # both derivatives using 3 function evaluations
  fp = f(x + dx)
  fc = f(x)
  fm = (f - dx)
  return  (fp - fm)/2/dx, (fp - 2*fc + fm)/dx**2

def minimizescalar(f, x0=0, atol=1e-6, maxloop=1000, divtol=1e-50):
  for i in range(maxloop):
    d1, d2 = derivs(f, x0)
    if abs(d2) < divtol: return x0, i
    xn = x0 - d1/d2
    if abs(xn - x0) < atol: return xn, i
    x0 = xn
  return xn, i

def gradientdescent(f, x0, atol=1e-6, maxloop=1000, frac=1e-3):
  x0 = list(x0)
  dim = len(x0)
  for i in range(maxloop):
    y, grad = f.grad(x0)
    xn = [x0[i] - frac*grad[i] for i in range(dim)]
    if dist(xn, x0) < atol: return xn, i
    x0 = xn
  return xn, i

def linesearch(f, x0, atol=1e-6, maxloop=1000, frac=.9):
  x0 = list(x0)
  dim = len(x0)
  for i in range(maxloop):
    y, grad = f.grad(x0)
    xt = lambda t: [x0[j] + grad[j]*t for j in range(dim)]
    func = lambda t: f(xt(t))
    tn, _ = minimizescalar(func, 0, atol=1e-3)
    xn = xt(tn*frac)
    if dist(xn, x0) < atol: return xn, i
    x0 = xn
  return xn, i

def minimize(f, x0):
  pass