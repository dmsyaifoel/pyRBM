from sym import symlist, sin, cos, parsedic
from opt import root_scalar, minimize_scalar, line_search

x = symlist('x', 10)

y = sum([100*(x[i+1] - x[i]**2)**2 + (1 - x[i])**2 for i in range(9)])

def g(x):
  return y.valgrad(parsedic({'x':list(x)}), 10, {})[1]

def f(x):
  return y.val(parsedic({'x':list(x)}))

from scipy.optimize import minimize
print(minimize(f, 10*[0], jac=g))
