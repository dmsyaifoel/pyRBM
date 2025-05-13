backend = 'n'

if backend == 'n':
  import numpy as np
  array = np.array
  sin = np.sin
  cos = np.cos
  acos = np.acos
  import scipy.optimize as so
  minimize = so.minimize

if backend == 'c':
  import sym as sy
  sin = sy.sin
  cos = sy.cos
  acos = sy.acos
  from matrix import matrix as array
  from opt import minimize

if backend == 's':
  import sympy as sm
  array = sm.Matrix
  sin = sm.sin
  cos = sm.cos
  acos = sm.acos