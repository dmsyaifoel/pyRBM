backend = 'n'

if backend == 'n':
  import numpy as np
  zeros = np.zeros
  array = np.array
  sin = np.sin
  cos = np.cos
  acos = np.acos
  dot = np.dot

  import numpy.linalg as nl
  norm = nl.norm
  import scipy.optimize as so
  minimize = so.minimize

if backend == 'c':
  import sym as sy
  zeros = sy.zeros
  sin = sy.sin
  cos = sy.cos
  acos = sy.acos
  from matrix import matrix as array
  from opt import linesearch as minimize

if backend == 's':
  import sympy as sm
  zeros = lambda n: sm.Matrix([0]*n)
  array = sm.Matrix
  sin = sm.sin
  cos = sm.cos
  acos = sm.acos