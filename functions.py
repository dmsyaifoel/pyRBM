import numpy as np
import numpy.linalg as nl

sin = np.sin
cos = np.cos
acos = np.acos
array = np.array
dot = np.dot
zeros = np.zeros
norm = nl.norm

def angle2(a, b, atol=1e-6):  # angle between vector a and b
  if norm(a - b) < atol:
    return 0
  return acos(dot(a, b) / norm(a) / norm(b))

angle = angle2

import scipy.optimize as so
minimize = so.minimize

def use_custom_instead():
  global sin, cos, acos, array, dot, zeros, norm, angle, minimize

  import sym
  sin = sym.sin
  cos = sym.cos
  acos = sym.acos

  import matrix
  array = matrix.matrix

  def dot2(a, b):
    return sum([a[i]*b[i] for i in range(len(a))])

  def zeros2(n):
    return n*[0]

  dot = dot2

  def norm2(v):
    return dot(v, v)**.5

  zeros = zeros2
  norm = norm2

  def angle2(a, b):
    return acos(dot(a, b)/norm(a)/norm(b))

  angle = angle2

# 2d rotation matrix for angle t
def R(t):
  return array([[cos(t), -sin(t)], [sin(t), cos(t)]])

# rotation matrix for list of angles in 3d, or one angle in 2d
def rotmat(angles, dim):
  if dim == 2: return R(angles)
  if dim == 3: return RX(angles[0])@RY(angles[1])@RZ(angles[2])

def RX(t): # 3d rotation matrix around x
  return array([[1, 0, 0], [0, cos(t), -sin(t)], [0, sin(t), cos(t)]])

def RY(t): # 3d rotation matrix around y
  return array([[cos(t), 0, sin(t)], [0, 1, 0], [-sin(t), 0, cos(t)]])

def RZ(t): # 3d rotation matrix around z
  return array([[cos(t), -sin(t), 0], [sin(t), cos(t), 0], [0, 0, 1]])