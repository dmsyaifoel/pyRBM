from matrix import Matrix
from auto import sin, cos, acos

def dot(a, b):
  return a.dot(b)

def norm(v):
  return dot(v, v)

def R(t):
  return Matrix([[cos(t), -sin(t)], [sin(t), cos(t)]])

def angle(a, b):
  return acos(dot(a, b) / norm(a) / norm(b))