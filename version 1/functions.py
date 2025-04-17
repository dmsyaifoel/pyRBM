import numpy as np
import scipy.linalg as sl
import scipy.optimize as so
import matplotlib.pyplot as mp
import plotly.graph_objects as pg

def rod(v, theta, k):
  return v*np.cos(theta) + np.cross(k, v)*np.sin(theta) + k*np.dot(k, v)*(1 - np.cos(theta))

def Rab(a, b):
  v = np.cross(a, b)
  s = sl.norm(v)
  c = np.dot(a, b)
  v_ = np.array([[0, -v[2], v[1]],
                 [v[2], 0, -v[0]],
                 [-v[1], v[0], 0]])
  return np.eye(3) + v_ + v_**2/(1 + c)

def R(t):
  r = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
  return r

def rotmat(angles, dim):
  if dim == 2:
    return R(angles)
  if dim == 3:
    return RX(angles[0])@RY(angles[1])@RZ(angles[2])

def RX(t):
  return np.array([[1, 0, 0],
                   [0, np.cos(t), -np.sin(t)],
                   [0, np.sin(t), np.cos(t)]])

def RY(t):
  return np.array([[np.cos(t), 0, np.sin(t)],
                   [0, 1, 0],
                   [-np.sin(t), 0, np.cos(t)]])

def RZ(t):
  return np.array([[np.cos(t), -np.sin(t), 0],
                   [np.sin(t), np.cos(t), 0],
                   [0, 0, 1]])

def angle(a, b, atol=1e-6):
  if sl.norm(a - b) < atol: return 0
  return np.arccos(np.clip(np.dot(a, b)/sl.norm(a)/sl.norm(b), -1, 1))

def normalize(v):
  return v/sl.norm(v)
