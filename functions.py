import numpy as np
import scipy.linalg as sl
import scipy.optimize as so
import matplotlib.pyplot as mp
import plotly.graph_objects as pg
import multiprocessing as mu
import copy as cp

# rotate vector v by an angle theta over a unit vector k
def rod(v, theta, k):
  return v*np.cos(theta) + np.cross(k, v)*np.sin(theta) + k*np.dot(k, v)*(1 - np.cos(theta))

# find the rotation matrix which rotates a to b
def Rab(a, b):
  v = np.cross(a, b)
  s = sl.norm(v)
  c = np.dot(a, b)
  v_ = np.array([[0, -v[2], v[1]],
                 [v[2], 0, -v[0]],
                 [-v[1], v[0], 0]])
  return np.eye(3) + v_ + v_**2/(1 + c)

# 2d rotation matrix for angle t
def R(t): return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])

# rotation matrix for list of angles in 3d, or one angle in 2d
def rotmat(angles, dim):
  if dim == 2: return R(angles)
  if dim == 3: return RX(angles[0])@RY(angles[1])@RZ(angles[2])

def RX(t): # 3d rotation matrix around x
  return np.array([[1, 0, 0],
                   [0, np.cos(t), -np.sin(t)],
                   [0, np.sin(t), np.cos(t)]])

def RY(t): # 3d rotation matrix around y
  return np.array([[np.cos(t), 0, np.sin(t)],
                   [0, 1, 0],
                   [-np.sin(t), 0, np.cos(t)]])

def RZ(t): # 3d rotation matrix around z
  return np.array([[np.cos(t), -np.sin(t), 0],
                   [np.sin(t), np.cos(t), 0],
                   [0, 0, 1]])

def angle(a, b, atol=1e-6): # angle between vector a and b
  if sl.norm(a - b) < atol: return 0
  return np.arccos(np.clip(np.dot(a, b)/sl.norm(a)/sl.norm(b), -1, 1))

# turn vector into unit vector
def normalize(v): return v/sl.norm(v)

# wrapper for the energy function for multiprocessing
def poolenergy(prbm, end_effector, free_bodies, position, displacement):
  copied_prbm = cp.deepcopy(prbm)
  dmove = displacement[:copied_prbm.dim]
  if copied_prbm.dim == 2: drot = displacement[-1]
  else: drot = displacement[copied_prbm.dim:]
  copied_prbm.move(end_effector, position + dmove, drot)
  copied_prbm.solve_pose(free_bodies)
  return copied_prbm.energy()
