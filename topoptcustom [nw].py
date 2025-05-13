import numpy as np
from scipy.sparse import coo_array
from scipy.sparse.linalg import spsolve

# MAKE THIS STIFFNESS MATRIX OF A SINGLE FEM ELEMENT
# DON'T QUESTION THIS
# IT IS ANCIENT MAGIC AND I DON'T KNOW HOW IT WORKS
E = 1.0
nu = 0.3
k = np.array([1 / 2 - nu / 6, 1 / 8 + nu / 8, -1 / 4 - nu / 12, -1 / 8 + 3 * nu / 8,
              -1 / 4 + nu / 12, -1 / 8 - nu / 8, nu / 6, 1 / 8 - 3 * nu / 8])
K_elem = E / (1 - nu ** 2) * np.array([
  [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
  [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
  [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
  [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
  [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
  [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
  [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
  [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]
])

# ASSEMBLE K_elem INTO A BIG STIFFNESS MATRIX FOR THE WHOLE THING
# IT IS A SPARSE MATRIX
# WHICH IS ZERO EVERYWHERE, EXCEPT AT (i, j), where it has value s
# so we find these triples to set up the K
def assemble_K(nelx, nely, x, penal):
  ndof = 2*(nelx + 1)*(nely + 1)
  iK = np.zeros(8 * 8 * nelx * nely, dtype=int)
  jK = np.zeros(8 * 8 * nelx * nely, dtype=int)
  sK = np.zeros(8 * 8 * nelx * nely)
  # we will fill these vectors to make this matrix

  count = 0
  for elx in range(nelx):
    for ely in range(nely):
      n1 = (nely + 1) * elx + ely
      n2 = (nely + 1) * (elx + 1) + ely
      edof = np.array([2 * n1, 2 * n1 + 1, 2 * n2, 2 * n2 + 1, 2 * n2 + 2, 2 * n2 + 3, 2 * n1 + 2, 2 * n1 + 3])
      # edof is an array of 8 long which denotes which degree of freedom of the
      # corners of the (elx, ely)th square
      # in other words, the first square has degrees of freedom at points
      # 0  1 82 83 84 85  2  3
      # in the state vector (the vector which contains all the volume fractions)
      for i in range(8):
        for j in range(8):
          iK[count] = edof[i]
          jK[count] = edof[j]
          sK[count] = K_elem[i, j] * x[ely, elx] ** penal
          count += 1
          # this fills the iK, jK, and sK arrays
  print('worked')

  # then we have all these triples so we can fill the K matrix
  K = coo_array((sK, (iK, jK)), shape=(ndof, ndof)).tocsc()
  # it is converted to CSC form
  return K


nelx, nely = 20, 20
volfrac = .5
penal = 3
rmin = 1.5

ndof = 2 * (nelx + 1) * (nely + 1)

x = volfrac * np.ones((nely, nelx))
loop = 0
change = 1

f = np.zeros(ndof)
f[1] = -1.0
# this is the array that contains all the forces
# now it looks like [0, -1, 0, ..., 0]
# with order x0, y0, x1, y1, x2, y2, ...
# so the first node has a downward force and here are no other forces

K = assemble_K(nelx, nely, penal)

fixed_dofs = np.union1d(np.arange(0, 2 * (nely + 1), 2), [2 * (nelx + 1) * (nely + 1) - 1])
# this is just a list of all the dof numbers that are fixed
# in this case this is a lot of small even numbers and one big odd numbers
# this means that all nodes along the left edge are constrained only in x direction but not in y(?)
# whereas one node is constrained only in y direction (so loads of rollers)

free_dofs = np.setdiff1d(np.arange(ndof), fixed_dofs)
# all the other dofs are free

while change > .01:
  loop += 1
  xold = x.copy()

  u = spsolve(K[free_dofs, :][:, free_dofs], f[free_dofs])
  # we solve Ku = f for u
  break




  # U = FE(nelx, nely, x, penal)