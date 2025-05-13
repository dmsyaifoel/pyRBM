import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def top(nelx, nely, volfrac, penal, rmin):
  x = volfrac * np.ones((nely, nelx))
  loop = 0
  change = 1.0

  # START ITERATION
  while change > 0.01:
    loop += 1
    xold = x.copy()

    # FE-ANALYSIS
    U = FE(nelx, nely, x, penal)

    # OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    KE = lk()
    c = 0.0
    dc = np.zeros((nely, nelx))

    for elx in range(nelx):
      for ely in range(nely):
        n1 = (nely + 1) * elx + ely
        n2 = (nely + 1) * (elx + 1) + ely
        edof = np.array([2 * n1, 2 * n1 + 1, 2 * n2, 2 * n2 + 1, 2 * n2 + 2, 2 * n2 + 3, 2 * n1 + 2, 2 * n1 + 3])
        Ue = U[edof].flatten()
        c += x[ely, elx] ** penal * Ue.dot(KE.dot(Ue))
        dc[ely, elx] = -penal * x[ely, elx] ** (penal - 1) * Ue.dot(KE.dot(Ue))

    # FILTERING OF SENSITIVITIES
    dc = check(nelx, nely, rmin, x, dc)

    # DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    x = OC(nelx, nely, x, volfrac, dc)

    # PRINT RESULTS
    change = np.max(np.abs(x - xold))
    print(f" It.: {loop:4d} Obj.: {c:10.4f} Vol.: {np.sum(x) / (nelx * nely):6.3f} ch.: {change:6.3f}")

    # PLOT DENSITIES
    plt.figure(1)
    plt.clf()
    plt.imshow(-x, cmap='gray', interpolation='none')
    plt.axis('equal')
    plt.axis('off')
    plt.pause(0.001)

  plt.ioff()
  plt.show()


def OC(nelx, nely, x, volfrac, dc):
  """
  Optimality criteria update
  """
  l1 = 0
  l2 = 1e5
  move = 0.2

  while (l2 - l1 > 1e-4):
    lmid = 0.5 * (l2 + l1)
    xnew = np.maximum(0.001, np.maximum(x - move, np.minimum(1.0, np.minimum(x + move, x * np.sqrt(-dc / lmid)))))

    if np.sum(xnew) - volfrac * nelx * nely > 0:
      l1 = lmid
    else:
      l2 = lmid

  return xnew


def check(nelx, nely, rmin, x, dc):
  """
  Mesh-independency filter
  """
  dcn = np.zeros((nely, nelx))

  for i in range(nelx):
    for j in range(nely):
      sum_val = 0.0
      for k in range(max(i - int(rmin), 0), min(i + int(rmin) + 1, nelx)):
        for l in range(max(j - int(rmin), 0), min(j + int(rmin) + 1, nely)):
          fac = rmin - np.sqrt((i - k) ** 2 + (j - l) ** 2)
          sum_val += max(0.0, fac)
          dcn[j, i] += max(0.0, fac) * x[l, k] * dc[l, k]

      dcn[j, i] /= (x[j, i] * sum_val)

  return dcn


def FE(nelx, nely, x, penal):
  """
  Finite Element Analysis
  """
  KE = lk()
  ndof = 2 * (nelx + 1) * (nely + 1)

  # Initialize global stiffness matrix
  iK = np.zeros(8 * 8 * nelx * nely, dtype=int)
  jK = np.zeros(8 * 8 * nelx * nely, dtype=int)
  sK = np.zeros(8 * 8 * nelx * nely)

  # Element stiffness matrices and assembly
  count = 0
  for elx in range(nelx):
    for ely in range(nely):
      n1 = (nely + 1) * elx + ely
      n2 = (nely + 1) * (elx + 1) + ely
      edof = np.array([2 * n1, 2 * n1 + 1, 2 * n2, 2 * n2 + 1, 2 * n2 + 2, 2 * n2 + 3, 2 * n1 + 2, 2 * n1 + 3])
      # print(edof)
      # input()
      # Insert values into the global stiffness matrix
      for i in range(8):
        for j in range(8):
          iK[count] = edof[i]
          jK[count] = edof[j]
          sK[count] = KE[i, j] * x[ely, elx] ** penal
          count += 1

  # Construct the global stiffness matrix
  K = coo_matrix((sK, (iK, jK)), shape=(ndof, ndof)).tocsc()

  # Define loads and supports (Half MBB-Beam)
  F = np.zeros(ndof)
  F[1] = -1.0  # Force at node 1 in y-direction

  # Define fixed DOFs
  fixed_dofs = np.union1d(np.arange(0, 2 * (nely + 1), 2), [2 * (nelx + 1) * (nely + 1) - 1])
  # print(fixed_dofs)
  # input()
  free_dofs = np.setdiff1d(np.arange(ndof), fixed_dofs)

  # Solve the system
  U = np.zeros(ndof)
  U[free_dofs] = spsolve(K[free_dofs, :][:, free_dofs], F[free_dofs])

  return U


def lk():
  """
  Element stiffness matrix
  """
  E = 1.0
  nu = 0.3

  k = np.array([1 / 2 - nu / 6, 1 / 8 + nu / 8, -1 / 4 - nu / 12, -1 / 8 + 3 * nu / 8,
                -1 / 4 + nu / 12, -1 / 8 - nu / 8, nu / 6, 1 / 8 - 3 * nu / 8])

  KE = E / (1 - nu ** 2) * np.array([
    [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
    [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
    [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
    [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
    [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
    [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
    [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
    [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]
  ])

  return KE


if __name__ == "__main__":
  # Default parameters
  nelx = 60  # Number of elements in x direction
  nely = 40  # Number of elements in y direction
  volfrac = 0.2  # Volume fraction
  penal = 3.0  # Penalty factor
  rmin = 1.5  # Filter radius

  # Run the optimization
  top(nelx, nely, volfrac, penal, rmin)