#!/usr/bin/env python3
# Python translation of the MATLAB 3D topology optimization code by Liu and Tovar (Jul 2013)
# Original MATLAB code citation:
# K. Liu and A. Tovar, "An efficient 3D topology optimization code written in Matlab",
# Struct Multidisc Optim, 50(6): 1175-1196, 2014, doi:10.1007/s00158-014-1107-x

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.sparse import coo_matrix, bmat
from scipy.sparse.linalg import spsolve
import time


def top3d(nelx, nely, nelz, volfrac, penal, rmin):
  """
  3D topology optimization function

  Parameters:
  nelx : int
      Number of elements in x direction
  nely : int
      Number of elements in y direction
  nelz : int
      Number of elements in z direction
  volfrac : float
      Volume fraction constraint (0 < volfrac < 1)
  penal : float
      Penalization power for SIMP method
  rmin : float
      Filter radius for sensitivity filtering
  """

  # USER-DEFINED LOOP PARAMETERS
  maxloop = 200  # Maximum number of iterations
  tolx = 0.01  # Termination criterion (maximum change in design variables)
  displayflag = 0  # Display structure flag (0: only final design, 1: every iteration)

  # USER-DEFINED MATERIAL PROPERTIES
  E0 = 1.0  # Young's modulus of solid material
  Emin = 1e-9  # Young's modulus of void-like material (non-zero to avoid singularity)
  nu = 0.3  # Poisson's ratio

  # PREPARE FINITE ELEMENT ANALYSIS
  # Total number of elements and degrees of freedom
  nele = nelx * nely * nelz
  ndof = 3 * (nelx + 1) * (nely + 1) * (nelz + 1)

  # USER-DEFINED LOAD DOFs
  # Create a meshgrid for the load nodes (at the top face)
  # Using 0-based indexing for Python
  il, jl, kl = np.meshgrid(range(nelx + 1), [0], range(nelz + 1))
  # Calculate node IDs for the load application nodes
  loadnid = kl * (nelx + 1) * (nely + 1) + il * (nely + 1) + (nely - jl)
  # Define the DOFs where loads are applied (y-direction, hence the -1)
  loaddof = 3 * loadnid.flatten() - 1

  # USER-DEFINED SUPPORT FIXED DOFs (fixed at x=0 face)
  # Create a meshgrid for the support nodes
  iif, jf, kf = np.meshgrid([0], range(nely + 1), range(nelz + 1))
  # Calculate node IDs for the fixed support nodes
  fixednid = kf * (nelx + 1) * (nely + 1) + iif * (nely + 1) + (nely - jf)
  # Define all DOFs that are fixed (all directions, hence -2, -1, 0)
  fixeddof = np.concatenate([
    3 * fixednid.flatten(),  # z-direction
    3 * fixednid.flatten() - 1,  # y-direction
    3 * fixednid.flatten() - 2  # x-direction
  ])

  # PREPARE FORCE VECTOR AND INITIALIZE DISPLACEMENT VECTOR
  # Create sparse force vector with -1 at each load DOF (negative y-direction)
  F = coo_matrix(([-1.0] * len(loaddof), (loaddof, [0] * len(loaddof))), shape=(ndof, 1)).tocsr()
  # Initialize displacement vector
  U = np.zeros((ndof, 1))

  # Define free DOFs (all DOFs except fixed ones)
  freedofs = np.setdiff1d(np.arange(1, ndof + 1), fixeddof)

  # GET ELEMENT STIFFNESS MATRIX
  KE = lk_H8(nu)

  # PREPARE FINITE ELEMENT ANALYSIS SYSTEM MATRICES
  # Create node grid for element connectivity
  nodegrd = np.reshape(np.arange(1, (nely + 1) * (nelx + 1) + 1), (nely + 1, nelx + 1))
  # Get nodeids for the first 2D layer
  nodeids = np.reshape(nodegrd[0:nely, 0:nelx], (nely * nelx, 1))
  # Get z-layer offsets
  nodeidz = np.arange(0, (nely + 1) * (nelx + 1) * nelz, (nely + 1) * (nelx + 1))
  # Create node IDs for all elements by repeating 2D pattern for each z-layer
  nodeids = np.tile(nodeids, (len(nodeidz), 1)) + np.repeat(nodeidz, nely * nelx).reshape(-1, 1)

  # Calculate DOF indices for each element
  edofVec = 3 * nodeids.flatten() + 1
  edofMat = np.tile(edofVec, (1, 24)) + \
            np.tile(np.array([0, 1, 2, 3 * nely + 3, 3 * nely + 4, 3 * nely + 5,
                              3 * nely + 0, 3 * nely + 1, 3 * nely + 2, -3, -2, -1,
                              3 * (nely + 1) * (nelx + 1) + np.array([0, 1, 2, 3 * nely + 3, 3 * nely + 4,
                                                                      3 * nely + 5, 3 * nely + 0, 3 * nely + 1,
                                                                      3 * nely + 2, -3, -2, -1])]), (nele, 1))

  # Create indices for global stiffness matrix assembly
  iK = np.kron(edofMat, np.ones((24, 1))).flatten()
  jK = np.kron(edofMat, np.ones((1, 24))).flatten()

  # PREPARE FILTER
  # Initialize filter matrices
  # Calculate maximum number of neighboring elements within filter radius
  max_neighbors = (2 * int(np.ceil(rmin)) + 1) ** 3
  iH = np.zeros(nele * max_neighbors, dtype=int)
  jH = np.zeros(nele * max_neighbors, dtype=int)
  sH = np.zeros(nele * max_neighbors)

  k = 0
  # Loop over all elements to compute filter weights
  for k1 in range(1, nelz + 1):
    for i1 in range(1, nelx + 1):
      for j1 in range(1, nely + 1):
        # Calculate current element index (0-based)
        e1 = (k1 - 1) * nelx * nely + (i1 - 1) * nely + j1 - 1
        # Loop over neighboring elements within filter radius
        for k2 in range(max(k1 - int(np.ceil(rmin)), 1), min(k1 + int(np.ceil(rmin)), nelz + 1)):
          for i2 in range(max(i1 - int(np.ceil(rmin)), 1), min(i1 + int(np.ceil(rmin)), nelx + 1)):
            for j2 in range(max(j1 - int(np.ceil(rmin)), 1), min(j1 + int(np.ceil(rmin)), nely + 1)):
              # Calculate neighbor element index (0-based)
              e2 = (k2 - 1) * nelx * nely + (i2 - 1) * nely + j2 - 1
              # Calculate distance between elements
              dist = np.sqrt((i1 - i2) ** 2 + (j1 - j2) ** 2 + (k1 - k2) ** 2)
              # Only include elements within filter radius
              if dist <= rmin:
                # Store filter data
                iH[k] = e1  # Current element index
                jH[k] = e2  # Neighbor element index
                sH[k] = max(0, rmin - dist)  # Weight based on distance
                k += 1

  # Trim arrays to only used entries
  iH = iH[:k]
  jH = jH[:k]
  sH = sH[:k]

  # Create sparse filter matrix
  H = coo_matrix((sH, (iH, jH)), shape=(nele, nele)).tocsr()
  Hs = np.array(H.sum(axis=1)).flatten()

  # INITIALIZE ITERATION
  # Initialize design variables with uniform distribution = volume fraction
  x = np.ones((nely, nelx, nelz)) * volfrac
  xPhys = x.copy()  # Physical density field (after filtering)
  loop = 0
  change = 1.0

  # START ITERATION
  print("Starting optimization loop...")
  print("Iter.   Objective    Volume    Change")
  print("-" * 40)

  # Main optimization loop
  while change > tolx and loop < maxloop:
    loop += 1

    # FE-ANALYSIS
    # Compute element stiffness matrices based on current density
    sK = np.zeros(24 * 24 * nele)

    # Create a vector with all element stiffness values based on density
    xPhys_penal = xPhys.flatten() ** penal
    for i in range(24 * 24):
      # Apply SIMP (Solid Isotropic Material with Penalization)
      # Stiffness is a function of physical density
      sK[i::24 * 24] = KE.flatten()[i] * (Emin + xPhys_penal * (E0 - Emin))

    # Assemble global stiffness matrix (convert from 1-indexed to 0-indexed)
    K = coo_matrix((sK, (iK.astype(int), jK.astype(int))), shape=(ndof, ndof)).tocsr()
    K = (K + K.T) / 2  # Ensure symmetry

    # Solve FE system KU=F for displacements
    # Note: we use a direct 0-indexing now
    U_free = spsolve(K[freedofs - 1, :][:, freedofs - 1], F[freedofs - 1])
    U[freedofs - 1] = U_free.reshape(-1, 1)

    # OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    # Calculate strain energy (compliance) for each element
    ce = np.zeros((nely, nelx, nelz))
    for i in range(nele):
      # Get element coordinates
      ez = i // (nelx * nely)
      exy = i % (nelx * nely)
      ex = exy // nely
      ey = exy % nely

      # Extract element DOFs
      edof = edofMat[i, :]
      # Get element displacement vector (accounting for 0-indexing)
      Ue = U[edof - 1].flatten()
      # Calculate element strain energy: Ue^T * KE * Ue
      ce[ey, ex, ez] = np.dot(Ue, np.dot(KE, Ue))

    # Calculate total compliance (objective function)
    c = np.sum((Emin + xPhys ** penal * (E0 - Emin)) * ce)

    # Sensitivity analysis
    dc = -penal * (E0 - Emin) * xPhys ** (penal - 1) * ce  # Derivative of compliance w.r.t. density
    dv = np.ones((nely, nelx, nelz))  # Derivative of volume w.r.t. density

    # FILTERING AND MODIFICATION OF SENSITIVITIES
    dc_vec = dc.flatten()
    dv_vec = dv.flatten()

    # Apply filter to sensitivity
    dc_vec = H.dot(dc_vec / Hs)
    dv_vec = H.dot(dv_vec / Hs)

    # Reshape back to 3D arrays
    dc = dc_vec.reshape((nely, nelx, nelz))
    dv = dv_vec.reshape((nely, nelx, nelz))

    # OPTIMALITY CRITERIA UPDATE
    l1 = 0
    l2 = 1e9
    move = 0.2  # Limit for design change in each iteration

    # Bisection algorithm to find Lagrange multiplier
    while (l2 - l1) / (l1 + l2) > 1e-3:
      lmid = 0.5 * (l2 + l1)

      # Calculate new design variables based on optimality criteria
      # sqrt(-dc/dv) is the optimal condition for density
      xnew = np.maximum(0.0, np.maximum(x - move,
                                        np.minimum(1.0, np.minimum(x + move,
                                                                   x * np.sqrt(-dc / dv / lmid)))))

      # Apply density filter to get physical densities
      xPhys_vec = H.dot(xnew.flatten()) / Hs
      xPhys_new = xPhys_vec.reshape((nely, nelx, nelz))

      # Update Lagrange multiplier based on volume constraint
      if np.mean(xPhys_new) > volfrac:
        l1 = lmid
      else:
        l2 = lmid

    # Update physical densities
    xPhys = xPhys_new.copy()

    # Calculate maximum change in design variables
    change = np.max(np.abs(xnew - x))
    x = xnew.copy()

    # Print results for this iteration
    print(f"{loop:5d}   {c:11.4f}   {np.mean(xPhys):7.3f}   {change:7.3f}")

    # PLOT DENSITIES
    if displayflag:
      display_3D(xPhys)
      plt.title(f"Iteration {loop}")
      plt.draw()
      plt.pause(0.01)

  # PLOT FINAL DESIGN
  print("\nOptimization completed!")
  print(f"Final compliance: {c:.4f}")
  print(f"Final volume fraction: {np.mean(xPhys):.3f}")

  display_3D(xPhys)
  plt.title("Final design")
  plt.show()

  return xPhys


def lk_H8(nu):
  """
  Generate the element stiffness matrix for a hexahedral element (H8)

  Parameters:
  nu : float
      Poisson's ratio

  Returns:
  KE : numpy.ndarray
      24x24 element stiffness matrix
  """
  # Define coefficient matrix A for material property calculations
  A = np.array([
    [32, 6, -8, 6, -6, 4, 3, -6, -10, 3, -3, -3, -4, -8],
    [-48, 0, 0, -24, 24, 0, 0, 0, 12, -12, 0, 12, 12, 12]
  ])

  # Calculate material property coefficients
  k = 1 / 144 * A.T @ np.array([1, nu])

  # Define the submatrices of the element stiffness matrix
  K1 = np.array([
    [k[0], k[1], k[1], k[2], k[4], k[4]],
    [k[1], k[0], k[1], k[3], k[5], k[6]],
    [k[1], k[1], k[0], k[3], k[6], k[5]],
    [k[2], k[3], k[3], k[0], k[7], k[7]],
    [k[4], k[5], k[6], k[7], k[0], k[1]],
    [k[4], k[6], k[5], k[7], k[1], k[0]]
  ])

  K2 = np.array([
    [k[8], k[7], k[11], k[5], k[3], k[6]],
    [k[7], k[8], k[11], k[4], k[2], k[4]],
    [k[9], k[9], k[12], k[6], k[3], k[5]],
    [k[5], k[4], k[10], k[8], k[1], k[9]],
    [k[3], k[2], k[4], k[1], k[8], k[11]],
    [k[10], k[3], k[5], k[11], k[9], k[12]]
  ])

  K3 = np.array([
    [k[5], k[6], k[3], k[8], k[11], k[7]],
    [k[6], k[5], k[3], k[9], k[12], k[9]],
    [k[4], k[4], k[2], k[7], k[11], k[8]],
    [k[8], k[9], k[1], k[5], k[10], k[4]],
    [k[11], k[12], k[9], k[10], k[5], k[3]],
    [k[1], k[11], k[8], k[3], k[4], k[2]]
  ])

  K4 = np.array([
    [k[13], k[10], k[10], k[12], k[9], k[9]],
    [k[10], k[13], k[10], k[11], k[8], k[7]],
    [k[10], k[10], k[13], k[11], k[7], k[8]],
    [k[12], k[11], k[11], k[13], k[6], k[6]],
    [k[9], k[8], k[7], k[6], k[13], k[10]],
    [k[9], k[7], k[8], k[6], k[10], k[13]]
  ])

  K5 = np.array([
    [k[0], k[1], k[7], k[2], k[4], k[3]],
    [k[1], k[0], k[7], k[3], k[5], k[10]],
    [k[7], k[7], k[0], k[4], k[10], k[5]],
    [k[2], k[3], k[4], k[0], k[7], k[1]],
    [k[4], k[5], k[10], k[7], k[0], k[7]],
    [k[3], k[10], k[5], k[1], k[7], k[0]]
  ])

  K6 = np.array([
    [k[13], k[10], k[6], k[12], k[9], k[11]],
    [k[10], k[13], k[6], k[11], k[8], k[1]],
    [k[6], k[6], k[13], k[9], k[1], k[8]],
    [k[12], k[11], k[9], k[13], k[6], k[10]],
    [k[9], k[8], k[1], k[6], k[13], k[6]],
    [k[11], k[1], k[8], k[10], k[6], k[13]]
  ])

  # Assemble the full 24x24 element stiffness matrix
  KE = 1 / ((nu + 1) * (1 - 2 * nu)) * np.block([
    [K1, K2, K3, K4],
    [K2.T, K5, K6, K3.T],
    [K3.T, K6, K5.T, K2.T],
    [K4, K3, K2, K1]
  ])

  return KE


def display_3D(rho):
  """
  Display 3D topology optimization result (ISO-VIEW)

  Parameters:
  rho : numpy.ndarray
      3D array of densities
  """
  # Get dimensions of the design
  nely, nelx, nelz = rho.shape

  # User-defined unit element size
  hx, hy, hz = 1, 1, 1

  # Define element face connectivity
  face = np.array([
    [0, 1, 2, 3],  # Bottom face
    [1, 5, 6, 2],  # Side face
    [3, 2, 6, 7],  # Side face
    [0, 4, 7, 3],  # Side face
    [0, 1, 5, 4],  # Side face
    [4, 5, 6, 7]  # Top face
  ])

  # Create new figure
  fig = plt.figure(figsize=(10, 8))
  ax = fig.add_subplot(111, projection='3d')

  # Loop through all elements
  for k in range(nelz):
    z = k * hz
    for i in range(nelx):
      x = i * hx
      for j in range(nely):
        y = nely * hy - j * hy

        # Only display elements with density above threshold
        if rho[j, i, k] > 0.5:  # User-defined display density threshold
          # Define vertices of the element
          vert = np.array([
            [x, y, z],
            [x, y - hy, z],
            [x + hx, y - hy, z],
            [x + hx, y, z],
            [x, y, z + hz],
            [x, y - hy, z + hz],
            [x + hx, y - hy, z + hz],
            [x + hx, y, z + hz]
          ])

          # Swap y and z coordinates for ISO view
          vert[:, [1, 2]] = vert[:, [2, 1]]
          vert[:, 1] = -vert[:, 1]

          # Create list of faces
          faces = []
          for f in face:
            faces.append([vert[f[0]], vert[f[1]], vert[f[2]], vert[f[3]]])

          # Create collection of polygons
          poly = Poly3DCollection(faces)
          gray_level = 0.2 + 0.8 * (1 - rho[j, i, k])
          poly.set_facecolor((gray_level, gray_level, gray_level))
          poly.set_edgecolor('black')
          poly.set_linewidth(0.1)
          ax.add_collection3d(poly)

  # Set axis properties
  ax.set_xlim([0, nelx])
  ax.set_ylim([-nelz, 0])
  ax.set_zlim([0, nely])
  ax.set_xlabel('X')
  ax.set_ylabel('Z')
  ax.set_zlabel('Y')
  ax.set_box_aspect([nelx, nelz, nely])  # Equal aspect ratio
  ax.view_init(elev=30, azim=30)

  # Remove axis ticks
  ax.set_xticks([])
  ax.set_yticks([])
  ax.set_zticks([])

  # Update the plot
  plt.tight_layout()
  plt.draw()
  plt.pause(0.001)

  return fig


# Run the program as a standalone script
if __name__ == "__main__":
  # Example usage
  print("3D Topology Optimization Example")
  print("--------------------------------")

  # Get user inputs or use defaults
  nelx = int(input("Enter number of elements in x-direction (default 30): ") or 30)
  nely = int(input("Enter number of elements in y-direction (default 10): ") or 10)
  nelz = int(input("Enter number of elements in z-direction (default 10): ") or 10)
  volfrac = float(input("Enter volume fraction (default 0.3): ") or 0.3)
  penal = float(input("Enter penalization parameter (default 3.0): ") or 3.0)
  rmin = float(input("Enter filter radius (default 1.5): ") or 1.5)

  # Run the optimization
  start_time = time.time()
  result = top3d(nelx, nely, nelz, volfrac, penal, rmin)
  end_time = time.time()

  print(f"\nOptimization completed in {end_time - start_time:.2f} seconds")
  print("\nPress any key to exit...")
  plt.show()  # Keep the plot window open until closed