import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def top_compliant(nelx, nely, volfrac, penal, rmin, input_point, output_point):
  """
  Topology optimization for compliant mechanisms

  Parameters:
  nelx - Number of elements in x direction
  nely - Number of elements in y direction
  volfrac - Volume fraction
  penal - Penalty factor
  rmin - Filter radius
  input_point - Tuple (x, y, direction) for input force (direction: 0=x, 1=y)
  output_point - Tuple (x, y, direction) for output displacement (direction: 0=x, 1=y)
  """
  # INITIALIZE
  x = volfrac * np.ones((nely, nelx))
  loop = 0
  change = 1.0

  # PREPARE FINITE ELEMENT ANALYSIS
  ndof = 2 * (nelx + 1) * (nely + 1)
  KE = lk()

  # INPUT AND OUTPUT LOCATIONS
  input_x, input_y, input_dir = input_point
  output_x, output_y, output_dir = output_point

  # Convert to node numbers and DOFs
  input_node = (nely + 1) * input_x + input_y
  input_dof = 2 * input_node + input_dir

  output_node = (nely + 1) * output_x + output_y
  output_dof = 2 * output_node + output_dir

  # PREPARE FORCE VECTORS
  # Input force vector (for calculating displacement due to input force)
  F_in = np.zeros(ndof)
  F_in[input_dof] = 1.0

  # Output "dummy" force vector (for calculating mutual energy)
  F_out = np.zeros(ndof)
  F_out[output_dof] = -1.0  # Negative for displacement maximization

  # START ITERATION
  max_iter = 150  # Maximum number of iterations
  plt.figure(figsize=(10, 5))

  # Keep track of objective history for plotting
  obj_history = []

  while change > 0.01 and loop < max_iter:
    loop += 1
    xold = x.copy()

    # FE-ANALYSIS
    U_in = FE(nelx, nely, x, penal, F_in)  # Displacement from input force
    U_out = FE(nelx, nely, x, penal, F_out)  # Displacement from "dummy" force

    # OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    # For compliant mechanisms, we want to maximize output displacement
    # in the direction of interest due to input force
    obj = U_in[output_dof]  # Output displacement due to input force
    obj_history.append(obj)

    # Calculate compliance for scaling sensitivity
    compliance = 0.0
    dc = np.zeros((nely, nelx))
    dv = np.ones((nely, nelx))

    # Element-wise sensitivity calculation
    for elx in range(nelx):
      for ely in range(nely):
        n1 = (nely + 1) * elx + ely
        n2 = (nely + 1) * (elx + 1) + ely
        edof = np.array([2 * n1, 2 * n1 + 1, 2 * n2, 2 * n2 + 1, 2 * n2 + 2, 2 * n2 + 3, 2 * n1 + 2, 2 * n1 + 3])

        Ue_in = U_in[edof].flatten()
        Ue_out = U_out[edof].flatten()

        # Sensitivity for compliant mechanism
        # Use mutual potential energy (MPE) as sensitivity
        dc[ely, elx] = -penal * x[ely, elx] ** (penal - 1) * Ue_in.dot(KE.dot(Ue_out))
        compliance += x[ely, elx] ** penal * Ue_in.dot(KE.dot(Ue_in))

    # FILTERING OF SENSITIVITIES
    dc = check(nelx, nely, rmin, x, dc)

    # DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    x = OC(nelx, nely, x, volfrac, dc)

    # PRINT RESULTS
    change = np.max(np.abs(x - xold))
    print(f" It.: {loop:4d} Obj.: {obj:10.6f} Vol.: {np.sum(x) / (nelx * nely):6.3f} ch.: {change:6.3f}")

    # PLOT DENSITIES
    if loop % 5 == 0 or loop == 1:  # Plot every 5 iterations to speed up
      plt.clf()
      plt.imshow(-x, cmap='gray', interpolation='none')
      plt.axis('equal')
      plt.axis('off')

      # Mark input and output points
      plt.plot(input_x, input_y, 'ro', markersize=10)  # Input force
      plt.plot(output_x, output_y, 'bo', markersize=10)  # Output location

      plt.title(f"Iteration {loop}, Objective: {obj:.6f}")
      plt.pause(0.01)

  # FINAL PLOT WITH ANNOTATIONS
  plt.figure(figsize=(12, 8))

  # Plot final design
  plt.subplot(2, 1, 1)
  plt.imshow(-x, cmap='gray', interpolation='none')
  plt.axis('equal')
  plt.axis('tight')

  # Annotate input and output
  dir_labels = ['x', 'y']
  plt.plot(input_x, input_y, 'ro', markersize=10)
  plt.plot(output_x, output_y, 'bo', markersize=10)

  plt.annotate(f"Input ({dir_labels[input_dir]})",
               xy=(input_x, input_y), xytext=(input_x - 10, input_y - 10),
               arrowprops=dict(facecolor='red', shrink=0.05))
  plt.annotate(f"Output ({dir_labels[output_dir]})",
               xy=(output_x, output_y), xytext=(output_x + 10, output_y + 10),
               arrowprops=dict(facecolor='blue', shrink=0.05))

  plt.title(f"Compliant Mechanism Design (Output Displacement: {obj:.6f})")

  # Plot objective history
  plt.subplot(2, 1, 2)
  plt.plot(range(1, len(obj_history) + 1), obj_history)
  plt.xlabel('Iteration')
  plt.ylabel('Output Displacement')
  plt.title('Convergence History')
  plt.grid(True)

  plt.tight_layout()
  plt.show()

  return x, obj


def OC(nelx, nely, x, volfrac, dc):
  """
  Optimality criteria update
  """
  l1 = 0
  l2 = 1e5
  move = 0.2

  # Normalize sensitivities
  dc = dc / np.max(np.abs(dc))

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

      if sum_val > 0:
        dcn[j, i] /= (x[j, i] * sum_val)

  return dcn


def FE(nelx, nely, x, penal, f):
  """
  Finite Element Analysis

  Parameters:
  nelx, nely - Element dimensions
  x - Design variable field
  penal - Penalty factor
  f - Force vector
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

      # Add lower bound on stiffness for stability
      # Using SIMP method (Solid Isotropic Material with Penalization)
      xe = max(x[ely, elx], 0.001) ** penal

      # Insert values into the global stiffness matrix
      for i in range(8):
        for j in range(8):
          iK[count] = edof[i]
          jK[count] = edof[j]
          sK[count] = KE[i, j] * xe
          count += 1

  # Construct the global stiffness matrix
  K = coo_matrix((sK, (iK, jK)), shape=(ndof, ndof)).tocsc()

  # Define fixed DOFs (boundary conditions)
  # Fix the left edge nodes (all DOFs)
  fixed_dofs = np.arange(0, 2 * (nely + 1), 1)
  free_dofs = np.setdiff1d(np.arange(ndof), fixed_dofs)

  # Solve the system
  # Add small value to diagonal for stability
  for i in range(K.shape[0]):
    K[i, i] += 1e-9

  u = np.zeros(ndof)
  u[free_dofs] = spsolve(K[free_dofs, :][:, free_dofs], f[free_dofs])

  return u


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


def create_inverter(nelx=60, nely=30):
  """
  Creates a compliant inverter mechanism
  """
  volfrac = 0.3
  penal = 3.0
  rmin = 1.5

  # Define input and output points
  # (x, y, direction) where direction is 0 for x and 1 for y
  input_point = (0, nely // 2, 0)  # Middle of left edge, x-direction
  output_point = (nelx, nely // 2, 0)  # Middle of right edge, x-direction

  print("Designing an inverter mechanism...")
  return top_compliant(nelx, nely, volfrac, penal, rmin, input_point, output_point)


def create_gripper(nelx=60, nely=30):
  """
  Creates a compliant gripper mechanism
  """
  volfrac = 0.3
  penal = 3.0
  rmin = 1.5

  # Define input and output points
  # (x, y, direction) where direction is 0 for x and 1 for y
  input_point = (0, nely // 2, 0)  # Middle of left edge, x-direction
  output_point = (nelx - 1, nely // 4, 1)  # Upper quarter of right edge, y-direction

  print("Designing a gripper mechanism...")
  return top_compliant(nelx, nely, volfrac, penal, rmin, input_point, output_point)


def create_amplifier(nelx=60, nely=30):
  """
  Creates a compliant displacement amplifier
  """
  volfrac = 0.3
  penal = 3.0
  rmin = 1.5

  # Define input and output points
  input_point = (0, nely // 2, 0)  # Middle of left edge, x-direction
  output_point = (nelx, nely // 2, 0)  # Middle of right edge, x-direction

  print("Designing a displacement amplifier...")
  return top_compliant(nelx, nely, volfrac, penal, rmin, input_point, output_point)


def create_custom_mechanism(nelx=60, nely=30):
  """
  Creates a custom compliant mechanism based on user inputs
  """
  volfrac = float(input("Enter volume fraction (0.1-0.5): "))
  penal = float(input("Enter penalization factor (default 3.0): ") or "3.0")
  rmin = float(input("Enter filter radius (default 1.5): ") or "1.5")

  # Get input location
  print("\nInput force location:")
  input_x = int(input(f"X coordinate (0-{nelx}): "))
  input_y = int(input(f"Y coordinate (0-{nely}): "))
  input_dir = int(input("Direction (0 for x, 1 for y): "))

  # Get output location
  print("\nOutput displacement location:")
  output_x = int(input(f"X coordinate (0-{nelx}): "))
  output_y = int(input(f"Y coordinate (0-{nely}): "))
  output_dir = int(input("Direction (0 for x, 1 for y): "))

  input_point = (input_x, input_y, input_dir)
  output_point = (output_x, output_y, output_dir)

  print("\nDesigning your custom mechanism...")
  return top_compliant(nelx, nely, volfrac, penal, rmin, input_point, output_point)


if __name__ == "__main__":
  print("Compliant Mechanism Design Tool")
  print("===============================")
  print("1. Inverter mechanism")
  print("2. Gripper mechanism")
  print("3. Displacement amplifier")
  print("4. Custom mechanism")

  choice = input("\nSelect mechanism type (1-4): ")

  if choice == '1':
    create_inverter()
  elif choice == '2':
    create_gripper()
  elif choice == '3':
    create_amplifier()
  elif choice == '4':
    create_custom_mechanism()
  else:
    print("Invalid choice. Running inverter mechanism as default.")
    create_inverter()