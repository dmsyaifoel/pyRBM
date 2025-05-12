from body import *

class PRBM:
  def __init__(self, properties):
    self.bodies = []
    self.flexures = []
    self.properties = properties
    self.dim = None

  def add_body(self, name, position=None):
    if self.dim is None: self.dim = len(position)
    if position is None: position = self.dim*[0]
    self.bodies.append(Body(name, position, self.dim))
    self.add_reference_point(name, np.zeros(self.dim))

  def add_bodies(self, names):
    assert self.dim == 2 or self.dim == 3, "Set PRBM.dim to 2 or 3 first"
    for name in names:
      self.add_body(name)

  def list_bodies(self):
    for body in self.bodies:
      print(body.name)

  def add_flexure(self, bodynameA, attachpoint_localA, bodynameB, attachpoint_localB, properties=None):

    # This first part organizes the naming of the flexures
    if bodynameA > bodynameB:
      bodynameA, bodynameB = bodynameB, bodynameA
      attachpoint_localA, attachpoint_localB = attachpoint_localB, attachpoint_localA

    name = bodynameA + bodynameB

    for body in self.bodies:
      if body.name == bodynameA: bodyA = body
      if body.name == bodynameB: bodyB = body

    occurrence = 0
    for flexure in self.flexures:
      if name in flexure.name: occurrence += 1

    name = name + str(occurrence)

    # This last part actually creates and links the flexures
    if properties is None: properties = self.properties

    flexure = Flexure(name, bodyA, attachpoint_localA, bodyB, attachpoint_localB, properties, self.dim)
    self.flexures.append(flexure)

    bodyA.flexures.append(flexure)
    bodyA.letters.append('A')
    bodyB.flexures.append(flexure)
    bodyB.letters.append('B')

  def move(self, bodyname, position, angles=None): 
    # move a body.py
    for body in self.bodies:
      if body.name == bodyname:
        if angles is None: body.move(position)
        else: body.move(position, angles)

  def moveback(self, bodyname=None): 
    # move bodies back to their original position
    for body in self.bodies:
      if body.name == bodyname or bodyname is None: body.moveback()

  def show(self, args=None): 
    # plot the prbm
    lines = [flexure.line() for flexure in self.flexures]
    names = [flexure.name for flexure in self.flexures]

    lines2 = [body.line() for body in self.bodies]
    names2 = [body.name for body in self.bodies]

    if self.dim == 2:

      if args == 'random':
        for line, name in zip(lines, names): mp.plot(line[0], line[1], 'o-', label=name, c=np.random.rand(3,))
        for line2, name2 in zip(lines2, names2): mp.plot(line2[0], line2[1], 'o-', label=name2, c=np.random.rand(3,))
      else:
        for line, name in zip(lines, names): mp.plot(line[0], line[1], 'o-', label=name)
        for line2, name2 in zip(lines2, names2): mp.plot(line2[0], line2[1], 'o-', label=name2)

      mp.xlabel('x (m)')
      mp.ylabel('y (m)')
      mp.legend()
      mp.axis('equal')
      mp.show()

    if self.dim == 3:
      fig = pg.Figure()
      for line, name in zip(lines, names): fig.add_trace(pg.Scatter3d(x=line[0], y=line[1], z=line[2], mode='lines+markers', name=name))
      for line2, name2 in zip(lines2, names2): fig.add_trace(pg.Scatter3d(x=line2[0], y=line2[1], z=line2[2], mode='lines', name=name2))

      fig.update_layout(scene=dict(aspectmode='cube'))
      fig.show()

  def energy(self): 
    # calculate and return the potential energy of the prbm
    flexure_energies = [flexure.energy() for flexure in self.flexures]
    body_energies = [body.energy() for body in self.bodies]
    return sum(flexure_energies) + sum(body_energies)

  # NOTE: added forces and torques are fixed in the global reference frame
  def add_force(self, bodyname, vector, attachpoint_local=None):
    if attachpoint_local is None: attachpoint_local = np.zeros(self.dim)
    for body in self.bodies:
      if body.name == bodyname: body.forces.append((np.array(attachpoint_local), np.array(vector)))

  def add_torque_2D(self, bodyname, torque):
    self.add_force(bodyname, (0, -torque))
    self.add_force(bodyname, (0, torque), (0, 1))

  def add_torque_x(self, bodyname, Mx):
    self.add_force(bodyname, (0, 0, -Mx))
    self.add_force(bodyname, (0, 0, Mx), (0, 1, 0))

  def add_torque_y(self, bodyname, My):
    self.add_force(bodyname, (-My, 0, 0))
    self.add_force(bodyname, (My, 0, 0), (0, 0, 1))

  def add_torque_z(self, bodyname, Mz):
    self.add_force(bodyname, (0, -Mz, 0))
    self.add_force(bodyname, (0, Mz, 0), (1, 0, 0))

  def add_torque(self, bodyname, torque):
    if self.dim == 2: self.add_torque_2D(bodyname, torque)
    else:
      self.add_torque_x(bodyname, torque[0])
      self.add_torque_y(bodyname, torque[1])
      self.add_torque_z(bodyname, torque[2])

  def clear_forces(self, bodyname=None):
    for body in self.bodies:
      if bodyname is None or body.name == bodyname: body.forces = []

  def add_reference_point(self, bodyname, attachpoint_local):
    # add points you want plotted
    for body in self.bodies:
      if body.name == bodyname: body.reference_points.append(np.array(attachpoint_local))

  def solve_pose(self, bodynames):
    # solve pose of all bodies in the list bodynames
    # using minimization of the total energy

    free_bodies = []
    x0 = []

    for body in self.bodies:
      if body.name in bodynames:
        free_bodies.append(body)
        if self.dim == 2: x0 = x0 + list(body.position) + [body.angles]

        if self.dim == 3: x0 = x0 + list(body.position) + list(body.angles)

    n_free_bodies = len(free_bodies)

    def optimize_function(x):
      for i in range(n_free_bodies):
        if self.dim == 2:
          s1 = 3*i
          s2 = 3*i + 2
          p = x[s1:s2]
          a = x[s2]
        if self.dim == 3:
          s1 = 6*i
          s2 = 6*i + 3
          s3 = 6*i + 6
          p = x[s1:s2]
          a = x[s2:s3]
        free_bodies[i].move(p, a)
      return self.energy()

    if optimizer == 'minimize': self.solution = so.minimize(optimize_function, x0)

    elif optimizer == 'shgo':
      # shgo is slower but may be better
      # bounds are needed for some optimization algorithms
      bounds = []
      for i in range(n_free_bodies):
        for j in range(self.dim): bounds.append((-.03, .03))
        for j in range(self.dim): bounds.append((-np.pi, np.pi))
      self.solution = so.shgo(optimize_function, bounds)

    else:
      print('Optimizer choice not supported')
      return

    x = self.solution.x
    optimize_function(x) # move bodies to the optimal pose

  def solve_reaction(self, end_effector, free_bodies, positions=None): 
    # solve the reaction forces using finite difference
    if positions is None:
      positions = []
      for body in self.bodies:
        if body.name == end_effector:
          positions.append(body.position)

    dim = self.dim

    if dim == 2: dim2 = 1
    else: dim2 = 3

    displacements = np.eye(dim+dim2)*finite_difference_delta

    star = []
    for position in positions:
      for displacement in displacements:
        for sign in [1, -1]: star.append((self, end_effector, free_bodies, position, sign*displacement))

    print(f"Solving for {len(star)} positions with {cpus} CPUs")
    with mu.Pool() as pool: energies = pool.starmap(poolenergy, star)

    i = 0
    self.reaction_list = []
    for position in positions:
      reactions = []
      for displacement in displacements:
        reactions.append((energies[i] - energies[i+1])/finite_difference_delta/2)
        i += 2
      self.reaction_list.append(reactions)

    if len(self.reaction_list) == 1: self.reaction = np.array(self.reaction_list[0])
