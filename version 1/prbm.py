from body import *

class PRBM:
  def __init__(self, properties):
    self.bodies = []
    self.flexures = []
    self.properties = properties
    self.dim = None
    self.optimizer = optimizer
    self.delta = finite_difference_delta

  def add_body(self, name, position):

    if self.dim is None:
      self.dim = len(position)

    self.bodies.append(Body(name, position, self.dim))
    self.add_reference_point(name, np.zeros(self.dim))

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
    if properties is None:
      properties = self.properties

    flexure = Flexure(name, bodyA, attachpoint_localA, bodyB, attachpoint_localB, properties, self.dim)
    self.flexures.append(flexure)

    bodyA.flexures.append(flexure)
    bodyA.letters.append('A')
    bodyB.flexures.append(flexure)
    bodyB.letters.append('B')

  def move(self, bodyname, position, angles=None):
    for body in self.bodies:
      if body.name == bodyname:
        if angles is None:
          body.move(position)
        else:
          body.move(position, angles)

  def moveback(self, bodyname=None):
    for body in self.bodies:
      if body.name == bodyname or bodyname is None:
        body.moveback()

  def show(self):
    lines = [flexure.line() for flexure in self.flexures]
    names = [flexure.name for flexure in self.flexures]

    lines2 = [body.line() for body in self.bodies]
    names2 = [body.name for body in self.bodies]

    if self.dim == 2:
      for line, name in zip(lines, names):
        mp.plot(line[0], line[1], 'o-', label=name)

      for line2, name2 in zip(lines2, names2):
        mp.plot(line2[0], line2[1], 'o-', label=name2)

      mp.xlabel('x (m)')
      mp.ylabel('y (m)')
      mp.legend()
      mp.show()

    if self.dim == 3:
      fig = pg.Figure()
      for line, name in zip(lines, names):
        fig.add_trace(pg.Scatter3d(x=line[0], y=line[1], z=line[2],
                                    mode='lines+markers', name=name))

      for line2, name2 in zip(lines2, names2):
        fig.add_trace(pg.Scatter3d(x=line2[0], y=line2[1], z=line2[2],
                                    mode='lines', name=name2))


      fig.update_layout(scene=dict(aspectmode='cube'))
      fig.show()

  def energy(self):
    flexure_energies = [flexure.energy() for flexure in self.flexures]
    body_energies = [body.energy() for body in self.bodies]
    return sum(flexure_energies) + sum(body_energies)

  # added forces and torques are fixed in the global reference frame
  def add_force(self, bodyname, vector, attachpoint_local=None):
    if attachpoint_local is None:
      attachpoint_local = np.zeros(self.dim)
    for body in self.bodies:
      if body.name == bodyname:
        body.forces.append((np.array(attachpoint_local), np.array(vector)))

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
    if self.dim == 2:
      self.add_torque_2D(bodyname, torque)
    else:
      self.add_torque_x(bodyname, torque[0])
      self.add_torque_y(bodyname, torque[1])
      self.add_torque_z(bodyname, torque[2])

  def clear_forces(self, bodyname=None):
    for body in self.bodies:
      if bodyname is None or body.name == bodyname:
        body.forces = []

  def add_reference_point(self, bodyname, attachpoint_local):
    for body in self.bodies:
      if body.name == bodyname:
        body.reference_points.append(np.array(attachpoint_local))

  def solve_pose(self, bodynames):
    free_bodies = []
    x0 = []

    for body in self.bodies:
      if body.name in bodynames:
        free_bodies.append(body)
        if self.dim == 2:
          x0 = x0 + list(body.position) + [body.angles]

        if self.dim == 3:
          x0 = x0 + list(body.position) + list(body.angles)

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

    if self.optimizer == 'minimize':
      self.solution = so.minimize(optimize_function, x0)

    elif self.optimizer == 'shgo':
      # shgo is slower but may be better
      # bounds are needed for some optimization algorithms
      bounds = []
      for i in range(n_free_bodies):
        for j in range(self.dim):
          bounds.append((-.03, .03))
        for j in range(self.dim):
          bounds.append((-np.pi, np.pi))

      self.solution = so.shgo(optimize_function, bounds)

    x = self.solution.x
    optimize_function(x)

  def solve_reaction(self, end_effector, free_bodies):
    for body in self.bodies:
      if body.name == end_effector:
        initial_position = body.position
        initial_angles = body.angles

    energies = [[], []]

    displacements = np.eye(self.dim)*self.delta

    for i in range(2):
      mult = i*2 - 1
      for displacement in displacements:
          self.move(end_effector, initial_position + mult*displacement)
          if len(free_bodies) > 0:
            self.solve_pose(free_bodies)
          energies[i].append(self.energy())

      if self.dim == 3:
        for displacement in displacements:
            self.move(end_effector, initial_position, initial_angles + mult*displacement)
            if len(free_bodies) > 0:
              self.solve_pose(free_bodies)
            energies[i].append(self.energy())
      else:
        self.move(end_effector, initial_position, initial_angles + self.delta)
        if len(free_bodies) > 0:
          self.solve_pose(free_bodies)
        energies[i].append(self.energy())

    self.reaction = (np.array(energies[1]) - np.array(energies[0]))/2/self.delta
    print(self.reaction)

    for body in self.bodies:
      if body.name == end_effector:
        body.move(initial_position, initial_angles)

    self.solve_pose(free_bodies)
