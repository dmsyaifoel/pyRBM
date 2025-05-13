from body import Body
from flexure import Flexure
from functions import zeros
from backend import array, minimize

class PRBM:
  '''
  Contains useful functions for setting up a pseudo rigid body model.
  These functions are just for convenience: it might be easier to set up something manual.
  '''
  def __init__(self, dim):
    self.dim = dim
    self.bodies = {}
    self.bodynames = []
    self.flexures = {}
    self.flexurenames = []

  def add_body(self, name, position=None):
    self.bodies[name] = Body(name, position, self.dim)
    self.bodynames.append(name)

  def add_flexure(self, bodynameA, attachpoint_localA, bodynameB, attachpoint_localB):
    bodyA = self.bodies[bodynameA]
    bodyB = self.bodies[bodynameB]

    name = bodynameA + bodynameB
    occurrence = 0
    for othername in self.flexurenames:
      if name in othername: occurrence += 1
    name = name + str(occurrence)
    self.flexurenames.append(name)

    flexure = Flexure(bodyA, attachpoint_localA, bodyB, attachpoint_localB)
    self.flexures[name] = flexure

    bodyA.flexures.append(flexure)
    bodyA.which.append(True)
    bodyA.points.append(attachpoint_localA)

    bodyB.flexures.append(flexure)
    bodyB.which.append(False)
    bodyB.points.append(attachpoint_localB)

  def move(self, bodyname, position, angles=None):
    # move a body
      if angles is None:
          self.bodies[bodyname].move(position)
      else:
          self.bodies[bodyname].move(position, angles)

  def show(self, args=None):
    # plot the prbm

    import numpy as np

    lines = [np.array([flexure.attachpoint_globalA,
                       flexure.springpoint_globalA,
                       flexure.springpoint_globalB,
                       flexure.attachpoint_globalB]).T for flexure in self.flexures.values()]

    lines2 = [np.array(body.line()).T for body in self.bodies.values()]

    if self.dim == 2:
      import matplotlib.pyplot as mp

      for line2, name2 in zip(lines2, self.bodynames):
        mp.plot(line2[0], line2[1], '-', label=name2)

      for line, name in zip(lines, self.flexurenames):
        mp.plot(line[0], line[1], 'o-', label=name)

      mp.xlabel('x (m)')
      mp.ylabel('y (m)')
      mp.legend()
      mp.axis('equal')
      mp.show()

    if self.dim == 3:
      import plotly.graph_objects as pg
      fig = pg.Figure()

      for line2, name2 in zip(lines2, self.bodynames):
        fig.add_trace(pg.Scatter3d(x=line2[0], y=line2[1], z=line2[2], mode='lines', name=name2))

      for line, name in zip(lines, self.flexurenames):
        fig.add_trace(pg.Scatter3d(x=line[0], y=line[1], z=line[2], mode='lines+markers', name=name))

      fig.update_layout(scene=dict(aspectmode='cube'))
      fig.show()

  def energy(self, A, E, I):
    # calculate and return the potential energy of the prbm
    flexure_energies = [flexure.energy(A, E, I) for flexure in self.flexures.values()]
    body_energies = [body.energy() for body in self.bodies.values()]
    return sum(flexure_energies) + sum(body_energies)

  # NOTE: added forces and torques are fixed in the global reference frame
  def add_force(self, bodyname, vector, attachpoint_local=None):
    body = self.bodies[bodyname]
    if attachpoint_local is None:
      attachpoint_local = zeros(self.dim)
    else:
      body.points.append(attachpoint_local)
    body.forces.append((attachpoint_local, vector))



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

  def solve_pose(self, bodynames, A, E, I):
    # solve pose of all bodies in the list bodynames
    # using minimization of the total energy

    free_bodies = []
    x0 = []

    for bodyname in bodynames:
      body = self.bodies[bodyname]
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
      return self.energy(A, E, I)

    self.solution = minimize(optimize_function, x0)

    x = self.solution.x
    optimize_function(x) # move bodies to the optimal pose