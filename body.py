from functions import rotmat
from numpy import array, zeros, dot

class Body:
  '''
  Class containing data and functions of a PRBM rigidbody

  Contains:
  - data, updated by self.move()
  - move(), updates the attached flexures
  - energy(), calculates and returns energy due to external forces
  '''
  def __init__(self, name, position0, dim):

    self.name = name
    self.position0 = array(position0)
    self.position = array(position0)
    self.dim = dim

    if self.dim == 2:
      self.angles0 = 0
      self.angles = 0

    if self.dim == 3:
      self.angles0 = zeros(3)
      self.angles = zeros(3)

    self.rotmat = rotmat(self.angles, self.dim)

    self.flexures = []
    self.which = []

    self.forces = []
    self.points = [zeros(dim)]

  def move(self, position, angles=None):
    # move a body to a pose
    self.position = array(position)
    if angles is None: self.angles = self.angles0
    else: self.angles = angles

    self.rotmat = rotmat(self.angles, self.dim)

    # then update all flexure attached
    for flexure, which in zip(self.flexures, self.which):
      if which:
        flexure.attachpoint_globalA = self.rotmat@flexure.attachpoint_localA + self.position
        flexure.springpoint_globalA = self.rotmat@flexure.springpoint_localA + self.position
        flexure.unitvector_globalA = self.rotmat@flexure.unitvector_localA

      else:
        flexure.attachpoint_globalB = self.rotmat@flexure.attachpoint_localB + self.position
        flexure.springpoint_globalB = self.rotmat@flexure.springpoint_localB + self.position
        flexure.unitvector_globalB = self.rotmat@flexure.unitvector_localB

  def energy(self):
    # calculate the potential energy of a body due to external forces
    energy = 0
    for attachpoint_local, vector in self.forces:
      attachpoint_global = self.rotmat@attachpoint_local + self.position
      energy -= dot(attachpoint_global, vector)
    return energy

  def line(self):
    return [self.rotmat@array(point) + self.position for point in self.points]
