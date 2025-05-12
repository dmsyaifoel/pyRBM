from functions import Matrix, R, dot

class Body:
  def __init__(self, position0):
    self.position0 = Matrix(position0)
    self.position = Matrix(position0)

    self.angle0 = 0
    self.angle = 0
    self.rotmat = R(self.angle)

    self.flexures = []
    self.which = []

    self.clear_force()

  def clear_force(self):
    self.force = Matrix((0, 0))

  def move(self, position, angle=None):
    self.position = Matrix(position)
    if angle is None: self.angle = self.angle0
    else: self.angle = angle

    self.rotmat = R(self.angle)

    for flexure, which in zip(self.flexures, self.which):
      if which:
        flexure.attachpoint_globalA = self.position
        flexure.springpoint_globalA = self.rotmat@flexure.springpoint_localA + self.position
        flexure.unitvector_globalA = self.rotmat@flexure.unitvector_localA

      else:
        flexure.attachpoint_globalB = self.position
        flexure.springpoint_globalB = self.rotmat@flexure.springpoint_localB + self.position
        flexure.unitvector_globalB = self.rotmat@flexure.unitvector_localB

  @property
  def energy(self):
    energy = -dot(self.position, self.force)
    return energy