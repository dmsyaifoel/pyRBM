from flexure import *

class Body:
  def __init__(self, name, position0, dim):
    self.name = name
    self.position0 = np.array(position0)
    self.position = np.array(position0)
    self.dim = dim

    if self.dim == 2:
      self.angles0 = 0
      self.angles = 0

    if self.dim == 3:
      self.angles0 = np.zeros(3)
      self.angles = np.zeros(3)

    self.rotmat = rotmat(self.angles, self.dim)

    self.flexures = []
    self.letters = []
    self.has_moved = False

    self.forces = []
    self.reference_points = []

  def move(self, position, angles=None): 
    # move a body.py to a pose
    self.has_moved = True
    self.position = np.array(position)
    if angles is None: self.angles = self.angles0
    else: self.angles = angles

    self.rotmat = rotmat(self.angles, self.dim)

    # then update all flexure attached
    for flexure, letter in zip(self.flexures, self.letters):
      if letter == 'A':
        flexure.attachpoint_globalA = self.rotmat@flexure.attachpoint_localA + self.position
        flexure.springpoint_globalA = self.rotmat@flexure.springpoint_localA + self.position
        flexure.unitvector_globalA = self.rotmat@flexure.unitvector_localA
        if self.dim == 3 and enabletorsion == True: flexure.orthvector_globalA = self.rotmat@flexure.orthvector_global

      if letter == 'B':
        flexure.attachpoint_globalB = self.rotmat@flexure.attachpoint_localB + self.position
        flexure.springpoint_globalB = self.rotmat@flexure.springpoint_localB + self.position
        flexure.unitvector_globalB = self.rotmat@flexure.unitvector_localB
        if self.dim == 3 and enabletorsion == True: flexure.orthvector_globalB = self.rotmat@flexure.orthvector_global

  def moveback(self): 
    # move bodies back to their original position
    self.move(self.position0, self.angles0)
    self.has_moved = False

  def energy(self): 
    # calculate the potential energy of a body.py due to external forces
    energy = 0
    for attachpoint_local, vector in self.forces:
      attachpoint_global = self.rotmat@attachpoint_local + self.position
      energy -= np.dot(attachpoint_global, vector)
    return energy

  def line(self): 
    # return a list of interesting points to plot
    line = []

    for point in self.reference_points: line.append(self.rotmat@point + self.position)

    for flexure, letter in zip(self.flexures, self.letters):
      if letter == 'A': line.append(flexure.attachpoint_globalA)
      else: line.append(flexure.attachpoint_globalB)

    return np.array(line).T
