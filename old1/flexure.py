from functions import *
from settings import *
from properties import *

class Flexure:
  def __init__(self, name, bodyA, attachpoint_localA, bodyB, attachpoint_localB, properties, dim):

    # constant
    self.name = name
    self.gamma = gamma
    self.kappa_theta = kappa_theta
    self.E = properties['E']
    self.I = properties['I']
    self.A = properties['A']
    self.dim = dim

    if self.dim == 3 and enabletorsion:
      self.G = properties['G']
      self.J = properties['J']

    # constant
    self.attachpoint_localA = np.array(attachpoint_localA)
    self.attachpoint_localB = np.array(attachpoint_localB)

    # do not need to be remembered besides len0, which is constant
    attachpoint_globalA0 = bodyA.position + attachpoint_localA
    attachpoint_globalB0 = bodyB.position + attachpoint_localB

    vector_global0 = attachpoint_globalB0 - attachpoint_globalA0
    self.len0 = sl.norm(vector_global0)
    unitvector_global0 = vector_global0/self.len0

    springpoint_globalA0 = attachpoint_globalA0 + (1 - self.gamma)/2*vector_global0
    springpoint_globalB0 = attachpoint_globalB0 - (1 - self.gamma)/2*vector_global0

    # constant
    self.springlen0 = sl.norm(springpoint_globalB0 - springpoint_globalA0)

    # variable, updated by Body.move()
    self.attachpoint_globalA = attachpoint_globalA0
    self.attachpoint_globalB = attachpoint_globalB0

    self.springpoint_globalA = springpoint_globalA0
    self.springpoint_globalB = springpoint_globalB0

    # constant
    self.springpoint_localA = springpoint_globalA0 - bodyA.position
    self.springpoint_localB = springpoint_globalB0 - bodyB.position

    # variable, updated by Body.move()
    self.unitvector_globalA = unitvector_global0
    self.unitvector_globalB = unitvector_global0

    # constant
    self.unitvector_localA = unitvector_global0
    self.unitvector_localB = unitvector_global0

    self.kappa = self.gamma*self.kappa_theta*self.E*self.I/self.len0
    self.k = self.E*self.A/self.springlen0

    if self.dim == 3 and enabletorsion:
      self.kappatorsion = self.J*self.G/self.springlen0
      orthvector = np.cross(np.random.random(3), unitvector_global0)
      self.orthvector_global = normalize(orthvector) # constant
      self.orthvector_globalA = normalize(orthvector) # updated by Body.move()
      self.orthvector_globalB = normalize(orthvector) # updated by Body.move()

  def line(self): 
    # return a list of points to plot
    return np.array((self.attachpoint_globalA, self.springpoint_globalA,
                     self.springpoint_globalB, self.attachpoint_globalB)).T

  def energy(self): 
    # calculate and return energy of the flexure
    vector_global = self.springpoint_globalB - self.springpoint_globalA
    springlen = sl.norm(vector_global)
    unitvector_global = vector_global/springlen

    thetaA = angle(self.unitvector_globalA, unitvector_global)
    thetaB = angle(self.unitvector_globalB, unitvector_global)

    energyA = self.kappa*thetaA**2/2
    energyB = self.kappa*thetaB**2/2

    deltaS = springlen - self.springlen0
    energyAB = self.k*deltaS**2/2

    energy = energyA + energyB + energyAB

    if self.dim == 3 and enabletorsion:
      RA = Rab(self.unitvector_globalA, unitvector_global)
      RB = Rab(self.unitvector_globalB, unitvector_global)
      rotatedvectorA = RA@self.orthvector_globalA
      rotatedvectorB = RB@self.orthvector_globalB
      torsionangle = angle(rotatedvectorA, rotatedvectorB)
      energytorsion = self.kappatorsion*torsionangle**2/2
      return energy + energytorsion

    return energy
