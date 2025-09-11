from constants import gamma, kappa_theta
from functions import angle
from numpy import array
from numpy.linalg import norm

class Flexure:
  '''
  Class containing the data of a PRBM flexure and the energy function

  Contains:
  - data, updated by the attached body.move()
  - energy(), returns the potential energy of the flexure
  '''
  def __init__(self, bodyA, attachpoint_localA, bodyB, attachpoint_localB):
    # constant
    self.attachpoint_localA = array(attachpoint_localA)
    self.attachpoint_localB = array(attachpoint_localB)

    # do not need to be remembered besides len0, which is constant
    attachpoint_globalA0 = bodyA.position + self.attachpoint_localA
    attachpoint_globalB0 = bodyB.position + self.attachpoint_localB

    vector_global0 = attachpoint_globalB0 - attachpoint_globalA0
    self.len0 = norm(vector_global0)
    unitvector_global0 = vector_global0/self.len0

    springpoint_globalA0 = attachpoint_globalA0 + (1 - gamma)/2*vector_global0
    springpoint_globalB0 = attachpoint_globalB0 - (1 - gamma)/2*vector_global0

    # constant
    self.springlen0 = norm(springpoint_globalB0 - springpoint_globalA0)

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

  def energy(self, A, E, I):
    # calculate and return energy of the flexure
    kappa = gamma*kappa_theta*E*I/self.len0
    k = E*A/self.springlen0
    vector_global = self.springpoint_globalB - self.springpoint_globalA
    springlen = norm(vector_global)
    unitvector_global = vector_global/springlen

    thetaA = angle(self.unitvector_globalA, unitvector_global)
    thetaB = angle(self.unitvector_globalB, unitvector_global)

    energyA = kappa*thetaA**2/2
    energyB = kappa*thetaB**2/2

    deltaS = springlen - self.springlen0
    energyAB = k*deltaS**2/2

    energy = energyA + energyB + energyAB

    return energy