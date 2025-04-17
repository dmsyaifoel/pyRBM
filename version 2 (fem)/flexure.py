from functions import norm, angle
from constants import gamma, kappa_theta

class Flexure:
  def __init__(self, bodyA, bodyB):
    attachpoint_globalA0 = bodyA.position
    attachpoint_globalB0 = bodyB.position

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

  @property
  def line(self):
    return [self.attachpoint_globalA,
            self.springpoint_globalA,
            self.springpoint_globalB,
            self.attachpoint_globalB]

  def energy(self, A, E, I):
    vector_global = self.springpoint_globalB - self.springpoint_globalA
    springlen = norm(vector_global)
    unitvector_global = vector_global/springlen

    thetaA = angle(self.unitvector_globalA, unitvector_global)
    thetaB = angle(self.unitvector_globalB, unitvector_global)

    kappa = gamma*kappa_theta*E*I/self.len0
    k = E * A / self.springlen0

    energyA = kappa*thetaA**2/2
    energyB = kappa*thetaB**2/2

    deltaS = springlen - self.springlen0
    energyAB = k*deltaS**2/2

    energy = energyA + energyB + energyAB
    return energy