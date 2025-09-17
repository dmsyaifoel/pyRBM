Code for modelling up 2D and 3D compliant mechanisms using PseudoRigidBody Models.

Check the examples folder for usage examples.

In general, place all .py files in one folder. Then run a new .py file within that folder with:

```python
from prbm import PRBM

p = PRBM(2) # for a 2D model

q = PRBM(3) # for a 3D model
```
Dependencies: numpy, scipy, matplotlib.

For 3D plots, you additionally need plotly.

Useful functions:
```python
# q.add_body(bodyname, position) adds a body with a given name to a given position. For example
q.add_body('A') # leaving out the position will place it at the origin
q.add_body('B', (0, 0, 2*cm)) # body placed 2 cm above body A
q.add_body('C', (0, 0, 4*cm)) # body placed 4 cm above body A

# q.add_flexure(bodynameA, attachpoint_localA, bodynameB, attachpoint_localB) adds a flexure between two bodies
# the attachment points are set in the bodies' local coordinate systems. This allows for easily making vertical flexure like this, for example
q.add_flexure('A', (-1*cm, 0, 0), 'B', (-1*cm, 0, 0))
q.add_flexure('A', (1*cm, 0, 0), 'B', (1*cm, 0, 0))
q.add_flexure('A', (0, 1*cm, 0), 'B', (0, 1*cm, 0))
q.add_flexure('C', (-1*cm, 0, 0), 'B', (-1*cm, 0, 0))
q.add_flexure('C', (1*cm, 0, 0), 'B', (1*cm, 0, 0))
q.add_flexure('C', (0, 1*cm, 0), 'B', (0, 1*cm, 0))

# q.move(bodyname, position, angles) moves a body, deforming any attached flexures, for example
q.move('B', (5*mm, 0, 2*cm)) # leaving out the angles will keep the body in the same orientation
q.move('C', (5*mm, 0, 4*cm), (.001, .001, .001) # otherwise the angles are in radians as XYZ euler angles

# you can move a body then attach new flexures to make a prestressed compliant mechanism

# q.add_force(bodyname, vector, attachpoint_local) # adds a force to a body, in Newton
q.add_force('B', (1, 0, 0)) # leaving out attachpoint_local will place the force at the body's origin

# q.solve_pose(bodynames, A, E, I, method) moves all bodies mentioned in bodynames to the pose with the least energy, for example
from numpy import pi
r = 1*mm
A = pi*r**2 # cross section area of a flexure, in meter**2
E = 1650e6 # E modulus of the material, in Pa
I = pi/2*r**4 # second moment of area of a flexure, in m**4
q.solve_pose('B', A, E, I) # will only move B. Leaving method empty will use the default method of scipy.optimize.minimize

# Note that if you have multiple 'free' bodies to solve, you can use p.solve_pose('BCDEF', A, E, I) if each bodyname is a single letter
# Otherwise use a list or tuple like p.solve_pose(['intermediate body 1', 'intermediate body 2', 'lever', 'end effector'], A, E, I)

# Lastly you can output the positions using
q.print()

# And visualize the PRBM
q.show()
```
