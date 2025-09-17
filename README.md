Code for modelling up 2D and 3D compliant mechanisms using PseudoRigidBody Models.

Check the examples folder for usage examples.

In general, place all .py files in one folder. Then run a new .py file within that folder with:

```python
from prbm import PRBM

p = PRBM(2) # for a 2D model

q = PRBM(3) # for a 3D model
```
Dependencies: numpy, scipy, matplotlib, and plotly for 3D plots.
