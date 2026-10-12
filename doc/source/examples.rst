.. _examples:

Examples
########

Direct Transfer
===============

Problem Statement
-----------------

Consider an interplanetary vehicle with the Sun as its main attractor. Its
initial and final position vectors are:

.. math::

   \vec{r_1} =
   \begin{bmatrix}
   0.159321004 \\
   0.579266185 \\
   0.052359607
   \end{bmatrix}
   \text{ [AU]}
   \qquad
   \vec{r_2} =
   \begin{bmatrix}
   0.057594337 \\
   0.605750797 \\
   0.068345246
   \end{bmatrix}
   \text{ [AU]}

The time of flight is :math:`\Delta t = 0.010794065` years. The orbit is
prograde and direct, so :math:`M = 0`. A direct transfer has only one solution,
which means the ``is_low_path`` option does not affect this problem.

Solution
--------

This example uses ``gooding1990``. Other solvers that support direct transfers
can solve the same problem.

.. testcode::

   import numpy as np

   from lamberthub import gooding1990

   mu_sun = 39.47692641
   r1 = np.array([0.159321004, 0.579266185, 0.052359607])
   r2 = np.array([0.057594337, 0.605750797, 0.068345246])
   tof = 0.010794065

   v1, v2 = gooding1990(mu_sun, r1, r2, tof, M=0, is_prograde=True)
   print(f"Initial velocity: {v1} [AU / years]")
   print(f"Final velocity:   {v2} [AU / years]")

Result
------

.. testoutput::

   Initial velocity: [-9.303608    3.01862016  1.53636008] [AU / years]
   Final velocity:   [-9.5111862   1.88884006  1.4213781 ] [AU / years]

This example comes from *An Introduction to the Mathematics and Methods of
Astrodynamics, revised edition*, by R. H. Battin, problem 7-12.
