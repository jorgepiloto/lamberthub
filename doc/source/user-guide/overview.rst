Overview
========

Once you have installed the library, it's time to learn how to use it. The main
goal of ``lamberthub`` is very simple: provide a collection of algorithms for
solving Lambert's problem.

All the routines are implemented in the form of Python functions, whose name is
given by the combination of the original author's name plus the year of
publication, that is: ``authorYYYY``.

Checking for available solvers
------------------------------

To answer this question, simply run the following code snippet or refer to the
official package API reference:

.. code-block:: python

    from lamberthub import ALL_SOLVERS
    print([solver.__name__ for solver in ALL_SOLVERS])

In addition, ``lamberthub`` provides other lists holding algorithms which
present particular features such as multi-revolutions or high-robustness. These
macros are listed down:

.. code-block:: python

    from lamberthub import ALL_SOLVERS, ZERO_REV_SOLVERS, MULTI_REV_SOLVERS, ROBUST_SOLVERS

Import a particular solver
--------------------------

If you are only interested in using a particular solver, you can easily import
it by running:

.. code-block:: python

    from lamberthub import authorYYYY

Where ``author`` is the name of the author who developed the solver and ``YYYY``
is the year of publication. Any of the solvers hosted in the ``ALL_SOLVERS``
macro can be used.

If you would like to use a solver which is not defined in ``lamberthub``, open a
``solver request`` in the ``issues board`` detailing all the information related
to the algorithm and any useful reference which can help to implement it.

Input and output parameters
---------------------------

Any of the routines in the library has the same number of input and output
parameters because they all solve the same astrodynamics problem.

As said before, any Lambert's problem algorithm implemented in ``lamberthub`` is
a Python function built with the following API architecture:

**Parameters**

- `mu`: the gravitational parameter, that is the mass of the attracting body times the gravitational constant.
- `r1`: initial position vector, needs to be a NumPy array instance.
- `r2`: final position vector, needs to be a NumPy array instance.
- `tof`: time of flight between initial and final vectors.

**Additional parameters**

- `M`: the number of desired revolutions. If zero, direct transfer is assumed.

- `prograde`: this parameter controls the inclination of the final orbit. If set
  to `True`, the transfer has an inclination between 0 and 90 degrees while if
  `False` inclinations between 90 and 180 are provided.

- `low_path`: selects the type of path when more than two solutions are
  available. There is no actual advantage on one or another solution, unless you
  have particular constrains on your mission. If `True`, the maximum value for
  the independent variable (when two solutions) is selected.

- `maxiter`: maximum number of iterations allowed when computing the solution.

- `atol`: absolute tolerance for the iterative method.

- `rtol`: relative tolerance for the iterative method.

- `full_output`: if `True`, it returns additional information such as the number
  of iterations.

**Returns**

- `v1`: initial velocity vector, it is a NumPy array instance.

- `v2`: final velocity vector, it is a NumPy array instance.

**Additional returns**

- `numiter`: number of iterations. Only if `full_output` has been set to `True`.

- `tpi`: time per iteration. Only if `full_output` has been set to `True`.

A real example
--------------

The following section presents a real example[^1]. Suppose you want to solve for
the orbit of an interplanetary vehicle (that is the Sun is the main attractor)
for which you know that the initial and final positions are given by:

.. math::

    \vec{r_1} = \begin{bmatrix} 0.159321004 \\ 0.579266185 \\ 0.052359607 \end{bmatrix} \text{[AU]}\;\;\;\;\;\;
    \vec{r_2} = \begin{bmatrix} 0.057594337 \\ 0.605750797 \\ 0.068345246 \end{bmatrix} \text{[AU]}

The dimension of previous vectors is astronomical units [AU], and the time of
flight, given in years, is known to be $\Delta t = 0.010794065 \text{[year]}$.
The orbit is prograde since the inclination is less than $90^{\circ}$, and
direct $M=0$. Remember that when $M=0$, there is only one possible solution, so
the ``low_path`` flag does not play any role in this problem.

To solve the problem, first import a solver. For this problem, ``gooding1990`` is
chosen:

.. code-block:: python

    from lamberthub import gooding1990
    import numpy as np
    
    
    # Initial conditions for the problem
    mu_sun = 39.47692641  # [AU ** 3 / year ** 2]
    r1 = np.array([0.159321004, 0.579266185, 0.052359607])  # [AU]
    r2 = np.array([0.057594337, 0.605750797, 0.068345246])  # [AU]
    tof = 0.010794065  # [year]
    
    # Solve the problem using Gooding's 1990 algorithm
    v1, v2 = gooding1990(mu_sun, r1, r2, tof)
    
    # Print the results
    print(f"Initial velocity: {v1} [AU / years]\nFinal velocity: {v2} [AU / years]")

.. code-block:: pycon

   Initial velocity: [-9.303608    3.01862016  1.53636008] [AU / years]
   Final velocity: [-9.5111862   1.88884006  1.4213781 ] [AU / years]

Previous values are the same ones stated in original example.
