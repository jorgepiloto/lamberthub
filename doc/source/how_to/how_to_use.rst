How to use
==========

Once that you have installed the library, its time for you to learn how to use
it. The main goal of ``lamberthub`` is very simple: provide a collection of
algorithms for solving the Lambert's problem. 

All the routines are implemented in the form of Python functions, who's name is
given by the combination of original author's name plus the year of publication,
that is ``authorYYYY``.

Checking for available solvers
------------------------------
To answer this question, simply run the following code snippet or refer to the
official package API reference:

.. jupyter-execute::

    from lamberthub import ALL_SOLVERS
    print([solver.__name__ for solver in ALL_SOLVERS])

In addition, `lamberthub` provides other lists holding algorithms which
present particular features such as multi-revolutions or high-robustness. These
macros are listed down:

.. jupyter-execute::

    from lamberthub import ALL_SOLVERS, ZERO_REV_SOLVERS, MULTI_REV_SOLVERS, ROBUST_SOLVERS

Import a particular solver
--------------------------
If you are only interested in using a particular solver, you can easily import
it by running:

.. code-block:: python

    from lamberthub import authorYYYY

where ``author`` is the name of the author which developed the solver and ``YYYY``
the year of publication. Any of the solvers hosted in the ``ALL_SOLVERS`` macro
can be used.

If you would like to use a solver which is not defined in ``lamberthub``, open
a ``solver request`` in the `issues board
<https://github.com/jorgepiloto/lamberthub/issues>`_ detailing all the
information related to the algorithm and any useful reference which can help to
implement it.

A real example
--------------
The following section presents a real example [1]_. Suppose you want to solve for
the orbit of an interplanetary vehicle (that is Sun is the main attractor) form
which you know that the initial and final positions are given by:

.. math::

    \vec{r_1} = \begin{bmatrix}
    0.159321004\\
    0.579266185\\
    0.052359607\\
    \end{bmatrix} \text{[AU]}\;\;\;\;\;\;
    \vec{r_2} = \begin{bmatrix}
    0.057594337\\
    0.605750797\\
    0.068345246\\
    \end{bmatrix} \text{[AU]}

the dimension of previous vectors is astronomical units [AU] and the time of
flight, given in years, is known to be :math:`\Delta t = 0.010794065
\text{[year]}`.  The orbit is prograde since inclination is less than
:math:`90^{\circ}`) and direct :math:`M=0`. Remember that when :math:`M=0`,
there is only one possible solution, so the ``low_path`` flag does not play any
role in this problem.

To solve for the problem, first import a solver. For this problem,
``gooding1990`` is chosen:

.. jupyter-execute::
    :hide-output:

    from lamberthub import gooding1990

Next, specify the initial conditions of the problem:

.. jupyter-execute::
    :hide-output:

    # Import NumPy for declaring position vectors
    import numpy as np

    # Initial conditions for the problem
    mu_sun = 39.47692641  # [AU ** 3 / year ** 2]
    r1 = np.array([0.159321004, 0.579266185, 0.052359607])  # [AU]
    r2 = np.array([0.057594337, 0.605750797, 0.068345246])  # [AU]
    tof = 0.010794065  # [year]

Finally, the problem can be solved. Notice that, as explained before, the
default value for the ``prograde`` flag is ``True``, which matches the one from
problem's statement.

.. jupyter-execute::

    # Solving the problem
    v1, v2 = gooding1990(mu_sun, r1, r2, tof)
    
    # Let us print the results
    print(f"Initial velocity: {v1} [AU / years]")
    print(f"Final velocity:   {v2} [AU / years]")

previous values are the same ones coming from the original example.


.. [1] Directly taken from *An Introduction to the Mathematics and Methods of
  Astrodynamics, revised edition*, by R.H. Battin, problem 7-12.

