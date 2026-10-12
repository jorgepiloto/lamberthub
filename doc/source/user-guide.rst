.. _user-guide:

User Guide
##########

Using a Solver
==============

Each Lambert's problem algorithm in ``lamberthub`` is a Python function with
the same interface:

.. code-block:: python

   from lamberthub import authorYYYY

   v1, v2 = authorYYYY(
       mu, r1, r2, tof, M=0, is_prograde=True, is_low_path=True,  # Type of solution
       maxiter=35, atol=1e-5, rtol=1e-7, full_output=False        # Iteration config
   )


Replace ``authorYYYY`` with the solver author's surname and its publication
year. For example, import the Izzo solver as ``izzo2015``. The
``lamberthub.ALL_SOLVERS`` list contains every available solver.

Parameters
==========

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Parameter
     - Type
     - Description
   * - ``mu``
     - ``float``
     - Gravitational parameter, which is the mass of the attracting body
       multiplied by the gravitational constant.
   * - ``r1``
     - ``numpy.ndarray``
     - Initial position vector.
   * - ``r2``
     - ``numpy.ndarray``
     - Final position vector.
   * - ``tof``
     - ``float``
     - Time of flight between the initial and final position vectors.
   * - ``M``
     - ``int``
     - Number of revolutions. A value of zero, the default, selects a direct
       transfer.
   * - ``is_prograde``
     - ``bool``
     - Controls the inclination of the final orbit. ``True`` selects an
       inclination between 0 and 90 degrees. ``False`` selects an inclination
       between 90 and 180 degrees.
   * - ``is_low_path``
     - ``bool``
     - Selects the path when more than two solutions are available. Neither
       path has an advantage unless the mission has specific constraints.
   * - ``maxiter``
     - ``int``
     - Maximum number of iterations allowed when computing the solution.
   * - ``atol``
     - ``float``
     - Absolute tolerance for the iterative method.
   * - ``rtol``
     - ``float``
     - Relative tolerance for the iterative method.
   * - ``full_output``
     - ``bool``
     - Returns iteration details when set to ``True``.

Returns
=======

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Return
     - Type
     - Description
   * - ``v1``
     - ``numpy.ndarray``
     - Initial velocity vector.
   * - ``v2``
     - ``numpy.ndarray``
     - Final velocity vector.
   * - ``numiter``
     - ``int``
     - Number of iterations. This value is returned only when ``full_output``
       is ``True``.
   * - ``tpi``
     - ``float``
     - Time per iteration. This value is returned only when ``full_output`` is
       ``True``.
