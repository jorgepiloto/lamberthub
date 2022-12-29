How to install
==============

This page contains the guide for installing the ``lamberthub`` package. If you
experience any kind of problem during one of the steps shown in the following
lines, please open an new issue (or select similar ones) in the `issues board
<https://github.com/jorgepiloto/lamberthub/issues>`_.

Install from PyPI using pip
---------------------------

The installation process is similar to other python packages, meaning that you
only need to run:

.. code-block:: text

    python -m pip install lamberthub


Previous command installs the latest stable version of the library. Once
done, you can open the Python terminal and import the package and verify its
version by running:

.. jupyter-execute::

    import lamberthub
    print(f"Current lamberthub version is {lamberthub.__version__}")
