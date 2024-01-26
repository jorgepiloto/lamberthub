lamberthub documentation |version|
##################################

The ``lamberthub`` library is a Python package which provides its users with
different algorithms for solving the Lambert's problem:

.. figure:: _static/lamberts_problem_geometry.png
   :width: 350px
   :align: center

.. grid:: 2

    .. grid-item-card:: Getting started :fa:`person-running`
        :link: getting-started/index
        :link-type: doc

        Step by step guidelines on how to set up your environment to work with
        lamberthub. Install the library and verify it works as expected.

    .. grid-item-card:: User guide :fa:`book-open-reader`
        :link: user-guide/index
        :link-type: doc

        Learn about the capabilities, features, and key topics in lamberthub.
        This guide provides useful background information and explanations.

Since the formulation of the problem, many different solutions have been
devised, being the latest ones in the form of computer routines. By collecting
and implementing all of them under a common programming language, it is possible
to carry out performance comparisons between those. Furthermore, ``lamberthub``
provides a framework for new authors to select the robustness and accuracy of
their routines.

.. jinja:: main_toctree

    .. toctree::
       :hidden:
       :maxdepth: 3
    
       getting-started/index
       user-guide/index
       {% if build_examples %}
       examples
       {% endif %}
       {% if build_api %}
       api/index
       {% endif %}
