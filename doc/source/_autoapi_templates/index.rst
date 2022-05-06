:noindex:
:orphan:

Package API
===========

This section holds `lamberthub's` source code documentation. You can navigate
through each of of the sub-packages and modules they contain to check the
different auxiliary routines and parameters they accept.

.. toctree::
   :maxdepth: 2

   {% for page in pages %}
   {% if page.top_level_object and page.display %}
   {{ page.include_path }}
   {% endif %}
   {% endfor %}
