API
===

This section holds the documentation for the source code of ``lamberthub``. You
can navigate through each of of the sub-packages and modules they contain to
check the different auxiliary routines and parameters they accept.


.. toctree::
   :titlesonly:
   :maxdepth: 2

   {% for page in pages %}
   {% if (page.top_level_object or page.name.split('.') | length == 3) and page.display %}
   {{ page.include_path }}
   {% endif %}
   {% endfor %}
