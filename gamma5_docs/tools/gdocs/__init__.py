"""Build code of the gamma5 documentation site (see ``DATA_FORMAT.md`` at the root of ``gamma5_docs``).

Modules:

* :mod:`gdocs.mathview`: the site notation (a registry of Symbolica symbols and function heads with MathML
  printers), the math object ``M(...)``, rich text ``prose(...)``, the converters from ``hepkit_renorm`` objects and
  the assets of the expression viewer (``site/assets/viewer``).
* :mod:`gdocs.demo_gallery`: the gallery ``site/assets/viewer/demo.html`` of the viewer.

Importing :mod:`gdocs.mathview` needs a licensed Symbolica (``SYMBOLICA_LICENSE`` in the environment).
"""

__all__ = ["mathview"]
