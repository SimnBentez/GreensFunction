"""
GreensFunction
==============

High-performance Green function beam solver implemented in C++
with Python bindings.

"""

from .greensfunction import BeamSolver as _BeamSolver

__version__ = "0.1.0"

class BeamSolver(_BeamSolver):
    """
    Python wrapper for the C++ BeamSolver.

    Parameters
    ----------
    E : float
        Young's modulus
    I : float
        Second moment of area
    L : float
        Beam length
    """

    pass


__all__ = ["BeamSolver"]