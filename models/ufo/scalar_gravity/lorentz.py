# type: ignore
from __future__ import absolute_import
from .object_library import all_lorentz, Lorentz

from .function_library import complexconjugate, re, im, csc, sec, acsc, asec
from . import parameters as Param

try:
    import form_factors as ForFac
except ImportError:
    pass

SSTmpart = Lorentz(
    name='SSTmpart',
    spins=[1, 1, 5],
    structure='-Metric(1003,2003)')

SST = Lorentz(
    name='SST',
    spins=[1, 1, 5],
    structure='-P(-1,1)*P(-1,2)*Metric(1003,2003) + P(1003,1)*P(2003,2) + P(2003,1)*P(1003,2)')

# TODO
SSTTmpart = Lorentz(
    name='SSTTmpart',
    spins=[1, 1, 5, 5],
    structure='Metric(1003,1004)*Metric(2003,2004)')

# TODO
SSTT = Lorentz(
    name='SSTT',
    spins=[1, 1, 5, 5],
    structure='Metric(1003,1004)*Metric(2003,2004)')

# TODO
SSTTTmpart = Lorentz(
    name='SSTTTmpart',
    spins=[1, 1, 5, 5, 5],
    structure='Metric(1003,1004)*Metric(2003,2004)*P(1005,1)*P(2005,2)')

# TODO
SSTTT = Lorentz(
    name='SSTTT',
    spins=[1, 1, 5, 5, 5],
    structure='Metric(1003,1004)*Metric(2003,2004)*P(1005,1)*P(2005,2)')

# TODO
TTT = Lorentz(
    name='TTT',
    spins=[5, 5, 5],
    structure='Metric(1001,1002)*Metric(2001,2002)*Metric(1003,2003)')

# TODO
TTTT = Lorentz(
    name='TTTT',
    spins=[5, 5, 5, 5],
    structure='Metric(1001,1002)*Metric(2001,2002)*Metric(1003,2003)*Metric(1004,2004)')
