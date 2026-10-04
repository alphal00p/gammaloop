# type: ignore
from __future__ import absolute_import
from .object_library import all_parameters, Parameter

from .function_library import complexconjugate, re, im, csc, sec, acsc, asec

import os

# Specify here how many scalars need to be considered
if os.environ.get('UFO_SCALARS_MODEL_N_SCALARS') is not None:
    try:
        N_SCALARS = int(os.environ['UFO_SCALARS_MODEL_N_SCALARS'])
    except ValueError:
        raise ValueError(
            "Environment variable UFO_SCALARS_MODEL_N_SCALARS must be an integer")
    print("Loading UFO scalars model with %d scalars" % N_SCALARS)
else:
    N_SCALARS = 3

# Add all possible three-point interactions. User can filter undesired ones in vertices.py if needed.
if os.environ.get('UFO_SCALARS_MODEL_N_POINT_INTERACTIONS') is not None:
    try:
        N_POINT_INTERACTIONS = [
            int(x) for x in os.environ['UFO_SCALARS_MODEL_N_POINT_INTERACTIONS'].split(',')]
    except ValueError:
        raise ValueError(
            "Environment variable UFO_SCALARS_MODEL_N_POINT_INTERACTIONS must be a comma-separated list of integers")
    print("Loading UFO scalars model with n-point interactions, n=[%s]" %
          '|'.join(str(x) for x in N_POINT_INTERACTIONS))
else:
    N_POINT_INTERACTIONS = [3, 4, 5, 6, 7, 8, 9, 10]
# N_POINT_INTERACTIONS = [3,4]

# This is a default parameter object representing 0.
ZERO = Parameter(name='ZERO',
                 nature='internal',
                 type='real',
                 value='0.0',
                 texname='0')

# Introduce the scalar coupling lambda
lam = Parameter(name='lam',
                nature='external',
                type='real',
                value=1.0,
                texname='\\lambda',
                lhablock='SCALARS',
                lhacode=[1337])

# Assign masses and widths to the scalars
scalar_mass_parameters = []
scalar_width_parameters = []
for i_scalar in range(N_SCALARS):
    # Set the first scalar with id 0 to always be massless
    if i_scalar == 0:
        scalar_mass_parameters.append(ZERO)
        scalar_width_parameters.append(ZERO)
        continue
    # remember that the constructor automatically adds the objects to all_parameters
    scalar_mass_parameters.append(Parameter(name='mass_scalar_%s' % i_scalar,
                                            nature='external',
                                            type='real',
                                            value=1.0*i_scalar,
                                            texname='m_%d' % i_scalar,
                                            lhablock='SCALARS',
                                            lhacode=[2001+i_scalar]))
    scalar_width_parameters.append(Parameter(name='width_scalar_%s' % i_scalar,
                                             nature='external',
                                             type='real',
                                             value=0.0,
                                             texname='w_%d' % i_scalar,
                                             lhablock='SCALARS',
                                             lhacode=[3001+i_scalar]))

# User-defined parameters. Keep key SM parameters as their use may be hardcoded in some generators
aEWM1 = Parameter(name='aEWM1',
                  nature='external',
                  type='real',
                  value=132.50698,
                  texname='\\text{aEWM1}',
                  lhablock='SMINPUTS',
                  lhacode=[1])

Gf = Parameter(name='Gf',
               nature='external',
               type='real',
               value=0.0000116639,
               texname='G_f',
               lhablock='SMINPUTS',
               lhacode=[2])

aS = Parameter(name='aS',
               nature='external',
               type='real',
               value=0.118,
               texname='\\alpha _s',
               lhablock='SMINPUTS',
               lhacode=[3])

aEW = Parameter(name='aEW',
                nature='internal',
                type='real',
                value='1/aEWM1',
                texname='\\alpha _{\\text{EW}}')

G = Parameter(name='G',
              nature='internal',
              type='real',
              value='2*cmath.sqrt(aS)*cmath.sqrt(cmath.pi)',
              texname='G')

ee = Parameter(name='ee',
               nature='internal',
               type='real',
               value='2*cmath.sqrt(aEW)*cmath.sqrt(cmath.pi)',
               texname='e')
