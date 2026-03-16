# %% [markdown]
# # Implementation of the IR-subtracted process $d \bar{d} \rightarrow \gamma \gamma \gamma$

# %%
# ENVIRONMENT SETUP

# Tested with gammaloop with rev c71701cefc11eb7184234d8d70a190222374c32b 
# Tested with symbolica-community v 1.1.0 (rev b796597474e3eee41a463b3b816df169220cc343)

import os
ROOT_PATH = os.path.dirname(__file__) 

SYMBOLICA_COMMUNITY_PATH = None # If not installed for the current python intepreter, specified e.g. "/Users/vjhirsch/HEP_programs/symbolica-community/python"
ROOT_DIR = ROOT_PATH

GAMMALOOP_CLI_PATH = os.path.join('../../../..','target/dev-optim/gammaloop')
SM_MODEL_PATH = 'sm-default.json'

#ROOT_DIR = os.getcwd()
import os
os.environ['SYMBOLICA_LICENSE'] = 'dcec4a5e#6a95649c#7dca8216-8afe-57c8-975e-03eb5e68e4ee'

# %% [markdown]
# ## Dependencies

# %%
import sys
from pprint import pprint, pformat
import pydot
import itertools
import subprocess
import shutil
import math
from math import sqrt, log10
from IPython.display import Markdown as md # type: ignore

if SYMBOLICA_COMMUNITY_PATH is not None:
    sys.path.insert(0, SYMBOLICA_COMMUNITY_PATH)
from symbolica import Expression, S,E, Replacement, PrintMode
from symbolica.community.spenso import TensorName as \
    N,LibraryTensor,TensorNetwork,Representation,TensorStructure,TensorIndices,Tensor,Slot,TensorLibrary, ExecutionMode
from symbolica.community.idenso import *
from symbolica.community.spenso import Tensor
#from symbolica.community.spenso import initialize as spenso_initialize
#spenso_initialize()

#print(E("spenso::{spenso::upper}::bis(4,x)").to_canonical_string())
#from symbolica_community import Expression, S,E, Replacement
#from symbolica_community.tensors import TensorName as \
#    N,LibraryTensor,TensorNetwork,Representation,TensorStructure,TensorIndices,Tensor,Slot,TensorLibrary, ExecutionMode
# from symbolica_community.algebraic_simplification import *
# import symbolica_community.tensors as tensors

# %%
def run_gammaloop_commands(cmds, capture_output=True, start_clean=True, debug=False, state_name='TMP_state', toml_run_card=False, turn_off_debug_logs=False):
    """Run a list of gammaloop commands"""

    if start_clean:
        if os.path.isdir(os.path.join(ROOT_DIR, state_name)):
            shutil.rmtree(os.path.join(ROOT_DIR, state_name))
    
    env = os.environ
    if turn_off_debug_logs:
        env['GL_DISPLAY_FILTER'] = 'info'
        env['GL_LOGFILE_FILTER'] = 'info'

    if not toml_run_card:
        for cmd in cmds:
            if cmd.startswith('#'):
                continue
            if cmd.strip() == '':
                continue
            if cmd.startswith('!'):
                c = cmd[1:].strip()
            else:
                c = f'{GAMMALOOP_CLI_PATH} -s {os.path.join(ROOT_DIR, state_name)} -o {cmd}'
            if debug:
                print(f'Running command: {c}')
            r = subprocess.run(c, capture_output=capture_output, text=True, cwd=ROOT_DIR, shell=True, env=env)
            if r.returncode != 0:
                raise RuntimeError(f"Command failed with error: {r.stderr}")        
    else:
        c = f'{GAMMALOOP_CLI_PATH} -s {os.path.join(ROOT_DIR, state_name)} {cmds}'
        if debug:
            print(f'Running command: {c}')
        r = subprocess.run(c, capture_output=capture_output, text=True, cwd=ROOT_DIR, shell=True, env=env)
        if r.returncode != 0:
            raise RuntimeError(f"Command failed with error: {r.stderr}")

# Work around expressions given as strings containing the wrapping quotes
def Es(expr: str) -> Expression:
    return E(expr.replace('"',''), default_namespace="gammalooprs")

hep_lib = TensorLibrary.hep_lib()

# %%
# import symbolica
# print(symbolica.__file__)
# print()

# %% [markdown]
# ## HELAS polarization vectors

# %%

class WaveFunction(list):
    """a objet for a WaveFunction"""
    
    spin_to_size={0:1,
                  1:3,
                  2:6,
                  3:6,
                  4:18,
                  5:18}
    
    def __init__(self, spin= None, size=None):
        """Init the list with zero value"""
        
        if spin:
            size = self.spin_to_size[spin]
        list.__init__(self, [0]*size)
        

def ixxxxx(p,fmass,nhel,nsf):
    """Defines an inflow fermion."""
    
    fi = WaveFunction(2)
    
    fi[0] = complex(-p[0]*nsf,-p[3]*nsf)
    fi[1] = complex(-p[1]*nsf,-p[2]*nsf) 
    nh = nhel*nsf 
    if (fmass != 0.):
        pp = min(p[0],sqrt(p[1]**2 + p[2]**2 + p[3]**2 ))
        if (pp == 0.): 
            sqm = [sqrt(abs(fmass))]
            sqm.append(sign(sqm[0],fmass)) 
            ip = (1+nh)//2 
            im = (1-nh)//2 

            fi[2] = ip*sqm[ip]
            fi[3] = im*nsf*sqm[ip]
            fi[4] = ip*nsf*sqm[im]
            fi[5] = im*sqm[im]

        else:
            sf = [(1+nsf+(1-nsf)*nh)*0.5,(1+nsf-(1-nsf)*nh)*0.5]
            omega = [sqrt(p[0]+pp),fmass/(sqrt(p[0]+pp))]
            ip = (1+nh)//2
            im = (1-nh)//2
            sfomeg = [sf[0]*omega[ip],sf[1]*omega[im]]
            pp3 = max(pp+p[3],0.)
            if (pp3 == 0.):
                chi1 = complex(-nh,0.) 
            else:
                chi1 = complex(nh*p[1]/sqrt(2.*pp*pp3),\
                p[2]/sqrt(2.*pp*pp3))
            chi = [complex(sqrt(pp3*0.5/pp)),chi1]

            fi[2] = sfomeg[0]*chi[im]
            fi[3] = sfomeg[0]*chi[ip]
            fi[4] = sfomeg[1]*chi[im]
            fi[5] = sfomeg[1]*chi[ip] 
    
    else: 
        sqp0p3 = sqrt(max(p[0]+p[3],0.))*nsf
        if (sqp0p3 == 0.):
            chi1 = complex(-nhel*sqrt(2.*p[0]),0.)
        else:
            chi1 = complex(nh*p[1]/sqp0p3,p[2]/sqp0p3)
        chi = [complex(sqp0p3,0.),chi1]
        if (nh == 1):
            fi[2] = complex(0.,0.)
            fi[3] = complex(0.,0.)
            fi[4] = chi[0]
            fi[5] = chi[1] 
        else:
            fi[2] = chi[1]
            fi[3] = chi[0]
            fi[4] = complex(0.,0.)
            fi[5] = complex(0.,0.) 
    
    return fi 

def oxxxxx(p,fmass,nhel,nsf):
    """ initialize an outgoing fermion"""
    
    fo = WaveFunction(2)
    fo[0] = complex(p[0]*nsf,p[3]*nsf)
    fo[1] = complex(p[1]*nsf,p[2]*nsf)
    nh = nhel*nsf
    if (fmass != 0.):
        pp = min(p[0],sqrt(p[1]**2 + p[2]**2 + p[3]**2 ))
        if (pp == 0.): 
            sqm = [sqrt(abs(fmass))]
            sqm.append( sign(sqm[0], fmass)) 
            ip = int(-((1-nh)//2) * nhel)
            im = int((1+nh)//2 * nhel)
            
            fo[2] = im*sqm[abs(im)]
            fo[3] = ip*nsf*sqm[abs(im)]
            fo[4] = im*nsf*sqm[abs(ip)]
            fo[5] = ip*sqm[abs(ip)]

        else:
            sf = [(1+nsf+(1-nsf)*nh)*0.5,(1+nsf-(1-nsf)*nh)*0.5]
            omega = [sqrt(p[0]+pp),fmass/(sqrt(p[0]+pp))]
            ip = (1+nh)//2
            im = (1-nh)//2
            sfomeg = [sf[0]*omega[ip],sf[1]*omega[im]]
            pp3 = max(pp+p[3],0.)
            if (pp3 == 0.):
                chi1 = complex(-nh,0.) 
            else:
                chi1 = complex(nh*p[1]/sqrt(2.*pp*pp3),\
                -p[2]/sqrt(2.*pp*pp3))
            chi = [complex(sqrt(pp3*0.5/pp)),chi1]

            fo[2] = sfomeg[1]*chi[im]
            fo[3] = sfomeg[1]*chi[ip]
            fo[4] = sfomeg[0]*chi[im]
            fo[5] = sfomeg[0]*chi[ip] 
            
    else: 
        sqp0p3 = sqrt(max(p[0]+p[3],0.))*nsf
        if (sqp0p3 == 0.):
            chi1 = complex(-nhel*sqrt(2.*p[0]),0.)
        else:
            chi1 = complex(nh*p[1]/sqp0p3,-p[2]/sqp0p3)
        chi = [complex(sqp0p3,0.),chi1]
        if (nh == 1):
            fo[2] = chi[0]
            fo[3] = chi[1]
            fo[4] = complex(0.,0.)
            fo[5] = complex(0.,0.) 
        else:
            fo[2] = complex(0.,0.)
            fo[3] = complex(0.,0.)
            fo[4] = chi[1]
            fo[5] = chi[0] 
    
    return fo

def vxxxxx(p,vmass,nhel,nsv):
    """ initialize a vector wavefunction. nhel=4 is for checking BRST"""
    
    vc = WaveFunction(3)
    
    sqh = sqrt(0.5)
    nsvahl = nsv*abs(nhel)
    pt2 = p[1]**2 + p[2]**2 
    pp = min(p[0],sqrt(pt2 + p[3]**2))
    pt = min(pp,sqrt(pt2))

    vc[0] = complex(p[0]*nsv,p[3]*nsv)
    vc[1] = complex(p[1]*nsv,p[2]*nsv)

    if (nhel == 4):
        if (vmass == 0.):
            vc[2] = 1.
            vc[3]=p[1]/p[0]
            vc[4]=p[2]/p[0]
            vc[5]=p[3]/p[0]
        else:
            vc[2] = p[0]/vmass
            vc[3] = p[1]/vmass
            vc[4] = p[2]/vmass
            vc[5] = p[3]/vmass
        
        return vc 

    if (vmass != 0.):
        hel0 = 1.-abs(nhel) 

        if (pp == 0.):
            vc[2] = complex(0.,0.)
            vc[3] = complex(-nhel*sqh,0.)
            vc[4] = complex(0.,nsvahl*sqh) 
            vc[5] = complex(hel0,0.)

        else:
            emp = p[0]/(vmass*pp)
            vc[2] = complex(hel0*pp/vmass,0.)
            vc[5] = complex(hel0*p[3]*emp+nhel*pt/pp*sqh)
            if (pt != 0.):
                pzpt = p[3]/(pp*pt)*sqh*nhel
                vc[3] = complex(hel0*p[1]*emp-p[1]*pzpt, \
                    -nsvahl*p[2]/pt*sqh)
                vc[4] = complex(hel0*p[2]*emp-p[2]*pzpt, \
                    nsvahl*p[1]/pt*sqh) 
            else:
                vc[3] = complex(-nhel*sqh,0.)
                vc[4] = complex(0.,nsvahl*sign(sqh,p[3]))
    else: 
        pp = p[0]
        pt = sqrt(p[1]**2 + p[2]**2)
        vc[2] = complex(0.,0.)
        vc[5] = complex(nhel*pt/pp*sqh)
        if (pt != 0.):
            pzpt = p[3]/(pp*pt)*sqh*nhel
            vc[3] = complex(-p[1]*pzpt,-nsv*p[2]/pt*sqh)
            vc[4] = complex(-p[2]*pzpt,nsv*p[1]/pt*sqh)
        else:
            vc[3] = complex(-nhel*sqh,0.)
            vc[4] = complex(0.,nsv*sign(sqh,p[3]))
    
    return vc

def sign(x,y):
    """Fortran's sign transfer function"""
    try:
        cmp = (y < 0.)
    except TypeError:
        # y should be complex
        if abs(y.imag) < 1e-6 * abs(y.real):
            y = y.real
        else:
            raise
    finally:
        if (y < 0.):
            return -abs(x) 
        else:
            return abs(x) 

def sxxxxx(p,nss):
    """initialize a scalar wavefunction"""
    
    sc = WaveFunction(1)
    
    sc[2] = complex(1.,0.)
    sc[0] = complex(p[0]*nss,p[3]*nss)
    sc[1] = complex(p[1]*nss,p[2]*nss)
    return sc


def txxxxx(p, tmass, nhel, nst):
    """ initialize a tensor wavefunction"""
    
    tc = WaveFunction(5)
    
    sqh = sqrt(0.5)
    sqs = sqrt(1/6)

    pt2 = p[1]**2 + p[2]**2
    pp = min(p[0],sqrt(pt2+p[3]**2))
    pt = min(pp,sqrt(pt2))

    ft = {}
    ft[(4,0)] = complex(p[0], p[3]) * nst
    ft[(5,0)] = complex(p[1], p[2]) * nst

    if ( nhel >= 0 ): 
        #construct eps+
        ep = [0] * 4
        
        if ( pp == 0 ):
            #ep[0] = 0
            ep[1] = -sqh
            ep[2] = complex(0, nst*sqh)
            #ep[3] = 0
        else:            
            #ep[0] = 0
            ep[3] = pt/pp*sqh
            if (pt != 0):
               pzpt = p[3]/(pp*pt)*sqh
               ep[1] = complex( -p[1]*pzpt , -nst*p[2]/pt*sqh )
               ep[2] = complex( -p[2]*pzpt ,  nst*p[1]/pt*sqh )
            else:
               ep[1] = -sqh 
               ep[2] = complex( 0 , nst*sign(sqh,p[3]) )
            
         
     
    if ( nhel <= 0 ): 
        #construct eps-
        em = [0] * 4
        if ( pp == 0 ):
            #em[0] = 0
            em[1] = sqh 
            em[2] = complex( 0 , nst*sqh )
            #em[3] = 0
        else:
            #em[0] = 0
            em[3] = -pt/pp*sqh
            if pt:
               pzpt = -p[3]/(pp*pt)*sqh
               em[1] = complex( -p[1]*pzpt , -nst*p[2]/pt*sqh )
               em[2] = complex( -p[2]*pzpt ,  nst*p[1]/pt*sqh )
            else:
               em[1] = sqh
               em[2] = complex( 0 , nst*sign(sqh,p[3]) )
            
    
    if ( abs(nhel) <= 1 ):  
        #construct eps0
        e0 = [0] * 4
        if ( pp == 0 ):
            #e0[0] = complex( rZero )
            #e0[1] = complex( rZero )
            #e0[2] = complex( rZero )
            e0[3] = 1
        else:
            emp = p[0]/(tmass*pp)
            e0[0] = pp/tmass 
            e0[3] = p[3]*emp
            if pt:
               e0[1] = p[1]*emp 
               e0[2] = p[2]*emp 
            #else:
            #   e0[1] = complex( rZero )
            #   e0[2] = complex( rZero )

    if nhel == 2:
        for j in range(4):
            for i in range(4):         
                ft[(i,j)] = ep[i]*ep[j]
    elif nhel == -2:
        for j in range(4):
            for i in range(4):         
                ft[(i,j)] = em[i]*em[j]
    elif tmass == 0:
        for j in range(4):
            for i in range(4):         
                ft[(i,j)] = 0
    elif nhel == 1:
        for j in range(4):
            for i in range(4): 
                ft[(i,j)] = sqh*( ep[i]*e0[j] + e0[i]*ep[j] )
    elif nhel == 0:
        for j in range(4):
            for i in range(4):       
                ft[(i,j)] = sqs*( ep[i]*em[j] + em[i]*ep[j] + 2 *e0[i]*e0[j] )
    elif nhel == -1:
        for j in range(4):
            for i in range(4): 
                ft[(i,j)] = sqh*( em[i]*e0[j] + e0[i]*em[j] )

    else:
        raise Exception('invalid helicity TXXXXXX') 



    tc[2] = ft[(0,0)]
    tc[3] = ft[(0,1)]
    tc[4] = ft[(0,2)]
    tc[5] = ft[(0,3)]
    tc[6] = ft[(1,0)]
    tc[7] = ft[(1,1)]
    tc[8] = ft[(1,2)]
    tc[9] = ft[(1,3)]
    tc[10] = ft[(2,0)]
    tc[11] = ft[(2,1)]
    tc[12] = ft[(2,2)]
    tc[13] = ft[(2,3)]
    tc[14] = ft[(3,0)]
    tc[15] = ft[(3,1)]
    tc[16] = ft[(3,2)]
    tc[17] = ft[(3,3)]
    tc[0] = ft[(4,0)]
    tc[1] = ft[(5,0)]

    return tc

def irxxxx(p, mass, nhel, nsr):
    """ initialize a incoming spin 3/2 wavefunction."""
    
    # This routines is a python translation of a routine written by
    # K.Mawatari in fortran dated from the 2008/02/26
    
    ri = WaveFunction(4)
    
    sqh = sqrt(0.5)
    sq2 = sqrt(2)
    sq3 = sqrt(3)
    
    pt2 = p[1]**2 + p[2]**2
    pp = min([p[0], sqrt(pt2+p[3]**2)])
    pt = min([pp,sqrt(pt2)])
    
    rc = {}
    rc[(4,0)] = -1*complex(p[0],p[3])*nsr
    rc[(5,0)] = -1*complex(p[1],p[2])*nsr
    

    nsv = -nsr # nsv=+1 for final, -1 for initial   
        
    # Construct eps+
    if nhel > 0:
        ep = [0] * 4
        if pp:
            #ep[0] = 0
            ep[3] = pt/pp*sqh
            if pt:
                pzpt = p[3]/(pp*pt)*sqh
                ep[1] = complex( -p[1]*pzpt , -nsv*p[2]/pt*sqh )
                ep[2] = complex( -p[2]*pzpt ,  nsv*p[1]/pt*sqh )
            else:
                ep[1] = -sqh 
                ep[2] = complex( 0 , nsv*sign(sqh,p[3]) )
        else:
            #ep[0] = 0d0
            ep[1] = -sqh
            ep[2] = complex(0, nsv * sqh)
            #ep[3] = 0d0

         
    if ( nhel < 0 ): 
        #construct eps-
        em = [0] * 4
        if ( pp == 0 ):
            #em[0] = 0
            em[1] = sqh 
            em[2] = complex( 0 , nsv*sqh )
            #em[3] = 0
        else:
            #em[0] = 0
            em[3] = -pt/pp*sqh
            if pt:
                pzpt = -p[3]/(pp*pt)*sqh
                em[1] = complex( -p[1]*pzpt , -nsv*p[2]/pt*sqh )
                em[2] = complex( -p[2]*pzpt ,  nsv*p[1]/pt*sqh )
            else:
                em[1] = sqh
                em[2] = complex( 0 , nsv*sign(sqh,p[3]) )            
                
    if ( abs(nhel) <= 1 ):  
        #construct eps0
        e0 = [0] * 4
        if ( pp == 0 ):
            #e0[0] = complex( rZero )
            #e0[1] = complex( rZero )
            #e0[2] = complex( rZero )
            e0[3] = 1
        else:
            emp = p[0]/(mass*pp)
            e0[0] = pp/mass 
            e0[3] = p[3]*emp
            if pt:
               e0[1] = p[1]*emp 
               e0[2] = p[2]*emp 
            #else:
            #   e0[1] = complex( rZero )
            #   e0[2] = complex( rZero )

    

    if ( nhel >= -1 ):
        # constract spinor+ 
        fip = [0] * 4
        sf, omega, sfomeg, chi = [0, 0], [0,0], [0,0], [0,0]
        nh = nsr
        if  mass:
            pp = min([p[0],sqrt(p[1]**2+p[2]**2+p[3]**2)])
            if pp == 0:
                sqm = sqrt(mass)
                ip = (1+nh)//2
                im = (1-nh)//2
                fip[0] = ip * sqm
                fip[1] = im * nsr * sqm
                fip[2] = ip * nsr * sqm
                fip[3] = im * sqm
            else:
                sf[0] = float(1+nsr+(1-nsr)*nh)*0.5
                sf[1] = float(1+nsr-(1-nsr)*nh)*0.5
                omega[0] = sqrt(p[0]+pp)
                omega[1] = mass/omega[0]
                ip = ((3+nh)//2) -1 # -1 since they are index 
                im = ((3-nh)//2) -1 # -1 since they are index
                sfomeg[0] = sf[0]*omega[ip]
                sfomeg[1] = sf[1]*omega[im]
                pp3 = max([pp+p[3],0])
                chi[0] = sqrt(pp3*0.5/pp)
                if  pp3 ==0:
                    chi[1] = -nh
                else:
                    chi[1] = complex( nh*p[1] , p[2] )/sqrt(2*pp*pp3)
            
                fip[0] = sfomeg[0]*chi[im]
                fip[1] = sfomeg[0]*chi[ip]
                fip[2] = sfomeg[1]*chi[im]
                fip[3] = sfomeg[1]*chi[ip]
        else:
            sqp0p3 = sqrt(max([p[0]+p[3],0])) * nsr
            chi[0] = sqp0p3
            if  sqp0p3 == 0:
                chi[1] = -nhel *  sqrt(2*p[0])
            else:
                chi[1] = complex( nh*p[1], p[2] )/sqp0p3
            if  nh == 1:
                #fip[0] = complex( rZero )
                #fip[1] = complex( rZero )
                fip[2] = chi[0]
                fip[3] = chi[1]
            else:
                fip[0] = chi[1]
                fip[1] = chi[0]
                #fip(3) = complex( rZero )
                #fip(4) = complex( rZero )
            
    if ( nhel <= 1 ):
        # constract spinor- 
        fim = [0] * 4
        sf, omega, sfomeg, chi = [0, 0], [0,0], [0,0], [0,0]
        nh = -nsr
        if  mass:
            pp = min([p[0],sqrt(p[1]**2+p[2]**2+p[3]**2)])
            if pp == 0:
                sqm = sqrt(mass)
                ip = (1+nh)/2
                im = (1-nh)/2
                fim[0] = ip * sqm
                fim[1] = im * nsr * sqm
                fim[2] = ip * nsr * sqm
                fim[3] = im * sqm
            else:
                sf[0] = float(1+nsr+(1-nsr)*nh)*0.5
                sf[1] = float(1+nsr-(1-nsr)*nh)*0.5
                omega[0] = sqrt(p[0]+pp)
                omega[1] = mass/omega[0]
                ip = (3+nh)//2 -1
                im = (3-nh)//2 -1
                sfomeg[0] = sf[0]*omega[ip]
                sfomeg[1] = sf[1]*omega[im]
                pp3 = max([pp+p[3],0])
                chi[0] = sqrt(pp3*0.5/pp)
                if  pp3 ==0:
                    chi[1] = -nh
                else:
                    chi[1] = complex( nh*p[1] , p[2] )/sqrt(2*pp*pp3)
            
                fim[0] = sfomeg[0]*chi[im]
                fim[1] = sfomeg[0]*chi[ip]
                fim[2] = sfomeg[1]*chi[im]
                fim[3] = sfomeg[1]*chi[ip]
        else:
            sqp0p3 = sqrt(max([p[0]+p[3],0])) * nsr
            chi[0] = sqp0p3
            if  sqp0p3 == 0:
                chi[1] = -nhel *  sqrt(2*p[0])
            else:
                chi[1] = complex( nh*p[1], p[2] )/sqp0p3
            if  nh == 1:
                #fip[0] = complex( rZero )
                #fip[1] = complex( rZero )
                fim[2] = chi[0]
                fim[3] = chi[1]
            else:
                fim[0] = chi[1]
                fim[1] = chi[0]
                #fip(3) = complex( rZero )
                #fip(4) = complex( rZero )        
      
    

    # recurent relation put her for optimization
    cond1  = (pt == 0 and p[3] >= 0)
    cond2  = (pt == 0 and p[3] < 0)
    
    # spin-3/2 fermion wavefunction
    if nhel == 3:
        for i,j in product(list(range(4)), list(range(4))):
            rc[(i, j)] = ep[i] *fip[j]
    
    elif nhel == 1:
        for i,j in product(list(range(4)), list(range(4))):
            if cond1:
                rc[(i,j)] = sq2/sq3*e0[i]*fip[j] +1/sq3*ep[i]*fim[j]
            elif cond2:
                rc[(i,j)] = sq2/sq3*e0[i]*fip[j] -1/sq3*ep[i]*fim[j]
            else:
                rc[(i,j)] = sq2/sq3*e0[i]*fip[j] + \
                                   1/sq3*ep[i]*fim[j] *complex(p[1],nsr*p[2])/pt  
    elif nhel == -1:
        for i,j in product(list(range(4)), list(range(4))):
            if cond1:
                rc[(i,j)] = 1/sq3*em[i]*fip[j] +sq2/sq3*e0[i]*fim[j]
            elif cond2:
                rc[(i,j)] = 1/sq3*em[i]*fip[j] -sq2/sq3*e0[i]*fim[j]
            else:
                rc[(i,j)] = 1/sq3*em[i]*fip[j] + \
                                sq2/sq3*e0[i]*fim[j] *complex(p[1],nsr*p[2])/pt  
    else:
        for i,j in product(list(range(4)), list(range(4))):
            if cond1:
                rc[(i, j)] = em[i] *fim[j]
            elif cond2:
                rc[(i, j)] = -em[i] *fim[j]
            else:
                rc[(i, j)] = em[i]*fim[j] *complex(p[1],nsr*p[2])/pt 
                
    ri[2] = rc[(0,0)]
    ri[3] = rc[(0,1)]
    ri[4] = rc[(0,2)]
    ri[5] = rc[(0,3)]
    ri[6] = rc[(1,0)]
    ri[7] = rc[(1,1)]
    ri[8] = rc[(1,2)]
    ri[9] = rc[(1,3)]
    ri[10] = rc[(2,0)]
    ri[11] = rc[(2,1)]
    ri[12] = rc[(2,2)]
    ri[13] = rc[(2,3)]
    ri[14] = rc[(3,0)]
    ri[15] = rc[(3,1)]
    ri[16] = rc[(3,2)]
    ri[17] = rc[(3,3)]
    ri[0] = rc[(4,0)]
    ri[1] = rc[(5,0)]              

    return ri

def orxxxx(p, mass, nhel, nsr):
    """ initialize a incoming spin 3/2 wavefunction."""
    
    # This routines is a python translation of a routine written by
    # K.Mawatari in fortran dated from the 2008/02/26

   
    ro = WaveFunction(spin=4)
    
    sqh = sqrt(0.5)
    sq2 = sqrt(2)
    sq3 = sqrt(3)
    
    pt2 = p[1]**2 + p[2]**2
    pp = min([p[0], sqrt(pt2+p[3]**2)])
    pt = min([pp,sqrt(pt2)])
    rc = {}
    rc[(4,0)] = complex(p[0],p[3])*nsr
    rc[(5,0)] = complex(p[1],p[2])*nsr
    

    nsv = nsr # nsv=+1 for final, -1 for initial   
        
    # Construct eps+
    if nhel > 0:
        ep = [0] * 4
        if pp:
            #ep[0] = 0
            ep[3] = pt/pp*sqh
            if pt:
                pzpt = p[3]/(pp*pt)*sqh
                ep[1] = complex( -p[1]*pzpt , -nsv*p[2]/pt*sqh )
                ep[2] = complex( -p[2]*pzpt ,  nsv*p[1]/pt*sqh )
            else:
                ep[1] = -sqh 
                ep[2] = complex( 0 , nsv*sign(sqh,p[3]) )
        else:
            #ep[0] = 0d0
            ep[1] = -sqh
            ep[2] = complex(0, nsv * sqh)
            #ep[3] = 0d0
         
    if ( nhel < 0 ): 
        #construct eps-
        em = [0] * 4
        if ( pp == 0 ):
            #em[0] = 0
            em[1] = sqh 
            em[2] = complex( 0 , nsv*sqh )
            #em[3] = 0
        else:
            #em[0] = 0
            em[3] = -pt/pp*sqh
            if pt:
                pzpt = -p[3]/(pp*pt)*sqh
                em[1] = complex( -p[1]*pzpt , -nsv*p[2]/pt*sqh )
                em[2] = complex( -p[2]*pzpt ,  nsv*p[1]/pt*sqh )
            else:
                em[1] = sqh
                em[2] = complex( 0 , nsv*sign(sqh,p[3]) )            
                
    if ( abs(nhel) <= 1 ):  
        #construct eps0
        e0 = [0] * 4
        if ( pp == 0 ):
            #e0[0] = complex( rZero )
            #e0[1] = complex( rZero )
            #e0[2] = complex( rZero )
            e0[3] = 1
        else:
            emp = p[0]/(mass*pp)
            e0[0] = pp/mass 
            e0[3] = p[3]*emp
            if pt:
               e0[1] = p[1]*emp 
               e0[2] = p[2]*emp 
            #else:
            #   e0[1] = complex( rZero )
            #   e0[2] = complex( rZero )

    if nhel >= -1:
        #constract spinor+ 
        nh = nsr
        sqm, fop, omega, sf, sfomeg = [0]*2,[0]*4,[0]*2,[0]*2,[0]*2
        chi = [0]*2
        if mass:
            pp = min([p[0],sqrt(p[1]**2+p[2]**2+p[3]**2)])
            if ( pp == 0):
                sqm[0] = sqrt(abs(mass)) # possibility of negative fermion masses
                sqm[1] = sign(sqm[0],mass) # possibility of negative fermion masses
                ip = -((1+nh)/2)
                im =  (1-nh)/2
                fop[0] = im     * sqm[im]
                fop[1] = ip*nsr * sqm[im]
                fop[2] = im*nsr * sqm[-ip]
                fop[3] = ip     * sqm[-ip]
            else:
                pp = min(p[0],sqrt(p[1]**2+p[2]**2+p[3]**2))
                sf[0] = (1+nsr+(1-nsr)*nh)*0.5
                sf[1] = (1+nsr-(1-nsr)*nh)*0.5
                omega[0] = sqrt(p[0]+pp)
                omega[1] = mass/omega[0]
                ip = (3+nh)//2  -1 # -1 since this is index
                im = (3-nh)//2  -1 # -1 since this is index 
                sfomeg[0] = sf[0]*omega[ip]
                sfomeg[1] = sf[1]*omega[im]
                pp3 = max([pp+p[3],0])
                chi[0] = sqrt(pp3*0.5/pp)
                if pp3 == 0:
                    chi[1] = -nh 
                else:
                    chi[1] = complex( nh*p[1] , -p[2] )/sqrt(2*pp*pp3)

            
                fop[0] = sfomeg[1]*chi[im]
                fop[1] = sfomeg[1]*chi[ip]
                fop[2] = sfomeg[0]*chi[im]
                fop[3] = sfomeg[0]*chi[ip]

        else:
            if(p[1] == 0 and p[2] == 0 and p[3] < 0):
                sqp0p3 = 0
            else:
                sqp0p3 = sqrt(max(p[0]+p[3], 0))*nsr
                
            chi[0] =  sqp0p3
            if ( sqp0p3 == 0 ):
                chi[1] = complex(-nhel )*sqrt(2*p[0])
            else:
                chi[1] = complex( nh*p[1], -p[2] )/sqp0p3
         
            if ( nh == 1 ):
                fop[0] = chi[0]
                fop[1] = chi[1]
                #fop[2] = 0
                #fop[3] = 0
            else:
                #fop[0] = 0
                #fop[1] = 0
                fop[2] = chi[1]
                fop[3] = chi[0]
         
    
    if ( nhel < 2 ):
        # constract spinor+ 
        sqm, fom, omega, sf, sfomeg = [0]*2,[0]*4,[0]*2,[0]*2,[0]*2
        chi = [0]*2

        
        nh = -nsr
        if mass:
            pp = min([p[0],sqrt(p[1]**2+p[2]**2+p[3]**2)])
            if ( pp == 0):
                sqm[0] = sqrt(abs(mass)) # possibility of negative fermion masses
                sqm[1] = sign(sqm[0],mass) # possibility of negative fermion masses
                ip = -((1+nh)/2)
                im =  (1-nh)/2
            
                fom[0] = im     * sqm[im]
                fom[1] = ip*nsr * sqm[im]
                fom[2] = im*nsr * sqm[-ip]
                fom[3] = ip     * sqm[-ip]
            
            else:
                pp = min([p[0],sqrt(p[1]**2+p[2]**2+p[3]**2)])
                sf[0] = (1+nsr+(1-nsr)*nh)*0.5
                sf[1] = (1+nsr-(1-nsr)*nh)*0.5
                omega[0] = sqrt(p[0]+pp)
                omega[1] = mass/omega[0]
                ip = (3+nh)//2 -1 #-1 since ip is an index
                im = (3-nh)//2 -1 
                sfomeg[0] = sf[0]*omega[ip]
                sfomeg[1] = sf[1]*omega[im]
                pp3 = max([pp+p[3], 0])
                chi[0] = sqrt(pp3*0.5/pp)
                if ( pp3 == 0):
                    chi[1] = -nh
                else:
                    chi[1] = complex( nh*p[1] , -p[2] )/sqrt(2*pp*pp3)
            
            
                fom[0] = sfomeg[1]*chi[im]
                fom[1] = sfomeg[1]*chi[ip]
                fom[2] = sfomeg[0]*chi[im]
                fom[3] = sfomeg[0]*chi[ip]
        else:
            if(p[1] == 0 == p[2] and p[3] < 0):
                sqp0p3 = 0
            else:
                sqp0p3 = sqrt(max([p[0]+p[3],0]))*nsr
            chi[0] = sqp0p3
            if ( sqp0p3 == 0):
                chi[1] = complex(-nhel )*sqrt(2*p[0])
            else:
                chi[1] = complex( nh*p[1], -p[2] )/sqp0p3
            if ( nh == 1 ):
                fom[0] = chi[0]
                fom[1] = chi[1]
                #fom[2] = 0
                #fom[3] = 0
            else:
                #fom[1] = 0
                #fom[2] = 0
                fom[2] = chi[1]
                fom[3] = chi[0]

    cond1 = ( pt==0 and p[3]>=0)
    cond2= (pt==0 and p[3]<0)

   
    # spin-3/2 fermion wavefunction
    if nhel == 3:
        for i,j in product(list(range(4)), list(range(4))):
            rc[(i, j)] = ep[i] *fop[j]  
    

    elif nhel == 1:
        for i,j in product(list(range(4)), list(range(4))):
            if cond1:
                rc[(i,j)] = sq2/sq3*e0[i]*fop[j] + 1/sq3*ep[i]*fom[j]
            elif cond2:
                rc[(i,j)] = sq2/sq3*e0[i]*fop[j] - 1/sq3*ep[i]*fom[j]
            else:
                rc[(i,j)] = sq2/sq3*e0[i]*fop[j] + 1/sq3*ep[i]*fom[j] * \
                                                      complex(p[1],-nsr*p[2])/pt  
                
    elif nhel == -1:
        for i,j in product(list(range(4)), list(range(4))):
            if cond1:
                rc[(i,j)] = 1/sq3*em[i]*fop[j]+sq2/sq3*e0[i]*fom[j]
            elif cond2:
                rc[(i,j)] =1/sq3*em[i]*fop[j]-sq2/sq3*e0[i]*fom[j]
            else:
                rc[(i,j)] =  1/sq3*em[i]*fop[j] + sq2/sq3*e0[i]*fom[j] *\
                                                      complex(p[1],-nsr*p[2])/pt              
    else:
        for i,j in product(list(range(4)), list(range(4))):
            if cond1:
                rc[(i,j)] = em[i] * fom[j]
            elif cond2:
                rc[(i,j)] = - em[i] * fom[j]
            else:
                rc[(i,j)] = em[i] * fom[j] * complex(p[1],-nsr*p[2])/pt 



    ro[2] = rc[(0,0)]
    ro[3] = rc[(0,1)]
    ro[4] = rc[(0,2)]
    ro[5] = rc[(0,3)]
    ro[6] = rc[(1,0)]
    ro[7] = rc[(1,1)]
    ro[8] = rc[(1,2)]
    ro[9] = rc[(1,3)]
    ro[10] = rc[(2,0)]
    ro[11] = rc[(2,1)]
    ro[12] = rc[(2,2)]
    ro[13] = rc[(2,3)]
    ro[14] = rc[(3,0)]
    ro[15] = rc[(3,1)]
    ro[16] = rc[(3,2)]
    ro[17] = rc[(3,3)]
    ro[0] = rc[(4,0)]
    ro[1] = rc[(5,0)]
    
    return ro

# %% [markdown]
# ## Toolbox

# %%
def expr_to_string(expr: Expression) -> str:
    """Convert a symbolica expression to string."""
    #return expr.to_canonical_string()
    return expr.format_plain()

# %%
def stepped_execution(tn, hep_lib, max_steps=None, t_delta=0.1):
    i_step = 0
    while True:
        i_step += 1
        print("Performing scalar step    #",i_step)
        import time
        if t_delta is not None:
            time.sleep(t_delta)
        print(tn)
        print(i_step)
        tn.execute(n_steps=1, mode=ExecutionMode.Scalar,library=hep_lib)
        print("DONE")

        i_step += 1
        print("Performing reduction step #",i_step)
    
        print(tn)
        print(i_step)
        tn.execute(n_steps=1, mode=ExecutionMode.Single,library=hep_lib)
        print("DONE")
    
        if max_steps and i_step > max_steps:
            break
        try:
            _ = tn.result_tensor(hep_lib)
        except:
            continue
        break

# %%
def get_numerator(g) -> Expression:
    num = E("1")
    for node in g.get_nodes():
        if node.get_name() not in ["edge", "node"]:
            n_num = node.get("num")
            if n_num:
                num *= Es(n_num)
    for edge in g.get_edges():
        e_num = edge.get("num")
        if e_num:
            num *= Es(e_num)
    
    g_attrs = g.get_attributes()
    if "num" in g_attrs:
        num *= Es(g_attrs["num"])
    if 'overall_factor' in g_attrs:
        num *= Es(g_attrs['overall_factor'])

    return num

def get_propagator_denominators(g) -> Expression:
    internal_nodes = [ n.get_name() for n in g.get_nodes() if not any(marker in n.get_name() for marker in ['graph', 'ext', 'edge', 'node']) ]
    den = E("1")
    for edge in g.get_edges():
        source = edge.get_source().split(':')[0]
        destination = edge.get_destination().split(':')[0]
        attrs = edge.get_attributes()
        if source in internal_nodes and destination in internal_nodes:
            a_den = E(f'gammalooprs::Q({edge.get("id")},spenso::cind(0))^2')
            a_den -= E(f'gammalooprs::Q({edge.get("id")},spenso::cind(1))^2')
            a_den -= E(f'gammalooprs::Q({edge.get("id")},spenso::cind(2))^2')
            a_den -= E(f'gammalooprs::Q({edge.get("id")},spenso::cind(3))^2')
            if 'mass' in attrs:
                a_den -= Es(f'{attrs.get("mass")}^2')
            den *= a_den
    return den

def get_projector(g) -> Expression:
    g_attrs = g.get_attributes()
    projector = None
    if "projector" in g_attrs:
        projector = Es(g_attrs["projector"])
    else:
        projector =  Es(g.get_graph_defaults()[0]["projector"])

    # TMPVH temporary fix to current issue in gammaloop when building external proectors
    # projector = projector.replace(E("gammalooprs::u(2,x__)"),E("gammalooprs::vbar(2,x__)"),repeat=True)

    return projector

def get_color_projector(_g) -> Expression:
    return E("(1/3)*spenso::g(spenso::dind(spenso::cof(3,gammalooprs::hedge(1))),spenso::cof(3,gammalooprs::hedge(2)))")

def get_emr_replacements(g) -> Expression:
    replacements = []
    for edge in g.get_edges():
        replacements.append((E(f"gammalooprs::Q({edge.get("id")},gammalooprs::a___)"),Es(edge.get("lmb_rep"))))
    return replacements

def tn_replace_multiple(tensor_network, replacements):
    res = tensor_network
    for lhs, rhs in replacements:
        res = res.replace(lhs, rhs)
    return res

# %%
def spdot(v1,v2):
    return v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3]

def to_complex(cmplx):
    return E(str(complex(cmplx.real, cmplx.imag))).replace(E("python::j"),E("1𝑖"))

def matrix_apply(mat, v, apply_right=False):
    """
    Compute (M[mu,nu] v[nu]) with Minkowski metric eta=diag(-1,1,1,1).
    v is taken as contravariant (v^0, v^1, v^2, v^3); we first lower it:
      v_0 = -v^0, v_i = v^i.
    mat is a 4x4 nested list representing M[mu,nu].
    """
    if len(mat) != 4 or any(len(row) != 4 for row in mat):
        raise ValueError("mat must be 4x4")
    if len(v) != 4:
        raise ValueError("v must have length 4")

    v_cov = (v[0], -v[1], -v[2], -v[3])  # apply eta_{nu rho} to v^rho
    if apply_right:
        return [sum(mat[nu][mu] * v_cov[nu] for nu in range(4)) for mu in range(4)]
    else:
        return [sum(mat[mu][nu] * v_cov[nu] for nu in range(4)) for mu in range(4)]

def function_map_for_evaluation(ks, hels, debug_pols=False, loop_mom=None, rotate_pols=False, mUV=None):
    
    function_map = {
        E('UFO::ee') : complex(0.30795376724436879,0.0),
        E('UFO::GC_1') : complex(0.0,-0.30795376724436879/3.),
    }

    if loop_mom is not None:
        function_map[E('UFO::G')] = complex(1.2177157847767197,0.0)
        function_map[E('UFO::GC_11')] = complex(0.0,1.2177157847767197)
        function_map[E('spenso::TR')] = complex(0.5,0.0)
        function_map[E('spenso::CF')] = complex(4./3.,0.0)
        function_map[E('vakint::EulerGamma')] = complex(0.577215664901533,0.0)
        function_map[E('𝜋')] = complex(math.pi,0.0)

    if mUV is not None:
        function_map[E('gammalooprs::mUV')] = complex(mUV)

    for ext_id, k in enumerate(ks):
        for lor_i, ki in enumerate(k):
            function_map[E(f"gammalooprs::P({ext_id+1},spenso::cind({lor_i}))")] = complex(ki, 0.0)

    if loop_mom is not None:
        for lor_i, ki in enumerate(loop_mom):            
            function_map[E(f"gammalooprs::K(0,spenso::cind({lor_i}))")] = complex(ki, 0.0)

    # rotation +𝜋/2 about z-axis
    vector_rotation = [
        [1, 0, 0, 0],
        [0, 0, 1, 0],
        [0, -1, 0, 0],
        [0, 0, 0, -1]
    ]
    inv_vector_rotation = [
        [1, 0, 0, 0],
        [0, 0, -1, 0],
        [0, 1, 0, 0],
        [0, 0, 0, -1]
    ]
    spinor_phase = complex(1./sqrt(2.),1./sqrt(2.))
    cmplx_zero = complex(0.0, 0.0)
    spinor_rotation = [
        [spinor_phase.conjugate(), cmplx_zero,cmplx_zero,cmplx_zero],
        [cmplx_zero,spinor_phase,cmplx_zero,cmplx_zero],
        [cmplx_zero,cmplx_zero,spinor_phase.conjugate(),cmplx_zero],
        [cmplx_zero,cmplx_zero,cmplx_zero,spinor_phase]
    ]
    inv_spinor_rotation = [
        [spinor_phase, cmplx_zero,cmplx_zero,cmplx_zero],
        [cmplx_zero,spinor_phase.conjugate(),cmplx_zero,cmplx_zero],
        [cmplx_zero,cmplx_zero,spinor_phase,cmplx_zero],
        [cmplx_zero,cmplx_zero,cmplx_zero,spinor_phase.conjugate()]
    ]

    M2 = sqrt(abs(spdot(ks[2],ks[2])))
    if M2 < 1.0e-5:
        M2 = 0.0
    M3 = sqrt(abs(spdot(ks[3],ks[3])))
    if M3 < 1.0e-5:
        M3 = 0.0
    M4 = sqrt(abs(spdot(ks[4],ks[4])))
    if M4 < 1.0e-5:
        M4 = 0.0

    if rotate_pols:
        vec = matrix_apply(inv_vector_rotation, ks[0])
        pol_1 = ixxxxx(vec,0.0,hels[0],1)[2::]
        pol_1 = matrix_apply(spinor_rotation,pol_1)
    else:
        pol_1 = ixxxxx(ks[0],0.0,hels[0],1)[2::]
    
    if debug_pols: print(pol_1)
    for lor_i, poli in enumerate(pol_1):
        function_map[E(f"gammalooprs::u(1,spenso::cind({lor_i}))")]    = poli

    if rotate_pols:
        vec = matrix_apply(inv_vector_rotation, ks[1])
        pol_2 = oxxxxx(vec,0.0,hels[1],-1)[2::]
        pol_2 = matrix_apply(inv_spinor_rotation,pol_2, apply_right=True)
    else:
        pol_2 = oxxxxx(ks[1],0.0,hels[1],-1)[2::]
    if debug_pols: print(pol_2)
    for lor_i, poli in enumerate(pol_2):    
        function_map[E(f"gammalooprs::vbar(2,spenso::cind({lor_i}))")] = poli
        
    if rotate_pols:
        vec = matrix_apply(inv_vector_rotation, ks[2])
        pol_3 = vxxxxx(vec,M2,hels[2],1)[2::]
        pol_3 = matrix_apply(vector_rotation,pol_3)
    else:
        pol_3 = vxxxxx(ks[2],M2,hels[2],1)[2::]
    if debug_pols: print(pol_3)
    for lor_i, poli in enumerate(pol_3):
        function_map[E(f"gammalooprs::ϵbar(3,spenso::cind({lor_i}))")] = poli
    if rotate_pols:
        vec = matrix_apply(inv_vector_rotation, ks[3])
        pol_4 = vxxxxx(vec,M3,hels[3],1)[2::]
        pol_4 = matrix_apply(vector_rotation,pol_4)
    else:
        pol_4 = vxxxxx(ks[3],M3,hels[3],1)[2::]
    if debug_pols: print(pol_4)
    for lor_i, poli in enumerate(pol_4):
        function_map[E(f"gammalooprs::ϵbar(4,spenso::cind({lor_i}))")] = poli
    if rotate_pols:
        vec = matrix_apply(inv_vector_rotation, ks[4])
        pol_5 = vxxxxx(vec,M4,hels[4],1)[2::]
        pol_5 = matrix_apply(vector_rotation,pol_5)
    else:
        pol_5 = vxxxxx(ks[4],M4,hels[4],1)[2::]
    if debug_pols: print(pol_5)
    for lor_i, poli in enumerate(pol_5):
        function_map[E(f"gammalooprs::ϵbar(0,spenso::cind({lor_i}))")] = poli

    return function_map

# %%
EpsExpansionReplacements = [
    Replacement(E("dim"),E("4 - 2*ε")),
    Replacement(E("𝚪(1-ε)"),E("1 + γₑ*ε + (1/12)*( 6*γₑ^2 + 𝜋^2)*ε^2 + ε^3*O(Gamma,eps^3)")),
    Replacement(E("𝚪(1+ε)"),E("1 - γₑ*ε + (1/12)*( 6*γₑ^2 + 𝜋^2)*ε^2 + ε^3*O(Gamma,eps^3)")),
    Replacement(E("𝚪(1-b_*ε)"),E("1 + γₑ*b_*ε + (1/12)*( 6*γₑ^2 + 𝜋^2)*b_^2*ε^2 + ε^3*O(Gamma,eps^3)")),
    Replacement(E("𝚪(1+b_*ε)"),E("1 - γₑ*b_*ε + (1/12)*( 6*γₑ^2 + 𝜋^2)*b_^2*ε^2 + ε^3*O(Gamma,eps^3)")),
    Replacement(E("𝚪(ε)"),E("1/ε - γₑ + (1/12)*( 6*γₑ^2 + 𝜋^2)*ε + ε^2*O(Gamma,eps^2)")),
    Replacement(E("𝚪(b_*ε)"),E("1/(b_*ε) - γₑ + (1/12)*( 6*γₑ^2 + 𝜋^2)*b_*ε + ε^2*O(Gamma,eps^2)")),
    Replacement(E("𝚪(2-ε)"),E("1 + (γₑ-1)*ε + (1/12)*( -12*γₑ + 6 * γₑ^2 + 𝜋^2)*ε^2 + ε^3*O(Gamma,eps^3)")),
    Replacement(E("𝚪(2+ε)"),E("1 + (1-γₑ)*ε + (1/12)*( -12*γₑ + 6 * γₑ^2 + 𝜋^2)*ε^2 + ε^3*O(Gamma,eps^3)")),
    Replacement(E("𝚪(2-b_*ε)"),E("1 + (γₑ-1)*b_*ε + (1/12)*( -12*γₑ + 6 * γₑ^2 + 𝜋^2)*b_^2*ε^2 + ε^3*O(Gamma,eps^3)")),
    Replacement(E("𝚪(2+b_*ε)"),E("1 + (1-γₑ)*b_*ε + (1/12)*( -12*γₑ + 6 * γₑ^2 + 𝜋^2)*b_^2*ε^2 + ε^3*O(Gamma,eps^3)")),
]

# %%
def eps_expansion_finite(expr, coeff_index=-1):
    expansion = expr.replace_multiple(EpsExpansionReplacements).series(E("ε"),0,0,depth_is_absolute=True).to_expression().coefficient_list(E("ε"))
    if coeff_index is None:
        return expansion
    else:
        return expansion[coeff_index][-1]

# %%
def prepare_integrated_ct_numerator_expression(expr):
    processed_expr = expr.replace(E("mUVsq"),E("gammalooprs::mUV^2"))
    processed_expr = processed_expr.replace(E("γₑ"),E("vakint::EulerGamma"))
    return processed_expr

# %% [markdown]
# ## Kinematics

# %%
kin_point = [
    [1.0, 0.0, 0.0, 1.0],
    [1.0, 0.0, 0.0, -1.0],
    [0.5, 0.3, 0.0, 0.0],
    [0.2, 0.0, 0.1, 0.1],
]
kin_point.append(
    [sum(k[i] for k in kin_point[:2])-sum(k[i] for k in kin_point[2:]) for i in range(4)]
)
helicities=[ 1, -1, -1, -1, 1 ]

kin_point2 = [
    [1.0, 0.0, 0.0, 1.0],
    [1.0, 0.0, 0.0, -1.0],
    [0.4, 0.1, 0.0, 0.0],
    [0.7, 0.0, 0.2, 0.6],
]
kin_point2.append(
    [sum(k[i] for k in kin_point2[:2])-sum(k[i] for k in kin_point2[2:]) for i in range(4)]
)
helicities2=[ 1, -1, -1, -1, 1 ]

kin_point_onshell = [
    [0.1000000000000000e+1,   0.0000000000000000e+0,   0.0000000000000000e+00,   0.1000000000000000e+01],
    [0.1000000000000000e+1,   0.0000000000000000e+0,   0.0000000000000000e+00,  -0.1000000000000000e+01],
    [0.9171575757708805e+0,   0.3389064406193597e+0,   0.7593073241563973e+00,  -0.3870049493005049e+00],
    [0.7281332414736353e+0,  -0.3665973858638370e-1,  -0.6954086026387343e+00,   0.2126992155174163e+00]
]
kin_point_onshell.append(
    [sum(k[i] for k in kin_point_onshell[:2])-sum(k[i] for k in kin_point_onshell[2:]) for i in range(4)]
)
helicities_onshell=[ 1, -1, -1, -1, 1 ]

kin_point2_onshell = [
    [0.1000000000000000e+1,   0.0000000000000000e+0,   0.0000000000000000e+0,   0.1000000000000000e+1],
    [0.1000000000000000e+1,   0.0000000000000000e+0,   0.0000000000000000e+0,  -0.1000000000000000e+1],
    [0.9903067841547669e+0,   0.3734458315384239e+0,   0.6393671561700484e+0,  -0.6576133949826120e+0],
    [0.2053558277500182e+0,  -0.1412408546154445e+0,  -0.8691190591393230e-1,   0.1211129951276980e+0]
]
kin_point2_onshell.append(
    [sum(k[i] for k in kin_point2_onshell[:2])-sum(k[i] for k in kin_point2_onshell[2:]) for i in range(4)]
)
helicities2_onshell=[ 1, -1, -1, -1, 1 ]

# rotation +𝜋/2 about z-axis
kin_point2_onshell_rotated = [
    [0.1000000000000000e+1,   0.0000000000000000e+0,   0.0000000000000000e+0,   0.1000000000000000e+1],
    [0.1000000000000000e+1,   0.0000000000000000e+0,   0.0000000000000000e+0,  -0.1000000000000000e+1],
    [0.9903067841547669e+0,   0.3734458315384239e+0,   0.6393671561700484e+0,  -0.6576133949826120e+0],
    [0.2053558277500182e+0,  -0.1412408546154445e+0,  -0.8691190591393230e-1,   0.1211129951276980e+0]
]
kin_point2_onshell_rotated.append(
    [sum(k[i] for k in kin_point2_onshell_rotated[:2])-sum(k[i] for k in kin_point2_onshell_rotated[2:]) for i in range(4)]
)
# Apply a pi2 rotation around the z-axis
kin_point2_onshell_rotated = [
    [k[0],-k[2],k[1],k[3]]
    for k in kin_point2_onshell_rotated
]
helicities2_onshell_rotated=[ 1, -1, -1, -1, 1 ]

def display_sample_point(kins, hels=None, name=None):
    if hels is None:
        hels = [None] * len(kins)

    print(f"Sample point{'' if name is None else ' '+name}:")
    print('-'*(148 if hels[0] is not None else 139))
    for i, (k, h) in enumerate(zip(kins, hels)):
        line = f"#{i}    {k[0]:-25.16e} {k[1]:-25.16e} {k[2]:-25.16e} {k[3]:-25.16e} | m={sqrt(abs(spdot(k,k))):-25.16e}"
        if h is not None:
            line += f" | hel={h:2d}"
        print(line)
    k_sum_in = [sum(k[i] for k in kins[:2]) for i in range(4)]
    k_sum_out = [sum(k[i] for k in kins[2:]) for i in range(4)]
    k_sum = [k_sum_in[i] - k_sum_out[i] for i in range(4)]
    print('-'*(148 if hels[0] is not None else 139))
    print(f"Sum:  {k_sum[0]:-25.16e} {k_sum[1]:-25.16e} {k_sum[2]:-25.16e} {k_sum[3]:-25.16e}")

# display_sample_point(kin_point, helicities, "PS point 1")
# print('')
# display_sample_point(kin_point2, helicities2, "PS point 2")
# print('')
#display_sample_point(kin_point_onshell, helicities_onshell, ' "onshell 1"')
#print('')
#display_sample_point(kin_point2_onshell, helicities2_onshell, ' "onshell 2"')
#print('')
#display_sample_point(kin_point2_onshell_rotated, helicities2_onshell_rotated, ' "rotated onshell 2"')

# This corresponds then to the ordering (p1,p3,p4,p5,p2) when reading from the quark line from d to d~, which is diagram #2 in MadGraph ordering
# Targets as obtained from /Users/vjhirsch/MG5/MG5_aMC_v3_5_8/DDX_AAA_from_loops/SubProcesses/P0_ddx_aaa
MGTarget_kinpoin1 = complex(1.78588223425167137e-3,-1.16298831119980095e-3)
MGTarget_kinpoin2 = complex(7.67583136957218854e-4, 6.61541093750165753e-5)
MGTarget_kinpoin1_onshell = complex(-6.19639559102668088e-4,8.41535549668894736e-4)
MGTarget_kinpoin2_onshell = complex(1.47276041641056167e-4,-1.15031393691302202e-3)
MGTarget_kinpoin2_onshell_rotated = complex(1.47276041641056167e-4,-1.15031393691302202e-3)

# %%
# Target result from MadLoop

MadLoopTargets = { 
    'onshell_point_1': { 'born': -6.1963955910266809e-04+8.4153554966889474e-04j,
                       'loop_amplitude': ( +9.3441903738697950e-04-3.2651826076211311e-04j, +1.5600762604558023e-04-1.0788527256794146e-04j, +1.1637006454444246e-05-1.5804276017048061e-05j ),
                       'loop_graphs': { 'qqx_aaa_box_A': { 'loop': ( +4.5804809636893807e-04+1.8047684780568232e-04j, +4.9024019250802263e-05+1.3364186853142279e-05j, -1.2186249490916134e-20-2.5370334496306547e-20j ),
                                                           'loop_no_r2': ( +4.5804809636893807e-04+1.8047684780568232e-04j, +4.9024019250802263e-05+1.3364186853142279e-05j, -1.2186249490916134e-20-2.5370334496306547e-20j ),
                                                           'r2': ( +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                        'qqx_aaa_box_B': { 'loop': ( +2.5095303182480980e-06-8.0722834798604608e-07j, +5.4473937339358428e-07-7.5599987780869753e-07j, +3.4191636510665646e-19-1.0641975219183083e-20j ),
                                                           'loop_no_r2': ( +2.5095303182480980e-06-8.0722834798604608e-07j, +5.4473937339358428e-07-7.5599987780869753e-07j, +3.4191636510665646e-19-1.0641975219183083e-20j ),
                                                           'r2': ( +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                        'qqx_aaa_bub_A': { 'loop': ( +5.0817672657145043e-05-6.9015732556447051e-05j, +5.8185032272218636e-06-7.9021380085242286e-06j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                           'loop_no_r2': ( +5.8575676960107526e-05-7.9551916567812694e-05j, +5.8185032272218636e-06-7.9021380085242286e-06j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                           'r2': ( -7.7580043029624842e-06+1.0536184011365638e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                        'qqx_aaa_bub_B': { 'loop': ( +5.6067706614605624e-05-7.6145829638321730e-05j, +5.8185032272218585e-06-7.9021380085242354e-06j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                           'loop_no_r2': ( +6.3825710917568107e-05-8.6682013649687373e-05j, +5.8185032272218585e-06-7.9021380085242354e-06j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                           'r2': ( -7.7580043029624808e-06+1.0536184011365638e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                        'qqx_aaa_pentagon': { 'loop': ( +5.9128949395222581e-05-3.7896897696869895e-04j, +7.0867876229831813e-05-1.1095348596919195e-04j, +1.1637006454443949e-05-1.5804276017048013e-05j ),
                                                              'loop_no_r2': ( +5.9128949395222581e-05-3.7896897696869895e-04j, +7.0867876229831813e-05-1.1095348596919195e-04j, +1.1637006454443949e-05-1.5804276017048013e-05j ),
                                                              'r2': ( +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                        'qqx_aaa_tri_A': { 'loop': ( +1.6864786944174371e-04-1.5899343657034765e-04j, +1.4766638547489225e-05-1.3425785950271334e-05j, -1.1290788709922817e-20+9.9463607474501254e-21j ),
                                                           'loop_no_r2': ( +1.5313186083581876e-04-1.3792106854761636e-04j, +1.4766638547489225e-05-1.3425785950271334e-05j, -1.1290788709922817e-20+9.9463607474501254e-21j ),
                                                           'r2': ( +1.5516008605924962e-05-2.1072368022731275e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                        'qqx_aaa_tri_B': { 'loop': ( -4.5286816488248621e-05+5.8856227407092242e-05j, -5.8185032272218610e-06+7.9021380085241016e-06j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                           'loop_no_r2': ( -6.0802825094173586e-05+7.9928595429823514e-05j, -5.8185032272218610e-06+7.9021380085241016e-06j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                           'r2': ( +1.5516008605924968e-05-2.1072368022731275e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                        'qqx_aaa_tri_C': { 'loop': ( +1.8448602907932503e-04+1.1807986810691383e-04j, +1.4985849416841507e-05+1.1787950384712611e-05j, -2.2689425663765589e-20-2.4897792632525507e-20j ),
                                                           'loop_no_r2': ( +1.6897002047340008e-04+1.3915223612964511e-04j, +1.4985849416841507e-05+1.1787950384712611e-05j, -2.2689425663765589e-20-2.4897792632525507e-20j ),
                                                           'r2': ( +1.5516008605924962e-05-2.1072368022731268e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )}}},
  'onshell_point_2': { 'born': +1.4727604164105617e-04-1.1503139369130220e-03j,
                       'loop_amplitude': ( -6.4224166840662377e-04+7.0413313456657925e-04j, -9.3147479931491922e-05+1.8875448630951172e-04j, -2.7658857830895706e-06+2.1603221601729886e-05j ),
                       'loop_graphs': { 'qqx_aaa_box_A': { 'loop': ( +4.3501666653595793e-05-5.3964229087703560e-05j, -2.0924449741914300e-07-2.9478686116650080e-06j, -5.2439418102948995e-20-2.5030196455125011e-20j ),
                                                           'loop_no_r2': ( +4.3501666653595793e-05-5.3964229087703560e-05j, -2.0924449741914300e-07-2.9478686116650080e-06j, -5.2439418102948995e-20-2.5030196455125011e-20j ),
                                                           'r2': ( +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                        'qqx_aaa_box_B': { 'loop': ( +3.9448460367809059e-05-4.3611872413422726e-05j, -1.2110960261350754e-07-9.3059904478270586e-08j, -8.3488012042107282e-19+7.9093906145154267e-19j ),
                                                           'loop_no_r2': ( +3.9448460367809059e-05-4.3611872413422726e-05j, -1.2110960261350754e-07-9.3059904478270586e-08j, -8.3488012042107282e-19+7.9093906145154267e-19j ),
                                                           'r2': ( +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                        'qqx_aaa_bub_A': { 'loop': ( -1.1754810054337015e-05+9.1812094354250658e-05j, -1.3829428915447305e-06+1.0801610800865160e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                           'loop_no_r2': ( -1.3598733909729990e-05+1.0621424208873755e-04j, -1.3829428915447305e-06+1.0801610800865160e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                           'r2': ( +1.8439238553929744e-06-1.4402147734486888e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                        'qqx_aaa_bub_B': { 'loop': ( -1.2039999654891499e-05+9.4039595640440967e-05j, -1.3829428915447326e-06+1.0801610800865163e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                           'loop_no_r2': ( -1.3883923510284476e-05+1.0844174337492784e-04j, -1.3829428915447326e-06+1.0801610800865163e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                           'r2': ( +1.8439238553929772e-06-1.4402147734486879e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                        'qqx_aaa_pentagon': { 'loop': ( -7.2767239131245711e-04+3.2550654135640023e-04j, -9.3698749587332969e-05+1.4913267100726384e-04j, -2.7658857830886825e-06+2.1603221601729161e-05j ),
                                                              'loop_no_r2': ( -7.2767239131245711e-04+3.2550654135640023e-04j, -9.3698749587332969e-05+1.4913267100726384e-04j, -2.7658857830886825e-06+2.1603221601729161e-05j ),
                                                              'r2': ( +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                        'qqx_aaa_tri_A': { 'loop': ( +1.2283616393238692e-06+1.8358075150912861e-04j, +4.0885145511176689e-07+1.5501247153951408e-05j, +1.0084704905163540e-21-3.9663093827305667e-21j ),
                                                           'loop_no_r2': ( +4.9162093501098194e-06+1.5477645604015485e-04j, +4.0885145511176689e-07+1.5501247153951408e-05j, +1.0084704905163540e-21-3.9663093827305667e-21j ),
                                                           'r2': ( -3.6878477107859502e-06+2.8804295468973768e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                        'qqx_aaa_tri_B': { 'loop': ( +8.4852866050047818e-06-8.9056696539317959e-05j, +1.3829428915439362e-06-1.0801610800864377e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                           'loop_no_r2': ( +1.2173134315790734e-05-1.1786099200829171e-04j, +1.3829428915439362e-06-1.0801610800864377e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                           'r2': ( -3.6878477107859523e-06+2.8804295468973755e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                        'qqx_aaa_tri_C': { 'loop': ( +1.6561757349328358e-05+1.9582694974680311e-04j, +1.8557151923074678e-06+1.6359885863573822e-05j, -1.6729854014464920e-21-3.8605033870764084e-20j ),
                                                           'loop_no_r2': ( +2.0249605060114308e-05+1.6702265427782936e-04j, +1.8557151923074678e-06+1.6359885863573822e-05j, -1.6729854014464920e-21-3.8605033870764084e-20j ),
                                                           'r2': ( -3.6878477107859502e-06+2.8804295468973762e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )}}},
  'point_1': { 'born': +1.7858822342516714e-03-1.1629883111998010e-03j,
               'loop_amplitude': ( -1.8198810659122121e-03-1.7869260147041730e-04j, -3.7515117422041994e-04+9.4252193393945936e-05j, -3.3539374272615706e-05+2.1841249941294190e-05j ),
               'loop_graphs': { 'qqx_aaa_box_A': { 'loop': ( -1.1270361743962122e-04-2.8279759755650742e-04j, -1.6856938837405495e-05-2.3416433486043445e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'loop_no_r2': ( -1.1270361743962122e-04-2.8279759755650742e-04j, -1.6856938837405495e-05-2.3416433486043445e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'r2': ( +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                'qqx_aaa_box_B': { 'loop': ( +4.0939928435764259e-04-3.5464221060181887e-04j, +3.3716667560516563e-05-3.6405539703372482e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'loop_no_r2': ( +4.0939928435764259e-04-3.5464221060181887e-04j, +3.3716667560516563e-05-3.6405539703372482e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'r2': ( +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                'qqx_aaa_bub_A': { 'loop': ( -1.6546426711789496e-04+1.0775235056861832e-04j, -1.6769687136307955e-05+1.0920624970647086e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'loop_no_r2': ( -1.8782384996630556e-04+1.2231318386281444e-04j, -1.6769687136307955e-05+1.0920624970647086e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'r2': ( +2.2359582848410602e-05-1.4560833294196113e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                'qqx_aaa_bub_B': { 'loop': ( -1.6586837551862686e-04+1.0801551089212402e-04j, -1.6769687136307955e-05+1.0920624970647078e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'loop_no_r2': ( -1.8822795836703746e-04+1.2257634418632012e-04j, -1.6769687136307955e-05+1.0920624970647078e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'r2': ( +2.2359582848410595e-05-1.4560833294196106e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                'qqx_aaa_pentagon': { 'loop': ( -1.2651108691430855e-03+5.0973333727282106e-04j, -3.1321692256751127e-04+1.6324315873432260e-04j, -3.3539374272615706e-05+2.1841249941294190e-05j ),
                                                      'loop_no_r2': ( -1.2651108691430855e-03+5.0973333727282106e-04j, -3.1321692256751127e-04+1.6324315873432260e-04j, -3.3539374272615706e-05+2.1841249941294190e-05j ),
                                                      'r2': ( +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                'qqx_aaa_tri_A': { 'loop': ( -1.4720169734556542e-04+3.5089836818087228e-04j, -9.0546320391435903e-06+2.9177426832762476e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'loop_no_r2': ( -1.0248253164874425e-04+3.2177670159248004e-04j, -9.0546320391435903e-06+2.9177426832762476e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'r2': ( -4.4719165696821190e-05+2.9121666588392216e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                'qqx_aaa_tri_B': { 'loop': ( +1.3730858133326711e-04-9.4857599105121907e-05j, +1.6769687136307606e-05-1.0920624970647232e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'loop_no_r2': ( +1.8202774703008831e-04-1.2397926569351412e-04j, +1.6769687136307606e-05-1.0920624970647232e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'r2': ( -4.4719165696821190e-05+2.9121666588392216e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                'qqx_aaa_tri_C': { 'loop': ( -5.1024010503832810e-04-5.2279476112140476e-04j, -5.2969661200567863e-05-4.9267043954370113e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'loop_no_r2': ( -4.6552093934150690e-04-5.5191642770979700e-04j, -5.2969661200567863e-05-4.9267043954370113e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'r2': ( -4.4719165696821183e-05+2.9121666588392212e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )}}},
  'point_2': { 'born': +7.6758313695721885e-04+6.6154109375016575e-05j,
               'loop_amplitude': ( -4.6454041054975819e-04-4.6649000332111597e-04j, -1.2784747198788324e-04-5.6642320047823324e-05j, -1.4415428756731519e-05-1.2423929145798689e-06j ),
               'loop_graphs': { 'qqx_aaa_box_A': { 'loop': ( +2.7419836240132268e-04-3.8442481267150623e-04j, +2.0833025996963404e-05-3.5538638114773706e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'loop_no_r2': ( +2.7419836240132268e-04-3.8442481267150623e-04j, +2.0833025996963404e-05-3.5538638114773706e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'r2': ( +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                'qqx_aaa_box_B': { 'loop': ( +1.7182397313821466e-04+5.8877561492207663e-06j, +1.4305178929600214e-05-7.9236889742559199e-07j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'loop_no_r2': ( +1.7182397313821466e-04+5.8877561492207663e-06j, +1.4305178929600214e-05-7.9236889742559199e-07j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'r2': ( +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                'qqx_aaa_bub_A': { 'loop': ( -7.2965827589493250e-05-6.2885557362927784e-06j, -7.2077143783627553e-06-6.2119645725966964e-07j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'loop_no_r2': ( -8.2576113427310259e-05-7.1168176793056727e-06j, -7.2077143783627553e-06-6.2119645725966964e-07j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'r2': ( +9.6102858378170059e-06+8.2826194301289448e-07j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                'qqx_aaa_bub_B': { 'loop': ( -8.1461236153753980e-05-7.0207320443516230e-06j, -7.2077143783627528e-06-6.2119645725966996e-07j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'loop_no_r2': ( -9.1071521991570989e-05-7.8489939873645173e-06j, -7.2077143783627528e-06-6.2119645725966996e-07j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'r2': ( +9.6102858378170059e-06+8.2826194301289448e-07j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                'qqx_aaa_pentagon': { 'loop': ( -6.9427037966633047e-04+1.4904202967785512e-04j, -1.4480229772107514e-04+8.0394431584783145e-07j, -1.4415428756731519e-05-1.2423929145798689e-06j ),
                                                      'loop_no_r2': ( -6.9427037966633047e-04+1.4904202967785512e-04j, -1.4480229772107514e-04+8.0394431584783145e-07j, -1.4415428756731519e-05-1.2423929145798689e-06j ),
                                                      'r2': ( +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                'qqx_aaa_tri_A': { 'loop': ( -1.4057746423562403e-04-4.6482650236712811e-06j, -1.1072050083994545e-05+1.6095251199582500e-07j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'loop_no_r2': ( -1.2135689255999002e-04-2.9917411376454909e-06j, -1.1072050083994545e-05+1.6095251199582500e-07j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'r2': ( -1.9220571675634015e-05-1.6565238860257904e-06j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                'qqx_aaa_tri_B': { 'loop': ( +6.8868625392449733e-05+1.1286719248963578e-05j, +7.2077143783627375e-06+6.2119645725961988e-07j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'loop_no_r2': ( +8.8089197068083751e-05+1.2943243134989368e-05j, +7.2077143783627375e-06+6.2119645725961988e-07j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'r2': ( -1.9220571675634015e-05-1.6565238860257896e-06j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )},
                                'qqx_aaa_tri_C': { 'loop': ( +9.8435361634564845e-06-2.3032414292133350e-04j, +9.6385268985607185e-08-2.0655013406207972e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'loop_no_r2': ( +2.9064107839090500e-05-2.2866761903530771e-04j, +9.6385268985607185e-08-2.0655013406207972e-05j, +0.0000000000000000e+00+0.0000000000000000e+00j ),
                                                   'r2': ( -1.9220571675634015e-05-1.6565238860257892e-06j, +0.0000000000000000e+00+0.0000000000000000e+00j, +0.0000000000000000e+00+0.0000000000000000e+00j )}}}}

print('')
display_sample_point(kin_point, helicities, "PS point 1")
print("Tree amplitude target         : {:+.16e}".format(MadLoopTargets['point_1']['born']))
print("Loop amplitude finite target  : {:+.16e}".format(MadLoopTargets['point_1']['loop_amplitude'][0]))
print("Loop amplitude 1/eps  target  : {:+.16e}".format(MadLoopTargets['point_1']['loop_amplitude'][1]))
print("Loop amplitude 1/eps^2 target : {:+.16e}".format(MadLoopTargets['point_1']['loop_amplitude'][2]))
print('')
display_sample_point(kin_point2, helicities2, "PS point 2")
print("Tree amplitude target         : {:+.16e}".format(MadLoopTargets['point_2']['born']))
print("Loop amplitude finite target  : {:+.16e}".format(MadLoopTargets['point_2']['loop_amplitude'][0]))
print("Loop amplitude 1/eps  target  : {:+.16e}".format(MadLoopTargets['point_2']['loop_amplitude'][1]))
print("Loop amplitude 1/eps^2 target : {:+.16e}".format(MadLoopTargets['point_2']['loop_amplitude'][2]))
print('')
display_sample_point(kin_point_onshell, helicities_onshell, "PS onshell point 1")
print("Tree amplitude target         : {:+.16e}".format(MadLoopTargets['onshell_point_1']['born']))
print("Loop amplitude finite target  : {:+.16e}".format(MadLoopTargets['onshell_point_1']['loop_amplitude'][0]))
print("Loop amplitude 1/eps  target  : {:+.16e}".format(MadLoopTargets['onshell_point_1']['loop_amplitude'][1]))
print("Loop amplitude 1/eps^2 target : {:+.16e}".format(MadLoopTargets['onshell_point_1']['loop_amplitude'][2]))
print('')
display_sample_point(kin_point2_onshell, helicities2_onshell, "PS onshell point 2")
print("Tree amplitude target         : {:+.16e}".format(MadLoopTargets['onshell_point_2']['born']))
print("Loop amplitude finite target  : {:+.16e}".format(MadLoopTargets['onshell_point_2']['loop_amplitude'][0]))
print("Loop amplitude 1/eps  target  : {:+.16e}".format(MadLoopTargets['onshell_point_2']['loop_amplitude'][1]))
print("Loop amplitude 1/eps^2 target : {:+.16e}".format(MadLoopTargets['onshell_point_2']['loop_amplitude'][2]))

# %%

# miaou = function_map_for_evaluation(kin_point2_onshell, helicities2_onshell)

# emr_replacements = get_emr_replacements(tree_qqx_aaa_graphs[0])
# for i, (lhs, rhs) in enumerate(emr_replacements):
#     print(f"{E(f"(gammalooprs::Q({i},spenso::cind(1))^2+gammalooprs::Q({i},spenso::cind(2))^2+gammalooprs::Q({i},spenso::cind(3))^2)^(1/2)").to_canonical_string()} -> {to_complex(E(f"(gammalooprs::Q({i},spenso::cind(1))^2+gammalooprs::Q({i},spenso::cind(2))^2+gammalooprs::Q({i},spenso::cind(3))^2)^(1/2)").replace(lhs, rhs).evaluate_complex(miaou,{})).to_canonical_string()}")

# %% [markdown]
# ## Tree-level ME evaluations

# %%
run_gammaloop_commands([
    f'import model {SM_MODEL_PATH}',
    f"import graphs {os.path.join(ROOT_DIR,'dot_graphs','qqx_aaa_tree_unprocessed.dot')}",
    f"save dot {os.path.join(ROOT_DIR,'TMP')}",
    f"!cp {os.path.join(ROOT_DIR,'TMP','processes','amplitudes','qqx_aaa_tree_unprocessed','default','qqx_aaa_tree_1.dot')} {os.path.join(ROOT_DIR,'dot_graphs','qqx_aaa_tree.dot')}",
    f"!rm -rf {os.path.join(ROOT_DIR,'TMP_state')}",
    f"!rm -rf {os.path.join(ROOT_DIR,'TMP')}"
],debug=False)

# %%
tree_qqx_aaa_graphs = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "qqx_aaa_tree.dot"))

# %%
#print(get_numerator(tree_qqx_aaa_graphs[0]).to_canonical_string())

# %% [markdown]
# Extract the full expression of the diagram:

# %%
tree_graph_0_expr = simplify_color(
    cook_indices(
        get_numerator(tree_qqx_aaa_graphs[0])
        *
        get_color_projector(tree_qqx_aaa_graphs[0])
    )
)
#tree_graph_0_expr= simplify_gamma(tree_graph_0_expr)
tree_graph_0_expr *= cook_indices(get_projector(tree_qqx_aaa_graphs[0]))
tree_graph_0_expr /= get_propagator_denominators(tree_qqx_aaa_graphs[0])

# Fix current spenso issue
#tree_graph_0_expr= tree_graph_0_expr.replace(Es("spenso::gamma(a_,b_,c_)"),Es("spenso::gamma(b_,a_,c_)")) 

# %%
#tree_graph_0_expr

# %%
# Experiment: compute a single wavefunction
# WF1expr= cook_indices(Es("""1
# *  (-1/3)*ee   
# *  u(1,spenso::bis(4,hedge(1)))
# *  spenso::gamma(spenso::bis(4,hedge_1),spenso::bis(4,hedge_5),spenso::mink(4,hedge_3))
# *  spenso::gamma(spenso::bis(4,hedge_5),spenso::bis(4,hedge_6),spenso::mink(4,edge_5_1))
# *  ϵbar(3,spenso::mink(4,hedge(3)))
# *  Q(5,spenso::mink(4,edge_5_1))
# *  (Q(5,spenso::cind(0))^2-Q(5,spenso::cind(1))^2-Q(5,spenso::cind(2))^2-Q(5,spenso::cind(3))^2)^-1
# """)).replace(Es("spenso::gamma(a_,b_,c_)"),Es("spenso::gamma(b_,a_,c_)"))
# WF1_tn = TensorNetwork(WF1expr, hep_lib)


# emr_replacements: Expression = get_emr_replacements(tree_qqx_aaa_graphs[0])
# WF1_tn = tn_replace_multiple(WF1_tn, emr_replacements)
# eval_function_map = function_map_for_evaluation(kin_point_onshell, helicities_onshell, debug_pols=True)
# WF1_tn = tn_replace_multiple(WF1_tn, [(lhs, to_complex(rhs)) for lhs, rhs in eval_function_map.items()] )
# WF1_tn.execute(hep_lib)
# wf_expr_tensor = WF1_tn.result_tensor()
# print(wf_expr_tensor)


# print(Es("Q(5,spenso::cind(3))").replace_multiple([Replacement(lhs, rhs) for lhs, rhs in emr_replacements]).evaluate_complex(eval_function_map,{}))



# %% [markdown]
# Prepare tensor networks for evaluation:

# %%
# The following fails execution because the mismatch of the namespace of the colour symbols, but it should not fail!
#test_expression = E("spenso::{symmetric,scalar,real}::g(spenso::{}::cof(3,gammalooprs::{}::hedge_2),spenso::{}::dind(spenso::{}::cof(3,gammalooprs::{}::hedge_1)))*spenso::{symmetric,scalar,real}::g(spenso::{}::cof(3,gammalooprs::{}::hedge_1),spenso::{}::dind(spenso::{}::cof(3,gammalooprs::{}::hedge_2)))")
#stepped_execution(tnnew, hep_lib)
#tn.execute(hep_lib)

# %%
#tn = TensorNetwork(tree_graph_0_expr, hep_lib)
#print(tn)
#stepped_execution(tn, hep_lib)
##tn.execute(hep_lib)
#print(tn)

# %%
#tree_graph_0_expr

# %%
tn = TensorNetwork(tree_graph_0_expr, hep_lib)
tn.execute(hep_lib)
hel_amp_expr = tn.result_scalar()
emr_replacements = get_emr_replacements(tree_qqx_aaa_graphs[0])
hel_amp_expr_lmb = hel_amp_expr.replace_multiple([Replacement(lhs, rhs) for lhs, rhs in emr_replacements])
tn_emr = TensorNetwork(tree_graph_0_expr, hep_lib)
tn_lmb = tn_replace_multiple(tn_emr, emr_replacements)

# %% [markdown]
# Select test points:

# %%
# function_map_for_evaluation(kin_point, helicities)

test_point_A = {'kin': kin_point_onshell, 'hel' : helicities_onshell, 'MG_target': MGTarget_kinpoin1_onshell}
test_point_B = {'kin' : kin_point2_onshell, 'hel' : helicities2_onshell, 'MG_target': MGTarget_kinpoin2_onshell}


# %%
#hel_amp_expr_lmb

# %% [markdown]
# Now run the evaluations:

# %%
eval_function_map = function_map_for_evaluation(kin_point2_onshell_rotated, helicities2_onshell_rotated, rotate_pols=True)
#atom_eval = hel_amp_expr_lmb.replace_multiple([Replacement(lhs, to_complex(rhs)) for lhs, rhs in eval_function_map.items()])
#pprint(eval_function_map)
#print('\n*'.join(hel_amp_expr_lmb.to_canonical_string().split('*')))
atom_eval=to_complex(hel_amp_expr_lmb.evaluate_complex(eval_function_map,{}))
print("Atom eval: ", atom_eval)

# %%
eval_function_map = function_map_for_evaluation(test_point_B['kin'], test_point_B['hel'])
#atom_eval = hel_amp_expr_lmb.replace_multiple([Replacement(lhs, to_complex(rhs)) for lhs, rhs in eval_function_map.items()])
atom_eval=to_complex(hel_amp_expr_lmb.evaluate_complex(eval_function_map,{}))
print("Atom eval: ", atom_eval)

# %%
eval_function_map = function_map_for_evaluation(test_point_B['kin'], test_point_B['hel'])
tn_eval = tn_replace_multiple(tn_lmb, [(lhs, to_complex(rhs)) for lhs, rhs in eval_function_map.items()] )
tn_eval.execute(hep_lib)
eager_eval = tn_eval.result_scalar()
print("Eager eval: ",eager_eval)

# %%
gL_PS1 = to_complex(hel_amp_expr_lmb.evaluate_complex(function_map_for_evaluation(test_point_A['kin'], test_point_A['hel']),{}))
print("gL_PSA: ", gL_PS1)
print("MG_PSA: ", to_complex(test_point_A['MG_target']))
print(gL_PS1 / to_complex(test_point_A['MG_target']))
gL_PS2 = to_complex(hel_amp_expr_lmb.evaluate_complex(function_map_for_evaluation(test_point_B['kin'], test_point_B['hel']),{}))
print("gL_PSB: ", gL_PS2)
print("MG_PSB: ", to_complex(test_point_B['MG_target']))
print(gL_PS2 / to_complex(test_point_B['MG_target']))

# %% [markdown]
# Now build an input graph with custom numerator

# %%
template_graph = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "qqx_aaa_tree_unprocessed.dot"))

attrs = template_graph[0].get_attributes()
attrs["num"] = f'"{expr_to_string(get_numerator(tree_qqx_aaa_graphs[0]))}"'
attrs["overall_factor"] = '"1"'
attrs["projector"] = f'"{expr_to_string(get_projector(tree_qqx_aaa_graphs[0])*get_color_projector(tree_qqx_aaa_graphs[0]))}"'

template_graph[0].set_edge_defaults(num='"1"')
template_graph[0].set_node_defaults(num='"1"')

#print(os.path.join(ROOT_DIR, "dot_graphs", "qqx_aaa_tree_user_numerator_unprocessed.dot"))
with open(os.path.join(ROOT_DIR, "dot_graphs", "qqx_aaa_tree_user_numerator_unprocessed.dot"),'w') as f:
    f.write(str(template_graph[0]))

# %%
run_gammaloop_commands([
    f'import model {SM_MODEL_PATH}',
    f"import graphs {os.path.join(ROOT_DIR,'dot_graphs','qqx_aaa_tree_user_numerator_unprocessed.dot')}",
    f"save dot {os.path.join(ROOT_DIR,'TMP')}",
    f"!cp {os.path.join(ROOT_DIR,'TMP','processes','amplitudes','qqx_aaa_tree_user_numerator_unprocessed','default','qqx_aaa_tree_1.dot')} {os.path.join(ROOT_DIR,'dot_graphs','qqx_aaa_tree_user_numerator.dot')}",
    f"!rm -rf {os.path.join(ROOT_DIR,'TMP_state')}",
    f"!rm -rf {os.path.join(ROOT_DIR,'TMP')}",
], debug=False)

# %%
# Try and evaluate directly this custom numerator version

testg = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "qqx_aaa_tree_user_numerator.dot"))

tree_graph_0_expr = simplify_color(cook_indices(get_numerator(testg[0])*get_projector(testg[0]))) / get_propagator_denominators(testg[0])
#tree_graph_0_expr = E("1") / get_propagator_denominators(testg[0])


tn = TensorNetwork(tree_graph_0_expr, hep_lib)
tn.execute(hep_lib)
hel_amp_expr = tn.result_scalar()
emr_replacements = get_emr_replacements(testg[0])
hel_amp_expr_lmb = hel_amp_expr.replace_multiple([Replacement(lhs, rhs) for lhs, rhs in emr_replacements])
tn_emr = TensorNetwork(tree_graph_0_expr, hep_lib)
tn_lmb = tn_replace_multiple(tn_emr, emr_replacements)

eval_function_map = function_map_for_evaluation(test_point_B['kin'], test_point_B['hel'])
atom_eval = to_complex(hel_amp_expr_lmb.evaluate_complex(eval_function_map,{}))
print(atom_eval)


# %% [markdown]
# ### Squared matrix element computation

# %% [markdown]
# Let us first compute the full squared matrix element for the process $d \bar{d} \rightarrow \gamma \gamma \gamma$:

# %%
input_kinematics = kin_point_onshell

avg_factor = E("3") # color trace
avg_factor *= E("1/3")**2 # color averaging factor
avg_factor *= E("1/2")**2 # helicity averaging factor
avg_factor *= E("1/(3*2*1)") # final-state symmetry factor

display_sample_point(input_kinematics, None, '"onshell 1"')

total_squared_me = E("0")
for hels in itertools.product(*[[-1,1],]*5):
#for hels in [helicities,]:
    me_for_hel = E("0")
    for external_permutation in itertools.permutations(range(2,5)):
        kinematics_for_this_graph = [
            input_kinematics[0],
            input_kinematics[1],
            input_kinematics[external_permutation[0]],
            input_kinematics[external_permutation[1]],
            input_kinematics[external_permutation[2]]
        ]
        helicities_for_this_graph = [
            hels[0],
            hels[1],
            hels[external_permutation[0]],
            hels[external_permutation[1]],
            hels[external_permutation[2]]
        ]
        eval_function_map_for_this_graph = function_map_for_evaluation(kinematics_for_this_graph, helicities_for_this_graph)
        me_for_this_hel_and_this_diag = to_complex(hel_amp_expr_lmb.evaluate_complex(eval_function_map_for_this_graph,{}))
        #print(f"ME amplitude for helicity: [{', '.join('%+2d'%h for h in helicities_for_this_graph)}] and this permutation [{', '.join('%d'%p for p in external_permutation)}]= ", me_for_this_hel_and_this_diag)
        me_for_hel += me_for_this_hel_and_this_diag

    squared_me_for_this_hel = me_for_hel * me_for_hel.conj() * avg_factor
    #print(f"ME squared for helicity: [{', '.join('%+2d'%h for h in helicities_for_this_graph)}] = ", squared_me_for_this_hel)

    total_squared_me += squared_me_for_this_hel
print('')
print("Total ME squared: %.16e"%total_squared_me.evaluate_complex({},{}).real)

# %% [markdown]
# ## Loop amplitude

# %% [markdown]
# ### Load original loop contributions

# %%
run_gammaloop_commands([
    f'import model {SM_MODEL_PATH}',
    f"import graphs {os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "qqx_aaa_loops_unprocessed.dot")}",
    f"save dot {os.path.join(ROOT_DIR,'TMP')} -c",
    f"!cp {os.path.join(ROOT_DIR,'TMP','processes','amplitudes','qqx_aaa_loops_unprocessed','default','default_graphs.dot')} {os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "qqx_aaa_loops.dot")}",
    f"!rm -rf {os.path.join(ROOT_DIR,'TMP_state')}",
    f"!rm -rf {os.path.join(ROOT_DIR,'TMP')}"
])

# %%
qqx_aaa_loop_graphs = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "qqx_aaa_loops.dot"))
loop_graphs_data = { g.get_name(): {'dot': g} for g in qqx_aaa_loop_graphs }

# %%
len(qqx_aaa_loop_graphs)

# %% [markdown]
# Now build the full integrand representation for each of those graphs

# %%
# g_name, g = list(loop_graphs_data.items())[0]

# loop_tn_expr = simplify_color(
#     cook_indices(
#         get_numerator(g['dot'])
#         *
#         get_color_projector(g['dot'])
#     )
# )

# loop_tn_expr *= cook_indices(get_projector(g['dot']))
# loop_tn_expr /= get_propagator_denominators(g['dot'])

# %%
# tn_emr = TensorNetwork(loop_tn_expr, hep_lib)

# %%
# emr_replacements = get_emr_replacements(g['dot'])

# %%
# tn_lmb = tn_replace_multiple(tn_emr, emr_replacements)

# %%
# tn_lmb.execute(hep_lib)

# %%
# for i_g, (g_name, g) in enumerate(loop_graphs_data.items()):
    
#     #print(f"Processing loop graph {g_name}")
#     #print(f"g:",g['dot'])
#     loop_tn_expr = simplify_color(
#         cook_indices(
#             get_numerator(g['dot'])
#             *
#             get_color_projector(g['dot'])
#         )
#     )
#     #tree_graph_0_expr= simplify_gamma(tree_graph_0_expr)
#     loop_tn_expr *= cook_indices(get_projector(g['dot']))
#     loop_tn_expr /= get_propagator_denominators(g['dot'])

#     g['tn_expr'] = loop_tn_expr

#     #print("Loop TN expr: ", loop_tn_expr)
#     tn_emr = TensorNetwork(loop_tn_expr, hep_lib)
#     g['tn_emr'] = TensorNetwork(loop_tn_expr, hep_lib)

#     emr_replacements = get_emr_replacements(g['dot'])
#     g['emr_replacements'] = emr_replacements
    
#     #print("Loop tn_emr: ", tn_emr)
#     tn_lmb = tn_replace_multiple(tn_emr, emr_replacements)
#     #print("Loop tn_lmb: ", tn_lmb)
#     g['tn_lmb'] = tn_replace_multiple(tn_emr, emr_replacements)

#     #print("Loop tn_lmb: ", tn_lmb)

#     tn_lmb.execute(hep_lib)
#     # print("prescalar: ", tn_lmb)
#     #print("prescalar scalar: ", tn_lmb.result_scalar())
#     #g['expr_lmb'] = tn_lmb.result_scalar()

# %%
for g_name, g in loop_graphs_data.items():
    #print(f"Processing loop graph {g_name}")
    #print(f"g:",g['dot'])
    # print("loop_tn_expr before: ", cook_indices(
    #         get_numerator(g['dot'])
    #         *
    #         get_color_projector(g['dot'])
    #     ).to_canonical_string())
    loop_tn_expr = simplify_color(
        cook_indices(
            get_numerator(g['dot'])
            *
            get_color_projector(g['dot'])
        )
    )
    # print("loop_tn_expr after: ", loop_tn_expr.to_canonical_string())
    #tree_graph_0_expr= simplify_gamma(tree_graph_0_expr)
    loop_tn_expr *= cook_indices(get_projector(g['dot']))
    loop_tn_expr /= get_propagator_denominators(g['dot'])

    # Fix current spenso issue
    #loop_tn_expr= loop_tn_expr.replace(Es("spenso::gamma(a_,b_,c_)"),Es("spenso::gamma(b_,a_,c_)")) 
    
    g['tn_expr'] = loop_tn_expr

    # print("Loop TN expr: ", loop_tn_expr)
    tn_emr = TensorNetwork(loop_tn_expr, hep_lib)
    g['tn_emr'] = TensorNetwork(loop_tn_expr, hep_lib)
    
    emr_replacements = get_emr_replacements(g['dot'])
    g['emr_replacements'] = emr_replacements
    
    # print("Loop tn_emr: ", tn_emr)
    tn_lmb = tn_replace_multiple(tn_emr, emr_replacements)
    # print("Loop tn_lmb: ", tn_lmb)
    g['tn_lmb'] = tn_replace_multiple(tn_emr, emr_replacements)

    tn_lmb.execute(hep_lib)
    # print("prescalar: ", tn_lmb)
    # print("prescalar scalar: ", tn_lmb.result_scalar())
    g['expr_lmb'] = tn_lmb.result_scalar()


# %% [markdown]
# Now run a test evaluation

# %%
for g_name, data in loop_graphs_data.items():
    eval_res = data["expr_lmb"].evaluate_complex(
        function_map_for_evaluation(kin_point_onshell,helicities_onshell,loop_mom=[0.1,0.2,0.3,0.4]),
        {}
    )
    print(f"Evaluation result for graph {g_name:20s}: {eval_res}")

# %% [markdown]
# Define the P1 and P2 projectors

# %%
P1Proj = E("spenso::gamma( spenso::bis(4,proj_bis_indices(left)), spenso::bis(4,dummy_bis_indices(1)), spenso::mink(4,dummy_lor_indices(1)) ) * gammalooprs::P(1,spenso::mink(4,dummy_lor_indices(1)))")
P1Proj *= E("spenso::gamma( spenso::bis(4,dummy_bis_indices(1)), spenso::bis(4,proj_bis_indices(right)), spenso::mink(4,dummy_lor_indices(2)) ) * gammalooprs::P(2,spenso::mink(4,dummy_lor_indices(2)))")
P1Proj /= E("2*(gammalooprs::P(1,spenso::cind(0))*gammalooprs::P(2,spenso::cind(0))-gammalooprs::P(1,spenso::cind(1))*gammalooprs::P(2,spenso::cind(1))-gammalooprs::P(1,spenso::cind(2))*gammalooprs::P(2,spenso::cind(2))-gammalooprs::P(1,spenso::cind(3))*gammalooprs::P(2,spenso::cind(3)))")

P2Proj = E("spenso::g(proj_bis_indices(left), proj_bis_indices(right))") - P1Proj

# %%
# Include the CF factor modified w.r.t Babis QED example
tree_amplitude_proj_emr = (simplify_color(E("spenso::CF")*get_color_projector(tree_qqx_aaa_graphs[0])*get_numerator(tree_qqx_aaa_graphs[0]))/get_propagator_denominators(tree_qqx_aaa_graphs[0]))
tree_amplitude_proj = tree_amplitude_proj_emr.replace_multiple([Replacement(lhs, rhs) for lhs, rhs in get_emr_replacements(tree_qqx_aaa_graphs[0])])

# %% [markdown]
# ### Build the IR triangle counterterm

# %%
ir_ct_template = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "ir_ct_template.dot"))[0]


# %%
tree_amplitude_proj_manipulated= tree_amplitude_proj.replace_multiple([
    Replacement(E("spenso::bis(4,gammalooprs::hedge(1))"),E("spenso::bis(4,tree_form_factor_spinor_1)")),
    Replacement(E("spenso::bis(4,gammalooprs::hedge(2))"),E("spenso::bis(4,tree_form_factor_spinor_2)"))
])

# Now add the actual polarization vectors
ir_ct_proj = tree_amplitude_proj_manipulated * get_projector(tree_qqx_aaa_graphs[0])

# %%
ir_ct_num = E("spenso::gamma( spenso::bis(4,gammalooprs::hedge(2)), spenso::bis(4,s1), spenso::mink(4,mu1) )")
# TODO investigate what this naive piece is not correct
#ir_ct_num *= E("spenso::gamma( spenso::bis(4,s1), spenso::bis(4,s2), spenso::mink(4,mu2) ) * ( gammalooprs::P(2,spenso::mink(4,mu2)) - gammalooprs::Q(7,spenso::mink(4,mu2)) ) ")
# And instead this one is correct
ir_ct_num *= E("spenso::gamma( spenso::bis(4,s1), spenso::bis(4,s2), spenso::mink(4,mu2) ) * ( -gammalooprs::P(2,spenso::mink(4,mu2)) - gammalooprs::Q(7,spenso::mink(4,mu2)) ) ")
ir_ct_num *= P1Proj.replace_multiple([
    Replacement(E("proj_bis_indices(left)"), E("s2")),
    Replacement(E("proj_bis_indices(right)"), E("tree_form_factor_spinor_2")),
    Replacement(E("dummy_bis_indices(1)"), E("s3")),
    Replacement(E("dummy_lor_indices(1)"), E("mu3")),
    Replacement(E("dummy_lor_indices(2)"), E("mu4")),
])

ir_ct_num *= P1Proj.replace_multiple([
    Replacement(E("proj_bis_indices(left)"), E("tree_form_factor_spinor_1")),
    Replacement(E("proj_bis_indices(right)"), E("s4")),
    Replacement(E("dummy_bis_indices(1)"), E("s5")),
    Replacement(E("dummy_lor_indices(1)"), E("mu5")),
    Replacement(E("dummy_lor_indices(2)"), E("mu6")),
])
# TODO investigate what this naive piece is not correct
#ir_ct_num *= E("spenso::gamma( spenso::bis(4,s4), spenso::bis(4,s6), spenso::mink(4,mu7) ) * ( gammalooprs::P(1,spenso::mink(4,mu7)) + ( -gammalooprs::Q(6,spenso::mink(4,mu7)) - gammalooprs::P(2,spenso::mink(4,mu7)) ) ) ")
# And instead this one is correct
ir_ct_num *= E("spenso::gamma( spenso::bis(4,s4), spenso::bis(4,s6), spenso::mink(4,mu7) ) * ( gammalooprs::P(1,spenso::mink(4,mu7)) - ( -gammalooprs::Q(6,spenso::mink(4,mu7)) - gammalooprs::P(2,spenso::mink(4,mu7)) ) ) ")
ir_ct_num *= E("spenso::gamma( spenso::bis(4,s6), spenso::bis(4,gammalooprs::hedge(1)), spenso::mink(4,mu1) )")
ir_ct_num *= E("UFO::G^2")

# %%

attrs = ir_ct_template.get_attributes()
attrs["num"] = f'"{expr_to_string(-ir_ct_num)}"'
attrs["overall_factor"] = '"1"'
attrs["projector"] = f'"{expr_to_string(ir_ct_proj)}"'
ir_ct_template.set_edge_defaults(num='"1"', dod=-100)
ir_ct_template.set_node_defaults(num='"1"', dod=-100)
with open(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "ir_ct_unprocessed.dot"),'w') as f:
    f.write(str(ir_ct_template))

# %% [markdown]
# Process this input dot file with gammaloop, in order to dress it with the emr to lmb replacement rules

# %% [markdown]
# Now evaluate the resulting triangle IR CT

# %%
run_gammaloop_commands([
    f'import model {SM_MODEL_PATH}',
    f"import graphs {os.path.join(ROOT_DIR,'dot_graphs','loop_graphs','ir_ct_unprocessed.dot')}",
    f"save dot {os.path.join(ROOT_DIR,'TMP')} -c",
    f"!cp {os.path.join(ROOT_DIR,'TMP','processes','amplitudes','ir_ct_unprocessed','default','default_graphs.dot')} {os.path.join(ROOT_DIR,'dot_graphs','loop_graphs','ir_ct.dot')}",
    f"!rm -rf {os.path.join(ROOT_DIR,'TMP_state')}",
    f"!rm -rf {os.path.join(ROOT_DIR,'TMP')}",
])

# %%
ir_ct = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "ir_ct.dot"))[0]

# %%
loop_graphs_data['ir_ct'] = {'dot': ir_ct}

loop_tn_expr = simplify_color(
    cook_indices(
        get_numerator(ir_ct)
        *
        get_projector(ir_ct)
    )
)
#tree_graph_0_expr= simplify_gamma(tree_graph_0_expr)
loop_tn_expr /= get_propagator_denominators(ir_ct)

loop_graphs_data['ir_ct']['tn_expr'] = loop_tn_expr

tn_emr = TensorNetwork(loop_tn_expr, hep_lib)
loop_graphs_data['ir_ct']['tn_emr'] = TensorNetwork(loop_tn_expr, hep_lib)

emr_replacements = get_emr_replacements(ir_ct)
loop_graphs_data['ir_ct']['emr_replacements'] = emr_replacements

tn_lmb = tn_replace_multiple(tn_emr, emr_replacements)
loop_graphs_data['ir_ct']['tn_lmb'] = tn_replace_multiple(tn_emr, emr_replacements)

tn_lmb.execute(hep_lib)
loop_graphs_data['ir_ct']['expr_lmb'] = tn_lmb.result_scalar()

# %%
#print('\n'.join([k.to_canonical_string()+' -> '+v.to_canonical_string() for k,v in loop_graphs_data['qqx_aaa_pentagon']['emr_replacements']]))

# %%
for g_name, data in loop_graphs_data.items():
    if g_name != 'ir_ct':
        continue
    eval_res = data["expr_lmb"].evaluate_complex(
        function_map_for_evaluation(kin_point_onshell,helicities_onshell,loop_mom=[0.1,0.2,0.3,0.4]),
        {}
    )
    print(f"Evaluation result for graph {g_name:20s}: {eval_res}")

# %% [markdown]
# ### Build UV counterterms

# %%
uv_counterterms = []

# %%
def gamma_builder(i,j,mu, hedge_indices=False):
    if hedge_indices:
        return E(f"spenso::gamma(spenso::bis(4,gammalooprs::hedge({i})),spenso::bis(4,gammalooprs::hedge({j})),spenso::mink(4,gammalooprs::hedge({mu})))")
    else:
        return E(f"spenso::gamma(spenso::bis(4,{i}),spenso::bis(4,{j}),spenso::mink(4,{mu}))")

def substitute_with_form_factor(expr, lhs, rhs_ff, indices):
    substituted_expr = expr.replace(lhs, rhs_ff)

    if len(indices) > 0:
        substituted_expr *= E(f"spenso::g(spenso::bis(4,{expr_to_string(indices[0])}),spenso::bis(4,uv_form_factor(1)))")
    if len(indices) > 1:
        substituted_expr *= E(f"spenso::g(spenso::bis(4,uv_form_factor(2)),spenso::bis(4,{expr_to_string(indices[1])}))")
    if len(indices) > 2:
        substituted_expr *= E(f"spenso::g(spenso::mink(4,uv_form_factor(3)),spenso::mink(4,{expr_to_string(indices[2])}))")
    # for i, idx in enumerate(indices):
        #substituted_expr = substituted_expr.replace(idx, E(f"uv_form_factor({i+1})"))
    return substituted_expr

# %%
# Factor -1 is for the overall orientation of the loop legs in the bubble
uv_triangle_form_factor = E("-1")*gamma_builder("uv_form_factor(1)","uv_bis(1)","uv_mink(1)") \
    * ( gamma_builder("uv_bis(1)","uv_bis(2)","uv_mink(2)") * E("gammalooprs::Q(5,spenso::mink(4,uv_mink(2)))")  ) \
    * gamma_builder("uv_bis(2)","uv_bis(3)","uv_form_factor(3)") \
    * ( gamma_builder("uv_bis(3)","uv_bis(4)","uv_mink(3)") * E("gammalooprs::Q(6,spenso::mink(4,uv_mink(3)))") ) \
    * gamma_builder("uv_bis(4)","uv_form_factor(2)","uv_mink(1)") \
    * E("UFO::G^2")

# %%
# print('\n*'.join(str(tree_amplitude_proj).split('*')))
# print('\n*'.join(str(get_projector(tree_qqx_aaa_graphs[0])).split('*')))

# %%
tree_amplitude_proj_uv_triangle_C = substitute_with_form_factor(
    tree_amplitude_proj,
    gamma_builder(2,8,0, hedge_indices=True),
    uv_triangle_form_factor,
    [ E("gammalooprs::hedge(2)"), E("gammalooprs::hedge(8)"), E("gammalooprs::hedge(0)") ]
)
tree_amplitude_proj_uv_triangle_B = substitute_with_form_factor(
    tree_amplitude_proj,
    gamma_builder(7,6,4, hedge_indices=True),
    uv_triangle_form_factor,
    [ E("gammalooprs::hedge(7)"), E("gammalooprs::hedge(6)"), E("gammalooprs::hedge(4)") ]
)
tree_amplitude_proj_uv_triangle_A = substitute_with_form_factor(
    tree_amplitude_proj,
    gamma_builder(5,1,3, hedge_indices=True),
    uv_triangle_form_factor,
    [ E("gammalooprs::hedge(5)"), E("gammalooprs::hedge(1)"), E("gammalooprs::hedge(3)") ]
)

# %%
# print('\n*'.join(str(tree_amplitude_proj_uv_triangle_A.to_canonical_string()).split('*')))

# %%
for uv_ct_name, uv_ct in zip(
    ['uv_triangle_A','uv_triangle_B','uv_triangle_C'],
    [tree_amplitude_proj_uv_triangle_A, tree_amplitude_proj_uv_triangle_B, tree_amplitude_proj_uv_triangle_C]
    ):
    uv_ct_template = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "uv_ct_template.dot"))[0]
    uv_ct_template.set_name(uv_ct_name)
    uv_ct_template.set_edge_defaults(num='"1"', dod=-100)
    uv_ct_template.set_node_defaults(num='"1"', dod=-100)
    attrs = uv_ct_template.get_attributes()
    attrs["num"] = f'"{expr_to_string(-uv_ct)}"' 
    attrs["overall_factor"] = '"1"'
    attrs["projector"] = f'"{expr_to_string(get_projector(tree_qqx_aaa_graphs[0]))}"'
    uv_counterterms.append(uv_ct_template)


# %%
def uv_bubble_form_factor(p):
    # Factor -1 is for the overall orientation of the loop legs in the bubble
    return E("-1")*E("UFO::G^2") * gamma_builder("uv_form_factor(1)","uv_bis(1)","uv_mink(1)") * (
        ( ( gamma_builder("uv_bis(1)","uv_bis(4)","uv_mink(2)") * E("gammalooprs::Q(5,spenso::mink(4,uv_mink(2)))")  ) * (
            E("-gammalooprs::Q(6,spenso::mink(4,uv_mink(1337)))") * E("-gammalooprs::Q(7,spenso::mink(4,uv_mink(1337)))") - E("gammalooprs::mUV^2")
        ) )
        - (
            ( gamma_builder("uv_bis(1)","uv_bis(2)","uv_mink(2)") * E("gammalooprs::Q(5,spenso::mink(4,uv_mink(2)))")  ) 
            * gamma_builder("uv_bis(2)","uv_bis(3)","uv_mink(5)") * p.replace(E('mu'),E('spenso::mink(4,uv_mink(5))'))
            * ( gamma_builder("uv_bis(3)","uv_bis(4)","uv_mink(3)") * E("gammalooprs::Q(6,spenso::mink(4,uv_mink(3)))") ) 
        )
    ) * gamma_builder("uv_bis(4)","uv_form_factor(2)","uv_mink(1)") 

# %%
tree_amplitude_proj_uv_bubble_A = substitute_with_form_factor(
    tree_amplitude_proj,
    gamma_builder("gammalooprs::hedge(6)","gammalooprs::hedge(5)","gammalooprs::edge(5,1)"),
    ( gamma_builder("gammalooprs::hedge(6)","uv_form_factor(1)","gammalooprs::edge(5,1)")
      * uv_bubble_form_factor(E("( gammalooprs::P(1, mu) - gammalooprs::P(3, mu) )"))
      * gamma_builder("uv_form_factor(2)","gammalooprs::hedge(5)","uv_mink(10)") 
      * E("( gammalooprs::P(1, spenso::mink(4,uv_mink(10))) - gammalooprs::P(3, spenso::mink(4,uv_mink(10))) )")
      * E("((gammalooprs::P(1,spenso::cind(0))-gammalooprs::P(3,spenso::cind(0)))^2-(gammalooprs::P(1,spenso::cind(1))-gammalooprs::P(3,spenso::cind(1)))^2-(gammalooprs::P(1,spenso::cind(2))-gammalooprs::P(3,spenso::cind(2)))^2-(gammalooprs::P(1,spenso::cind(3))-gammalooprs::P(3,spenso::cind(3)))^2)^-1")
    ),
    [ ]
)

tree_amplitude_proj_uv_bubble_B = substitute_with_form_factor(
    tree_amplitude_proj,
    gamma_builder("gammalooprs::hedge(8)","gammalooprs::hedge(7)","gammalooprs::edge(6,1)"),
    ( gamma_builder("gammalooprs::hedge(8)","uv_form_factor(1)","gammalooprs::edge(6,1)")
      * uv_bubble_form_factor(E("( gammalooprs::P(1, mu) - gammalooprs::P(3, mu) - gammalooprs::P(4, mu) )"))
      * gamma_builder("uv_form_factor(2)","gammalooprs::hedge(7)","uv_mink(10)") 
      * E("( gammalooprs::P(1,spenso::mink(4,uv_mink(10))) - gammalooprs::P(3,spenso::mink(4,uv_mink(10))) - gammalooprs::P(4,spenso::mink(4,uv_mink(10))) )")
      * E("((gammalooprs::P(1,spenso::cind(0))-gammalooprs::P(3,spenso::cind(0))-gammalooprs::P(4,spenso::cind(0)))^2-(gammalooprs::P(1,spenso::cind(1))-gammalooprs::P(3,spenso::cind(1))-gammalooprs::P(4,spenso::cind(1)))^2-(gammalooprs::P(1,spenso::cind(2))-gammalooprs::P(3,spenso::cind(2))-gammalooprs::P(4,spenso::cind(2)))^2-(gammalooprs::P(1,spenso::cind(3))-gammalooprs::P(3,spenso::cind(3))-gammalooprs::P(4,spenso::cind(3)))^2)^-1")
    ),
    [ ]
)

# %%
for uv_ct_name, uv_ct in zip(
    ['uv_bubble_A','uv_bubble_B'],
    [tree_amplitude_proj_uv_bubble_A, tree_amplitude_proj_uv_bubble_B]
    ):
    uv_ct_template = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "uv_ct_template.dot"))[0]
    uv_ct_template.set_name(uv_ct_name)
    uv_ct_template.set_edge_defaults(num='"1"', dod=-100)
    uv_ct_template.set_node_defaults(num='"1"', dod=-100)
    attrs = uv_ct_template.get_attributes()
    #print(expr_to_string(-uv_ct))
    attrs["num"] = f'"{expr_to_string(-uv_ct)}"' 
    attrs["overall_factor"] = '"1"'
    attrs["projector"] = f'"{expr_to_string(get_projector(tree_qqx_aaa_graphs[0]))}"'
    uv_counterterms.append(uv_ct_template)

# %%
ir_uv_proj = E("spenso::gamma( spenso::bis(4,gammalooprs::hedge(2)), spenso::bis(4,s1), spenso::mink(4,mu1) )")
ir_uv_proj *= E("spenso::gamma( spenso::bis(4,s1), spenso::bis(4,s2), spenso::mink(4,mu2) ) * ( gammalooprs::Q(5,spenso::mink(4,mu2)) ) ")
ir_uv_proj *= P1Proj.replace_multiple([
    Replacement(E("proj_bis_indices(left)"), E("s2")),
    Replacement(E("proj_bis_indices(right)"), E("tree_form_factor_spinor_2")),
    Replacement(E("dummy_bis_indices(1)"), E("s3")),
    Replacement(E("dummy_lor_indices(1)"), E("mu3")),
    Replacement(E("dummy_lor_indices(2)"), E("mu4")),
])

ir_uv_proj *= P1Proj.replace_multiple([
    Replacement(E("proj_bis_indices(left)"), E("tree_form_factor_spinor_1")),
    Replacement(E("proj_bis_indices(right)"), E("s4")),
    Replacement(E("dummy_bis_indices(1)"), E("s5")),
    Replacement(E("dummy_lor_indices(1)"), E("mu5")),
    Replacement(E("dummy_lor_indices(2)"), E("mu6")),
])
ir_uv_proj *= E("spenso::gamma( spenso::bis(4,s4), spenso::bis(4,s6), spenso::mink(4,mu7) ) * ( -gammalooprs::Q(6,spenso::mink(4,mu7)) )")
ir_uv_proj *= E("spenso::gamma( spenso::bis(4,s6), spenso::bis(4,gammalooprs::hedge(1)), spenso::mink(4,mu1) )")
ir_uv_proj *= E("UFO::G^2")

# %%
uv_ct_template = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "uv_ct_template.dot"))[0]
uv_ct_template.set_name("uv_ir_ct")
tree_amplitude_proj_manipulated= tree_amplitude_proj.replace_multiple([
    Replacement(E("spenso::bis(4,gammalooprs::hedge(1))"),E("spenso::bis(4,tree_form_factor_spinor_1)")),
    Replacement(E("spenso::bis(4,gammalooprs::hedge(2))"),E("spenso::bis(4,tree_form_factor_spinor_2)"))
])
attrs = uv_ct_template.get_attributes()
attrs["num"] = f'"{expr_to_string(ir_uv_proj)}"'
attrs["overall_factor"] = '"1"'
attrs["projector"] = f'"{expr_to_string(tree_amplitude_proj_manipulated * get_projector(tree_qqx_aaa_graphs[0]))}"'
uv_ct_template.set_edge_defaults(num='"1"', dod=-100)
uv_ct_template.set_node_defaults(num='"1"', dod=-100)
uv_counterterms.append(uv_ct_template)

# %% [markdown]
# ### Integrated CTs

# %%
integrated_uv_counterterms = []

# %%
#print('*\n'.join(s for s in str(ct_born_me_dressed).split('*')))
#print('*\n'.join(s for s in str(get_projector(tree_qqx_aaa_graphs[0])).split('*')))

# %%
# For the integrated CT directly proportional to the born
ct_born_me_bare = tree_amplitude_proj*(1/E("spenso::CF"))

ct_born_me_dressed = ct_born_me_bare.replace_multiple([
    Replacement(E("spenso::bis(4,gammalooprs::hedge(1))"),E("spenso::bis(4,ct_connection_spinor_1)")),
    Replacement(E("spenso::bis(4,gammalooprs::hedge(2))"),E("spenso::bis(4,ct_connection_spinor_2)"))
])
# ct_born_me_dressed *= E("spenso::g(spenso::bis(4,gammalooprs::hedge(1)),spenso::bis(4,ct_connection_spinor_1))")
# ct_born_me_dressed *= E("spenso::g(spenso::bis(4,ct_connection_spinor_2),spenso::bis(4,gammalooprs::hedge(2)))")

# For the integrated CT proportional to Gamma[mu] * Gamma[nu] * P1 * Born * P1 * Gamma[nu] * Gamma[mu]
ct_born_me_dressed *= E("spenso::gamma( spenso::bis(4,gammalooprs::hedge(2)), spenso::bis(4,s1), spenso::mink(4,mu1) )")
ct_born_me_dressed *= E("spenso::gamma( spenso::bis(4,s1), spenso::bis(4,s2), spenso::mink(4,mu2) ) ")
ct_born_me_dressed *= P1Proj.replace_multiple([
    Replacement(E("proj_bis_indices(left)"), E("s2")),
    Replacement(E("proj_bis_indices(right)"), E("ct_connection_spinor_2")),
    Replacement(E("dummy_bis_indices(1)"), E("s3")),
    Replacement(E("dummy_lor_indices(1)"), E("mu3")),
    Replacement(E("dummy_lor_indices(2)"), E("mu4")),
])

ct_born_me_dressed *= P1Proj.replace_multiple([
    Replacement(E("proj_bis_indices(left)"), E("ct_connection_spinor_1")),
    Replacement(E("proj_bis_indices(right)"), E("s4")),
    Replacement(E("dummy_bis_indices(1)"), E("s5")),
    Replacement(E("dummy_lor_indices(1)"), E("mu5")),
    Replacement(E("dummy_lor_indices(2)"), E("mu6")),
])
ct_born_me_dressed *= E("spenso::gamma( spenso::bis(4,s4), spenso::bis(4,s6), spenso::mink(4,mu2) ) ")
ct_born_me_dressed *= E("spenso::gamma( spenso::bis(4,s6), spenso::bis(4,gammalooprs::hedge(1)), spenso::mink(4,mu1) )")

triangle_normalization_factor = E("-2 *(4*𝜋)^2 * mUVsq")

# %% [markdown]
# Test of triangle normalization factor

# %%
uv_ct_template = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "uv_ct_template.dot"))[0]
uv_ct_template.set_name("tri_norm_test")
uv_ct_template.set_edge_defaults(num='"1"', dod=-100)
uv_ct_template.set_node_defaults(num='"1"', dod=-100)
attrs = uv_ct_template.get_attributes()
attrs["num"] = f'"{expr_to_string(prepare_integrated_ct_numerator_expression(triangle_normalization_factor))}"' 
attrs["overall_factor"] = '"1"'
attrs["projector"] = '"1"'
integrated_uv_counterterms.append(uv_ct_template)

# %% [markdown]
# Integrated UV CT

# %%
Z1 = E("UFO::G^2 / (4*𝜋)^(dim/2) * 𝚪(1+ε) * (1 - ε)^2 / ε * (mUVsq)^(-ε)")
Z1_expanded = eps_expansion_finite(Z1)
int_uv_ct_num = ct_born_me_bare * Z1_expanded * triangle_normalization_factor

# %%
uv_ct_template = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "uv_ct_template.dot"))[0]
uv_ct_template.set_name("int_uv_ct")
uv_ct_template.set_edge_defaults(num='"1"', dod=-100)
uv_ct_template.set_node_defaults(num='"1"', dod=-100)
attrs = uv_ct_template.get_attributes()
attrs["num"] = f'"{expr_to_string(prepare_integrated_ct_numerator_expression(int_uv_ct_num))}"' 
attrs["overall_factor"] = '"1"'
attrs["projector"] = f'"{expr_to_string(get_projector(tree_qqx_aaa_graphs[0]))}"'
integrated_uv_counterterms.append(uv_ct_template)

# %% [markdown]
# Integrated IR CT

# %%
mandelstamm_s = E("spenso::P(1,spenso::cind(0))*spenso::P(2,spenso::cind(0)) - spenso::P(1,spenso::cind(1))*spenso::P(2,spenso::cind(1)) - spenso::P(1,spenso::cind(2))*spenso::P(2,spenso::cind(2)) - spenso::P(1,spenso::cind(3))*spenso::P(2,spenso::cind(3))")

ir_int_ct_prefactor = E("(UFO::G^2) * -1 / ( 4*𝜋 )^(dim/2) * 𝚪(ε) * 𝚪(1-ε)^2 / 𝚪(2-2*ε) * (-s)^(-ε)")

ir_bare_int_ct_prefactor = E("2 / ε + ( ε / ( 1 - ε ) ) ") * ir_int_ct_prefactor
ir_dressed_int_ct_prefactor = E("-1 / (4*( 1 - ε ))") * ir_int_ct_prefactor

ir_bare_int_ct_prefactor_expanded = eps_expansion_finite(ir_bare_int_ct_prefactor).replace(E("s"),mandelstamm_s)
ir_dressed_int_ct_prefactor_expanded = eps_expansion_finite(ir_dressed_int_ct_prefactor).replace(E("s"),mandelstamm_s)


uv_ct_template = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "uv_ct_template.dot"))[0]
uv_ct_template.set_name("int_ir_bare_ct")
uv_ct_template.set_edge_defaults(num='"1"', dod=-100)
uv_ct_template.set_node_defaults(num='"1"', dod=-100)
attrs = uv_ct_template.get_attributes()
attrs["num"] = f'"{expr_to_string(prepare_integrated_ct_numerator_expression(ct_born_me_bare * ir_bare_int_ct_prefactor_expanded * triangle_normalization_factor))}"' 
attrs["overall_factor"] = '"1"'
attrs["projector"] = f'"{expr_to_string(get_projector(tree_qqx_aaa_graphs[0]))}"'
integrated_uv_counterterms.append(uv_ct_template)

uv_ct_template = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "uv_ct_template.dot"))[0]
uv_ct_template.set_name("int_ir_dressed_ct")
uv_ct_template.set_edge_defaults(num='"1"', dod=-100)
uv_ct_template.set_node_defaults(num='"1"', dod=-100)
attrs = uv_ct_template.get_attributes()
attrs["num"] = f'"{expr_to_string(prepare_integrated_ct_numerator_expression(ct_born_me_dressed * ir_dressed_int_ct_prefactor_expanded * triangle_normalization_factor))}"' 
attrs["overall_factor"] = '"1"'
attrs["projector"] = f'"{expr_to_string(get_projector(tree_qqx_aaa_graphs[0]))}"'
integrated_uv_counterterms.append(uv_ct_template)

# %% [markdown]
# Integrated UV IR CT

# %%
uv_ir_ct_prefactor = Z1 / E("4 * ( 1 - ε )^2")
uv_ir_ct_prefactor_expanded = eps_expansion_finite(uv_ir_ct_prefactor)

uv_ct_template = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "uv_ct_template.dot"))[0]
uv_ct_template.set_name("int_uv_ir_ct")
uv_ct_template.set_edge_defaults(num='"1"', dod=-100)
uv_ct_template.set_node_defaults(num='"1"', dod=-100)
attrs = uv_ct_template.get_attributes()
attrs["num"] = f'"{expr_to_string(prepare_integrated_ct_numerator_expression(ct_born_me_dressed * uv_ir_ct_prefactor_expanded * triangle_normalization_factor))}"' 
attrs["overall_factor"] = '"1"'
attrs["projector"] = f'"{expr_to_string(get_projector(tree_qqx_aaa_graphs[0]))}"'
integrated_uv_counterterms.append(uv_ct_template)

# %% [markdown]
# ### Write dot files to disk

# %% [markdown]
# Now process those dot files with gammaloop

# %%
with open(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "uv_cts_unprocessed.dot"),'w') as f:
    f.write('\n\n\n'.join(str(uv_ct_graph) for uv_ct_graph in uv_counterterms + integrated_uv_counterterms))

# %%
run_gammaloop_commands([
    f'import model {SM_MODEL_PATH}',
    f"import graphs {os.path.join(ROOT_DIR,'dot_graphs','loop_graphs','uv_cts_unprocessed.dot')}",
    f"save dot {os.path.join(ROOT_DIR,'TMP')} -c",
    f"!cp {os.path.join(ROOT_DIR,'TMP','processes','amplitudes','uv_cts_unprocessed','default','default_graphs.dot')} {os.path.join(ROOT_DIR,'dot_graphs','loop_graphs','uv_cts.dot')}",
    f"!rm -rf {os.path.join(ROOT_DIR,'TMP_state')}",
    f"!rm -rf {os.path.join(ROOT_DIR,'TMP')}",
])

# %%
qqx_aaa_uv_ct_graphs = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "uv_cts.dot"))
uv_ct_graphs_data = { g.get_name() : {'dot': g} for g in qqx_aaa_uv_ct_graphs }

for (g_name, g) in uv_ct_graphs_data.items():
    # if g_name != 'int_ir_dressed_ct':
    #     continue
    # import time
    # print("Processing: ", g_name)
    # time.sleep(0.5)
    # if g_name != "uv_ir_ct":
    #     continue
    loop_tn_expr = simplify_color(
        cook_indices(
            get_numerator(g['dot'])
            *
            get_projector(g['dot'])
        )
    )
    loop_tn_expr /= get_propagator_denominators(g['dot'])
    # if g_name == "uv_ir_ct":
    #     #loop_tn_expr = loop_tn_expr.replace(E('gammalooprs::Q(5,spenso::mink(4,mu2))'),E('-gammalooprs::Q(5,spenso::mink(4,mu2))+gammalooprs::P(2,spenso::mink(4,mu2))'))
    #     loop_tn_expr = loop_tn_expr.replace(E('gammalooprs::Q(5,spenso::mink(4,mu2))'),E('gammalooprs::Q(5,spenso::mink(4,mu2))+gammalooprs::P(2,spenso::mink(4,mu2))'))

    g['tn_expr'] = loop_tn_expr

    # print("A: ")
    # time.sleep(0.5)
    # print('loop_tn_expr=\n',loop_tn_expr)

    tn_emr = TensorNetwork(loop_tn_expr, hep_lib)
    g['tn_emr'] = TensorNetwork(loop_tn_expr, hep_lib)

    emr_replacements = get_emr_replacements(g['dot'])
    # print('emr_replacements=\n','\n'.join([f"{lhs} -> {rhs}" for lhs, rhs in emr_replacements]))

    g['emr_replacements'] = emr_replacements

    # print("B: ")
    # print("emr_replacements=",emr_replacements)
    # time.sleep(0.5)
    tn_lmb = tn_replace_multiple(tn_emr, emr_replacements)
    # print(get_propagator_denominators(g['dot']).replace_multiple([Replacement(lhs, rhs) for lhs, rhs in emr_replacements]))

    # print("B1: ")
    # time.sleep(0.5)
    g['tn_lmb'] = tn_replace_multiple(tn_emr, emr_replacements)
    # print("B2: ")
    # time.sleep(0.5)
    
    # print("C: ", tn_lmb)
    # time.sleep(0.5)
    tn_lmb.execute(hep_lib)
    #print(g_name, tn_lmb)
    #stepped_execution(tn_lmb, t_delta=0.1)
    g['expr_lmb'] = tn_lmb.result_scalar()

    # print(g['expr_lmb'])


# %%
# tn=TensorNetwork(E("""
# (
#      (-1*gammalooprs::{}::P(3,spenso::{}::cind(0))+-1*gammalooprs::{}::P(4,spenso::{}::cind(0))+gammalooprs::{}::P(1,spenso::{}::cind(0)))^2
#     +(-1*gammalooprs::{}::P(3,spenso::{}::cind(1))+-1*gammalooprs::{}::P(4,spenso::{}::cind(1))+gammalooprs::{}::P(1,spenso::{}::cind(1)))^2*-1
#     +(-1*gammalooprs::{}::P(3,spenso::{}::cind(2))+-1*gammalooprs::{}::P(4,spenso::{}::cind(2))+gammalooprs::{}::P(1,spenso::{}::cind(2)))^2*-1
#     +(-1*gammalooprs::{}::P(3,spenso::{}::cind(3))+-1*gammalooprs::{}::P(4,spenso::{}::cind(3))+gammalooprs::{}::P(1,spenso::{}::cind(3)))^2*-1
# )^(-1)*(
#      (-1*gammalooprs::{}::P(3,spenso::{}::cind(0))+gammalooprs::{}::P(1,spenso::{}::cind(0)))^2
#     +(-1*gammalooprs::{}::P(3,spenso::{}::cind(1))+gammalooprs::{}::P(1,spenso::{}::cind(1)))^2*-1
#     +(-1*gammalooprs::{}::P(3,spenso::{}::cind(2))+gammalooprs::{}::P(1,spenso::{}::cind(2)))^2*-1
#     +(-1*gammalooprs::{}::P(3,spenso::{}::cind(3))+gammalooprs::{}::P(1,spenso::{}::cind(3)))^2*-1)^(-2)*(
#         (-1*gammalooprs::{}::P(3,spenso::{}::mink(4,gammalooprs::{}::uv_mink_5))+gammalooprs::{}::P(1,spenso::{}::mink(4,gammalooprs::{}::uv_mink_5)))
#             *-1*gammalooprs::{}::Q(5,spenso::{}::mink(4,gammalooprs::{}::uv_mink_2))*gammalooprs::{}::Q(6,spenso::{}::mink(4,gammalooprs::{}::uv_mink_3))
#             *spenso::{}::gamma(spenso::{}::bis(4,gammalooprs::{}::uv_bis_1),spenso::{}::bis(4,gammalooprs::{}::uv_bis_2),spenso::{}::mink(4,gammalooprs::{}::uv_mink_2))
#             *spenso::{}::gamma(spenso::{}::bis(4,gammalooprs::{}::uv_bis_2),spenso::{}::bis(4,gammalooprs::{}::uv_bis_3),spenso::{}::mink(4,gammalooprs::{}::uv_mink_5))
#             *spenso::{}::gamma(spenso::{}::bis(4,gammalooprs::{}::uv_bis_3),spenso::{}::bis(4,gammalooprs::{}::uv_bis_4),spenso::{}::mink(4,gammalooprs::{}::uv_mink_3))
#         +(-1*gammalooprs::{}::mUV^2+gammalooprs::{}::Q(6,spenso::{}::mink(4,gammalooprs::{}::uv_mink_1337))
#             *gammalooprs::{}::Q(7,spenso::{}::mink(4,gammalooprs::{}::uv_mink_1337)))
#             *gammalooprs::{}::Q(5,spenso::{}::mink(4,gammalooprs::{}::uv_mink_2))
#             *spenso::{}::gamma(spenso::{}::bis(4,gammalooprs::{}::uv_bis_1),spenso::{}::bis(4,gammalooprs::{}::uv_bis_4),spenso::{}::mink(4,gammalooprs::{}::uv_mink_2))
#     )
#     *(
#         (gammalooprs::{}::Q(5,spenso::{}::cind(0)))^2+(gammalooprs::{}::Q(5,spenso::{}::cind(1)))^2*-1+(gammalooprs::{}::Q(5,spenso::{}::cind(2)))^2*-1       +(gammalooprs::{}::Q(5,spenso::{}::cind(3)))^2*-1+-1*gammalooprs::{}::mUV^2
#     )^(-1)*(
#         (gammalooprs::{}::Q(6,spenso::{}::cind(0)))^2
#         +(gammalooprs::{}::Q(6,spenso::{}::cind(1)))^2*-1
#         +(gammalooprs::{}::Q(6,spenso::{}::cind(2)))^2*-1
#         +(gammalooprs::{}::Q(6,spenso::{}::cind(3)))^2*-1
#         +-1*gammalooprs::{}::mUV^2
#     )^(-1)
#     *(
#         (gammalooprs::{}::Q(7,spenso::{}::cind(0)))^2
#         +(gammalooprs::{}::Q(7,spenso::{}::cind(1)))^2*-1
#         +(gammalooprs::{}::Q(7,spenso::{}::cind(2)))^2*-1
#         +(gammalooprs::{}::Q(7,spenso::{}::cind(3)))^2*-1
#         +-1*gammalooprs::{}::mUV^2
#     )^(-1)
#     *(
#         -1*gammalooprs::{}::P(3,spenso::{}::mink(4,gammalooprs::{}::edge_5_1))
#         +gammalooprs::{}::P(1,spenso::{}::mink(4,gammalooprs::{}::edge_5_1))
#     )*(
#         -1*gammalooprs::{}::P(3,spenso::{}::mink(4,gammalooprs::{}::edge_6_1))
#         +-1*gammalooprs::{}::P(4,spenso::{}::mink(4,gammalooprs::{}::edge_6_1))
#         +gammalooprs::{}::P(1,spenso::{}::mink(4,gammalooprs::{}::edge_6_1))
#     )*(
#         -1*gammalooprs::{}::P(3,spenso::{}::mink(4,gammalooprs::{}::uv_mink_10))
#         +gammalooprs::{}::P(1,spenso::{}::mink(4,gammalooprs::{}::uv_mink_10))
#     )
#     *1𝑖*UFO::{}::GC_1^3*UFO::{}::G^2
#     *gammalooprs::{}::u(1,spenso::{}::bis(4,gammalooprs::{}::hedge_1))
#     *gammalooprs::{}::vbar(2,spenso::{}::bis(4,gammalooprs::{}::hedge_2))
#     *gammalooprs::{}::ϵbar(0,spenso::{}::mink(4,gammalooprs::{}::hedge_0))
#     *gammalooprs::{}::ϵbar(3,spenso::{}::mink(4,gammalooprs::{}::hedge_3))
#     *gammalooprs::{}::ϵbar(4,spenso::{}::mink(4,gammalooprs::{}::hedge_4))
#     *spenso::{}::CF
#     *spenso::{}::gamma(spenso::{}::bis(4,gammalooprs::{}::hedge_2),spenso::{}::bis(4,gammalooprs::{}::hedge_8),spenso::{}::mink(4,gammalooprs::{}::hedge_0))
#     *spenso::{}::gamma(spenso::{}::bis(4,gammalooprs::{}::hedge_5),spenso::{}::bis(4,gammalooprs::{}::hedge_1),spenso::{}::mink(4,gammalooprs::{}::hedge_3))
#     *spenso::{}::gamma(spenso::{}::bis(4,gammalooprs::{}::hedge_6),spenso::{}::bis(4,gammalooprs::{}::uv_form_factor_1),spenso::{}::mink(4,gammalooprs::{}::edge_5_1))
#     *spenso::{}::gamma(spenso::{}::bis(4,gammalooprs::{}::hedge_7),spenso::{}::bis(4,gammalooprs::{}::hedge_6),spenso::{}::mink(4,gammalooprs::{}::hedge_4))
#     *spenso::{}::gamma(spenso::{}::bis(4,gammalooprs::{}::hedge_8),spenso::{}::bis(4,gammalooprs::{}::hedge_7),spenso::{}::mink(4,gammalooprs::{}::edge_6_1))
#     *spenso::{}::gamma(spenso::{}::bis(4,gammalooprs::{}::uv_bis_4),spenso::{}::bis(4,gammalooprs::{}::uv_form_factor_2),spenso::{}::mink(4,gammalooprs::{}::uv_mink_1))
#     *spenso::{}::gamma(spenso::{}::bis(4,gammalooprs::{}::uv_form_factor_1),spenso::{}::bis(4,gammalooprs::{}::uv_bis_1),spenso::{}::mink(4,gammalooprs::{}::uv_mink_1))
#     *spenso::{}::gamma(spenso::{}::bis(4,gammalooprs::{}::uv_form_factor_2),spenso::{}::bis(4,gammalooprs::{}::hedge_5),spenso::{}::mink(4,gammalooprs::{}::uv_mink_10))
# """, default_namespace="spenso"), hep_lib)
# tn.execute(hep_lib)
# ts = tn.result_scalar()
# print(ts.to_canonical_string())
# print('gammalooprs::{}::Q(7,spenso::{}::mink(4,gammalooprs::{}::uv_mink_1337))' in str(ts.to_canonical_string()))



# %%
all_graphs_data = { g_name: g_data for g_name, g_data in list(loop_graphs_data.items())+list(uv_ct_graphs_data.items()) }

# Pick a nice order
all_ordered_graphs_dat = {
    g_name : all_graphs_data[g_name] for g_name in
    [
        "qqx_aaa_pentagon",
        "qqx_aaa_box_A"   ,
        "qqx_aaa_box_B"   ,
        "qqx_aaa_tri_A"   ,
        "uv_triangle_A"   ,
        "qqx_aaa_tri_B"   ,
        "uv_triangle_B"   ,
        "qqx_aaa_tri_C"   ,
        "uv_triangle_C"   ,
        "qqx_aaa_bub_A"   ,
        "uv_bubble_A"     ,
        "qqx_aaa_bub_B"   ,
        "uv_bubble_B"     ,
        "ir_ct"           ,
        "uv_ir_ct"        ,
        "int_ir_bare_ct"  ,
        "int_ir_dressed_ct",
        "int_uv_ir_ct"    ,
        "int_uv_ct"       ,
        "tri_norm_test"
    ]
}

# %% [markdown]
# ### Now test collinear limits

# %%
scaling_factor = 0.1
mUV = 1.0
#test_point = (kin_point, helicities)
#test_point = (kin_point_onshell, helicities_onshell)
test_point = (kin_point2_onshell, helicities2_onshell)

#test_point = (kin_point2_onshell, [-1, 1, -1, -1, 1])

prefactor_power = 2.0
x = 0.1
for test_name, k_base, scaling, pref_power in [
    ('IR limit', [1.0*abs(x),-0.0001,0.0001,-1.0*x], 0.1, 2.0),
    #('UV limit', [100.0,100.0,100.0,100.0], 10.0, 4.0),
]:
    running_sum_non_ct = complex(0.0,0.0)
    running_sum_ct = complex(0.0,0.0)
    running_sum_non_ct_scaled = complex(0.0,0.0)
    running_sum_ct_scaled = complex(0.0,0.0)

    print('')
    print(f'Investigating {test_name}...')
    k = k_base
    scaling_factor = scaling
    prefactor_power = pref_power
    for g_name, data in all_ordered_graphs_dat.items():

        expr = data['expr_lmb']

        # if g_name in ['ir_ct', 'uv_ir_ct']:
        #     #print(f'{g_name}:\n{"\n*".join(str(data['tn_expr']).split('*'))}')
        #     print(f'{g_name}:\n{"\n*".join(f"{lhs} -> {rhs}" for lhs, rhs in data['emr_replacements'])}')
        # if g_name in ['qqx_aaa_tri_A', 'uv_triangle_A']:
        #         print(f'{g_name}:\n{"\n*".join(str(data['tn_expr']).split('*'))}')
        #         print(f'{g_name}:\n{"\n*".join(f"{lhs} -> {rhs}" for lhs, rhs in data['emr_replacements'])}')
        # if g_name in ['uv_bubble_A']:
        #     print("Doing", g_name)
        #     print("with:\n",data['tn_expr'].to_canonical_string())
        #     print("with:\n",expr)
        eval_res = expr.evaluate_complex(
            function_map_for_evaluation(test_point[0],test_point[1],loop_mom=k, mUV = mUV),
            {}
        )
        if test_name == 'IR limit':
            scaled_k = [k[0],k[1]*scaling_factor,k[2]*scaling_factor,k[3]]
        elif test_name == 'UV limit':
            scaled_k = [k[0]*scaling_factor,k[1]*scaling_factor,k[2]*scaling_factor,k[3]*scaling_factor]
        eval_res_scaled = expr.evaluate_complex(
            function_map_for_evaluation(test_point[0],test_point[1],loop_mom=scaled_k, mUV = mUV),
            {}
        )
        # if g_name == 'ir_ct':
        #     eval_res *= -(1.0+x)/(1.0-x)
        #     eval_res_scaled *= -(1.0+x)/(1.0-x)

        eval_res_scaled *= scaling_factor**pref_power
        if g_name in ['ir_ct', 'uv_triangle_A', 'uv_triangle_B', 'uv_triangle_C', 'uv_bubble_A', 'uv_bubble_B', 'uv_ir_ct']:
            running_sum_ct += eval_res
            running_sum_ct_scaled += eval_res_scaled
        else:
            running_sum_non_ct += eval_res
            running_sum_non_ct_scaled += eval_res_scaled
        print(f"{'Evaluation result for graph "'+g_name+'"':50s}: {eval_res:50.16e} | scaled: {eval_res_scaled:50.16e} | log10(abs(scaled / orig)) = {log10(abs(eval_res_scaled / eval_res)) if eval_res!=0. else float('nan'):.3f} ")
    print('-'*100)
    print(f"{'Total evaluation result from non-CT':50s}: {running_sum_non_ct:50.16e} | scaled: {running_sum_non_ct_scaled:50.16e}  | log10(abs(scaled / orig)) = {log10(abs(running_sum_non_ct_scaled / running_sum_non_ct)) if running_sum_non_ct!=0. else float('nan'):.3f} ")
    print(f"{'Total evaluation result from CT':50s}: {running_sum_ct:50.16e} | scaled: {running_sum_ct_scaled:50.16e}  | log10(abs(scaled / orig)) = {log10(abs(running_sum_ct_scaled / running_sum_ct)) if running_sum_non_ct!=0. else float('nan'):.3f} ")
    print(f"{'Total evaluation result':50s}: {running_sum_ct+running_sum_non_ct:50.16e} | scaled: {running_sum_ct_scaled+running_sum_non_ct_scaled:50.16e}  | log10(abs(scaled / orig)) = {log10(abs((running_sum_ct_scaled+running_sum_non_ct_scaled) / (running_sum_ct+running_sum_non_ct))) if running_sum_ct+running_sum_non_ct!=0. else float('nan'):.3f} ")
    print(f"{'Ratio non-CT / CT ':50s}: {running_sum_non_ct / running_sum_ct:50.16e} | scaled: {running_sum_non_ct_scaled / running_sum_ct_scaled:50.16e}  | log10(abs(scaled / orig)) = {log10(abs((running_sum_non_ct_scaled / running_sum_ct_scaled) / (running_sum_non_ct / running_sum_ct))) if running_sum_ct_scaled!=0. and running_sum_ct!=0. and running_sum_non_ct!=0. else float('nan'):.3f} ")

    print('')

# %% [markdown]
# ## Testing final diagrams on file

# %% [markdown]
# Combine all graphs and write them into a single dot file named 'qqx_aaa_subtracted.dot'

# %%
def get_fudge_factor(g_name):
    if g_name in     [
        "uv_triangle_A"   ,
        "uv_triangle_B"   ,
        "uv_triangle_C"   ,
        "uv_bubble_A"     ,
        "uv_bubble_B"     ,
        "ir_ct"           ,
        "uv_ir_ct"        ,
        "int_ir_bare_ct"  ,
        "int_ir_dressed_ct",
        "int_uv_ir_ct"    ,
        "int_uv_ct"
    ]:
        return complex(0.0,-1.0)
    elif g_name in ['tri_norm_test']:
        return complex(1.0,0.0)        
    else:
        return complex(1.0,0.0)

def preprocess_expression_for_gammaloop(expr):
    expr = expr.replace(E("gammalooprs::P(i_,a___)"),E("gammalooprs::Q(i_,a___)"))
    expr = expr.replace(E("vakint::EulerGamma"),to_complex(complex(0.577215664901533,0.0)))
    return cook_indices(simplify_metrics(expr)).replace(E('spenso::CF'),E('4/3'))

# %%
loop_graphs = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "qqx_aaa_loops.dot"))

# Add the color projectors to the loop_graphs generated from gammaloop
for lg in loop_graphs:
    lg_attrs = lg.get_attributes()
    lg_attrs["projector"] = f'{lg_attrs["projector"].replace('"','')}*({get_color_projector(lg).to_canonical_string()})'
    if lg.get_name() in ['qqx_aaa_pentagon',]:
        lg_attrs["group_master"] = '"true"'
    for e in lg.get_edges():
        e_attrs = e.get_attributes()
        if "dod" in e_attrs:
            e_attrs["dod"] = -100
    for n in lg.get_nodes():
        n_attrs = n.get_attributes()
        if "dod" in n_attrs:
            n_attrs["dod"] = -100

ir_ct = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "ir_ct.dot"))[0]
uv_cts = pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "uv_cts.dot"))

all_graphs= {g.get_name(): g for g in loop_graphs + [ir_ct,] + uv_cts}
all_ordered_graphs = {
    g_name : all_graphs[g_name] for g_name in
    [
        "qqx_aaa_pentagon",
        "qqx_aaa_box_A"   ,
        "qqx_aaa_box_B"   ,
        "qqx_aaa_tri_A"   ,
        "uv_triangle_A"   ,
        "qqx_aaa_tri_B"   ,
        "uv_triangle_B"   ,
        "qqx_aaa_tri_C"   ,
        "uv_triangle_C"   ,
        "qqx_aaa_bub_A"   ,
        "uv_bubble_A"     ,
        "qqx_aaa_bub_B"   ,
        "uv_bubble_B"     ,
        "ir_ct"           ,
        "uv_ir_ct"        ,
        "int_ir_bare_ct"  ,
        "int_ir_dressed_ct",
        "int_uv_ir_ct"    ,
        "int_uv_ct"       ,
        "tri_norm_test"
    ]
}
        

for g_name, g in all_ordered_graphs.items():
    g_attrs = g.get_attributes()

    fudge_factor = to_complex(get_fudge_factor(g_name))

    g_attrs["overall_factor"] = f'{expr_to_string(preprocess_expression_for_gammaloop(fudge_factor))}'
    g_attrs["num"] = f'{(expr_to_string(preprocess_expression_for_gammaloop(Es(g_attrs["num"]))))}'
    g_attrs["projector"] = f'{expr_to_string(preprocess_expression_for_gammaloop(Es(g_attrs["projector"])))}'
    if g_name != "tri_norm_test":
        g_attrs["group_id"] = 0

    for e in g.get_edges():
        e_attrs = e.get_attributes()
        if "num" in e_attrs:
            e_attrs["num"] = f'{expr_to_string(preprocess_expression_for_gammaloop(Es(e_attrs["num"])))}'
    for n in g.get_nodes():
        n_attrs = n.get_attributes()
        if "num" in n_attrs:
            n_attrs["num"] = f'{expr_to_string(preprocess_expression_for_gammaloop(Es(n_attrs["num"])))}'

with open(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "qqx_aaa_subtracted.dot"),'w') as f:
    f.write('\n\n\n'.join(str(all_graphs[g_name]) for g_name in all_ordered_graphs))

# %% [markdown]
# Now read the files and process the data for evaluation

# %%
all_graphs = {g.get_name(): {'dot': g} for g in pydot.graph_from_dot_file(os.path.join(ROOT_DIR, "dot_graphs", "loop_graphs", "qqx_aaa_subtracted.dot"))}

for (g_name, g) in all_graphs.items():

    loop_tn_expr = simplify_color(
        cook_indices(
            get_numerator(g['dot'])
            *
            get_projector(g['dot'])
        )
    )
    loop_tn_expr /= get_propagator_denominators(g['dot'])

    g['tn_expr'] = loop_tn_expr

    tn_emr = TensorNetwork(loop_tn_expr, hep_lib)
    g['tn_emr'] = TensorNetwork(loop_tn_expr, hep_lib)

    emr_replacements = get_emr_replacements(g['dot'])

    g['emr_replacements'] = emr_replacements
    tn_lmb = tn_replace_multiple(tn_emr, emr_replacements)

    g['tn_lmb'] = tn_replace_multiple(tn_emr, emr_replacements)
    tn_lmb.execute(hep_lib)

    g['expr_lmb'] = tn_lmb.result_scalar()

# %% [markdown]
# And we can now test the evaluation of the full subtracted amplitudes in the collinear and UV limits

# %%
scaling_factor = 0.1
mUV = 1.0
test_point = (kin_point, helicities)
#test_point = (kin_point_onshell, helicities_onshell)
#test_point = (kin_point2_onshell, helicities2_onshell)

#test_point = (kin_point2_onshell, [-1, 1, -1, -1, 1])

prefactor_power = 2.0
x = 0.1
for test_name, k_base, scaling, pref_power in [
    ('IR limit', [1.0*abs(x),-0.0001,0.0001,-1.0*x], 0.1, 2.0),
    ('UV limit', [100.0,100.0,100.0,100.0], 10.0, 4.0),
]:
    running_sum_non_ct = complex(0.0,0.0)
    running_sum_ct = complex(0.0,0.0)
    running_sum_non_ct_scaled = complex(0.0,0.0)
    running_sum_ct_scaled = complex(0.0,0.0)
    print('')
    print(f'Investigating {test_name}...')
    k = k_base
    scaling_factor = scaling
    prefactor_power = pref_power
    for g_name, data in all_graphs.items():

        expr = data['expr_lmb']

        # if g_name in ['ir_ct', 'uv_ir_ct']:
        #     #print(f'{g_name}:\n{"\n*".join(str(data['tn_expr']).split('*'))}')
        #     print(f'{g_name}:\n{"\n*".join(f"{lhs} -> {rhs}" for lhs, rhs in data['emr_replacements'])}')
        # if g_name in ['qqx_aaa_tri_A', 'uv_triangle_A']:
        #         print(f'{g_name}:\n{"\n*".join(str(data['tn_expr']).split('*'))}')
        #         print(f'{g_name}:\n{"\n*".join(f"{lhs} -> {rhs}" for lhs, rhs in data['emr_replacements'])}')

        eval_res = expr.evaluate_complex(
            function_map_for_evaluation(test_point[0],test_point[1],loop_mom=k, mUV = mUV),
            {}
        )
        if test_name == 'IR limit':
            scaled_k = [k[0],k[1]*scaling_factor,k[2]*scaling_factor,k[3]]
        elif test_name == 'UV limit':
            scaled_k = [k[0]*scaling_factor,k[1]*scaling_factor,k[2]*scaling_factor,k[3]*scaling_factor]
        eval_res_scaled = expr.evaluate_complex(
            function_map_for_evaluation(test_point[0],test_point[1],loop_mom=scaled_k, mUV = mUV),
            {}
        )
        fudge_factor = get_fudge_factor(g_name)
        eval_res /= fudge_factor
        eval_res_scaled /= fudge_factor

        eval_res_scaled *= scaling_factor**pref_power
        if g_name in ['ir_ct', 'uv_triangle_A', 'uv_triangle_B', 'uv_triangle_C', 'uv_bubble_A', 'uv_bubble_B', 'uv_ir_ct']:
            running_sum_ct += eval_res
            running_sum_ct_scaled += eval_res_scaled
        else:
            running_sum_non_ct += eval_res
            running_sum_non_ct_scaled += eval_res_scaled
        #print(f"{'Evaluation result for graph "'+g_name+'"':50s}: {eval_res:50.16e} | scaled: {eval_res_scaled:50.16e} | log10(abs(scaled / orig)) = {log10(abs(eval_res_scaled / eval_res)) if eval_res!=0. else float('nan'):.3f} ")
    print('-'*100)
    print(f"{'Total evaluation result from non-CT':50s}: {running_sum_non_ct:50.16e} | scaled: {running_sum_non_ct_scaled:50.16e}  | log10(abs(scaled / orig)) = {log10(abs(running_sum_non_ct_scaled / running_sum_non_ct)) if running_sum_non_ct!=0. else float('nan'):.3f} ")
    print(f"{'Total evaluation result from CT':50s}: {running_sum_ct:50.16e} | scaled: {running_sum_ct_scaled:50.16e}  | log10(abs(scaled / orig)) = {log10(abs(running_sum_ct_scaled / running_sum_ct)) if running_sum_non_ct!=0. else float('nan'):.3f} ")
    print(f"{'Total evaluation result':50s}: {running_sum_ct+running_sum_non_ct:50.16e} | scaled: {running_sum_ct_scaled+running_sum_non_ct_scaled:50.16e}  | log10(abs(scaled / orig)) = {log10(abs((running_sum_ct_scaled+running_sum_non_ct_scaled) / (running_sum_ct+running_sum_non_ct))) if running_sum_ct+running_sum_non_ct!=0. else float('nan'):.3f} ")
    # print(f"{'Ratio non-CT / CT ':50s}: {running_sum_non_ct / running_sum_ct:50.16e} | scaled: {running_sum_non_ct_scaled / running_sum_ct_scaled:50.16e}  | log10(abs(scaled / orig)) = {log10(abs((running_sum_non_ct_scaled / running_sum_ct_scaled) / (running_sum_non_ct / running_sum_ct))) if running_sum_ct_scaled!=0. and running_sum_ct!=0. and running_sum_non_ct!=0. else float('nan'):.3f} ")

    print('')

# %%

capture_output = False
print(">>>>> INPUT DOTS GENERATION SUCCESSFUL")

run_gammaloop_commands("./toml_cards/generate_qqx_aaa_tree.toml", capture_output=capture_output, start_clean=True, debug=True, state_name='./states/qqx_aaa_tree', toml_run_card=True)

print(">>>>> TREE PROCESS GENERATION SUCCESSFUL")

run_gammaloop_commands("./toml_cards/generate_qqx_aaa_euclidean.toml", capture_output=capture_output, start_clean=True, debug=True, state_name='./states/qqx_aaa_euclidean', toml_run_card=True)

print(">>>>> EUCLIDEAN PROCESS GENERATION SUCCESSFUL")

run_gammaloop_commands("./toml_cards/generate_qqx_aaa_physical.toml", capture_output=capture_output, start_clean=True, debug=True, state_name='./states/qqx_aaa_physical', toml_run_card=True)

print(">>>>> PHYSICAL PROCESS GENERATION SUCCESSFUL")

cmd = "./test_limits_qqx_aaa.py -d I -t IR UV -s ./states/qqx_aaa_euclidean"
print(f"Now running '{cmd}'...")
r = subprocess.run(cmd, capture_output=capture_output, text=True, cwd=ROOT_DIR, shell=True)
if r.returncode != 0:
    raise RuntimeError(f"Command failed with error: {r.stderr}")        

print(">>>>> LOCAL APPROACH TEST EUCLIDEAN COMPLETED")

cmd = "./test_limits_qqx_aaa.py -d I th.CT tot -t IR UV THRES -s ./states/qqx_aaa_physical"
print(f"Now running '{cmd}'...")
r = subprocess.run(cmd, capture_output=capture_output, text=True, cwd=ROOT_DIR, shell=True)
if r.returncode != 0:
    raise RuntimeError(f"Command failed with error: {r.stderr}")        

print(">>>>> LOCAL APPROACH TEST PHYSICAL COMPLETED")
run_gammaloop_commands([
    f'integrate --workspace-path ./integration_workspace --result-path ./integration_workspace/integration_results.txt -c 8 --restart',
], capture_output=capture_output, start_clean=False, debug=True, state_name='./states/qqx_aaa_physical', toml_run_card=False, turn_off_debug_logs=True)

print(">>>>> INTEGRATION COMPLETED (PHYSICAL)")

print(">>>>> ALL SUCCESSFUL")
