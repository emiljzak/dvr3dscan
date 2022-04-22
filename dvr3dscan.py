	
from distutils.util import run_2to3
import numpy as np
import os
import shutil
import sys
import time
import subprocess
import matplotlib.pyplot as plt


executable  = "run.sh"
inputfile   = "dvr.inp"

NALF        = 80
MAX3D       = 2000

De1         = 0.2
De2         = 0.2
NPNTfixed   = 50
omegafixed  = 0.03
refixed     = 0.35


rmin        = 0.25
rmax        = 0.6
Nr          = 1

omegamin    = 0.01
omegamax    = 0.2
Nomega      = 1

NPNT_max    = 30
NPNT_min    = 20
NPNT_incr   = 10    # increment for NPNT
thr         = 1.0   # convergence threshold in cm^-1
mode        = "rms" # "band_origin"
nlevels     = 10    # number of lowest J=0 energy levels taken in RMS calculation

scan_coord = "1" # which of the radial coordinates we take as active in the scan

def gen_params_dict(*args):

    params = {}

    params['NALF']      = NALF
    params['MAX3D']     = MAX3D
    params['De1']       = De1
    params['De2']       = De2
    params['NPNTfixed'] = NPNTfixed
    params['omegafixed']= omegafixed
    params['refixed']   = refixed
    params['rmin']      = rmin
    params['rmax']      = rmax
    params['Nr']        = Nr
    params['omegamin']  = omegamin
    params['omegamax']  = omegamax
    params['Nomega']    = Nomega
    params['NPNT_min']  = NPNT_min
    params['NPNT_max']  = NPNT_max
    params['NPNT_incr'] = NPNT_incr
    params['thr']       = thr
    params['mode']      = mode
    params['nlevels']   = nlevels


    return params

def gen_grid3D():
    """This function generates a 3D grid of basis set parameters, produces appropriate input files and submits DVR3D jobs.
    
        Note: for non-symmetric triatomic molecules (N2O) two separate sets of radial basis set parameters must be used. A 6D grid might be needed to find optimal radial basis parameters. 
    """
    print("____DVR3D calculations on a 3D grid of basis set parameters___")

    path        = os.getcwd()
    START_FILES = os.listdir(path+"/START")
    print ("The current working directory is %s" % path)

    params = gen_params_dict(   NALF,    
                                MAX3D,     
                                De1,
                                De2,
                                NPNT1,
                                omega2,
                                re2,  
                                rmin,   
                                rmax,   
                                Nr,     
                                omegamin,
                                omegamax,    
                                Nomega,   
                                NPNT_max,   
                                NPNT_min,   
                                NPNT_incr,
                                thr,         
                                mode,     
                                nlevels )

    rlist       = list(np.linspace(params['rmin'],params['rmax'],params['Nr'],endpoint=True))
    omegalist   = list(np.linspace(params['omegamin'],params['omegamax'],params['Nomega'],endpoint=True))
    npntlist   = list(np.arange(params['NPNT_min'],params['NPNT_max'],params['NPNT_incr'],dtype=int))


    for ir1,r1 in enumerate(rlist):

        for iw1,w1 in enumerate(omegalist):

            for inpnt1, npnt1 in enumerate(npntlist):

                if npnt1 >= 100:
                    dirname = "r%1d"%ir1+"w%1d"%iw1+"N%3d"%npnt1
                elif npnt1 < 100:
                    dirname = "r%1d"%ir1+"w%1d"%iw1+"N%2d"%npnt1
                else:
                    dirname = "r%1d"%ir1+"w%1d"%iw1+"N%1d"%npnt1

                print("dirname: " + dirname)

                try:
                    os.mkdir(path+"/runs/"+dirname)
                except OSError:
                    print ("Creation of the directory %s failed" % dirname)
                else:
                    print ("Successfully created the directory %s " % dirname)

        
                for f in START_FILES: #copy all files from START to the active directory
                    shutil.copy2(path+"/START/"+f, path+"/runs/"+dirname)


                os.chdir(path+"/runs/"+dirname)
                gen_input3D(params,r1,w1,npnt1)
                exit()
                outputfile = open('dvr.out','w')
                outputfile.write('Generated with dvr3dscan\n')
                outputfile.flush()  

                errorfile = open('dvr.err','w')
                errorfile.write('Generated with dvr3dscan\n')
                errorfile.flush()  
                print("executing command: "+ executable)
                proc_dvrrun = subprocess.Popen(["./run.sh"], stdout=outputfile, stderr=errorfile, shell=True)

                #p2 = subprocess.Popen(["ls", "-l"], stdout=subprocess.PIPE)
                #stdout, stderr = p2.communicate()
                #print(stdout)
                #print(stderr)
                os.chdir(path)
        

def gen_input3D(params,r2,w2,npnt2):
    with open(inputfile,'w') as inp:
        inp.write("&PRT zrot=.true.,ztran=.true.,zlin=.true.,zpfun=.true.,ztheta=.false.,zr2r1=.false.,zembed=.false. /"+"\n") 
        inp.write("    3"+"\n")
        #inp.write("   50    0   60   50 4000 4000    1    2   50    "+"\n") #kmin = 0 implicitly, see input explanation
        inp.write('{:5d}'.format(npnt2) + "    0   60   50 4000 4000    1    2" +'{:5d}'.format(params['NPNT1']) + "\n") #kmin = 0 implicitly, see input explanation
        inp.write(" N2O: JACOBI NON-SYMMETRISED COORDINATES, MASSES 14.003074 15.994915 NUCL 13.995394 15.986138"+"\n")
        inp.write("\n")
        inp.write("      14.003074          15.994915             14.003074"+"\n")
        inp.write("      14.003074          15.994915             14.003074"+"\n")                                                                      
        inp.write("      10005000.0"+"\n")
        if scan_coord == "1":
            inp.write('{:11.2f}'.format(r1) +     '{:19.1f}'.format(De1) +  '{:23.4f}'.format(w1) +"\n") #r1 coordinate re, De, we
            inp.write('{:11.2f}'.format(params['re2']) +  '{:19.1f}'.format(De2) +   '{:23.4f}'.format(params['omega2']) +"\n") #r2 coordinate re, De, we
        elif scan_coord == "2":

        else:
            raise ValueError("Incorrect name of the radial coordinate")
       


if __name__ == '__main__':
    gen_grid3D()