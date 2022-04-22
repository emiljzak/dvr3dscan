	
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
De          = 0.3

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


def gen_params_dict(*args):

    params = {}

    params['NALF']      = NALF
    params['MAX3D']     = MAX3D
    params['De']        = De        
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

def gen_grid2D():

    print("____DVR3D calculations on a 3D grid of basis set parameters___")

    path        = os.getcwd()
    START_FILES = os.listdir(path+"/START")
    print ("The current working directory is %s" % path)

    params = gen_params_dict(   NALF,    
                                MAX3D,     
                                De,  
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


    for ir,r in enumerate(rlist):

        for iw,w in enumerate(omegalist):

            for inpnt, npnt in enumerate(npntlist):

                if npnt >= 100:
                    dirname = "r%1d"%ir+"w%1d"%iw+"N%3d"%npnt 
                elif npnt < 100:
                    dirname = "r%1d"%ir+"w%1d"%iw+"N%2d"%npnt
                else:
                    dirname = "r%1d"%ir+"w%1d"%iw+"N%1d"%npnt

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
                gen_input(params,r,w,npnt)
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
        

def gen_input(params,r,w,npnt):
    with open(inputfile,'w') as inp:
        inp.write("&PRT zrot=.true.,ztran=.true.,zlin=.true.,zpfun=.true.,ztheta=.false.,zr2r1=.false.,zembed=.false. /"+"\n") 
        inp.write("    3"+"\n")
        #inp.write("   50    0   60   50 4000 4000    1    2   50    "+"\n") #kmin = 0 implicitly, see input explanation
        inp.write( '{:5d}'.format(npnt) + "    0   60   50 4000 4000    1    2" +'{:5d}'.format(npnt) + "\n") #kmin = 0 implicitly, see input explanation
        inp.write(" N2O: JACOBI NON-SYMMETRISED COORDINATES, MASSES 14.003074 15.994915 NUCL 13.995394 15.986138"+"\n")
        inp.write("\n")
        inp.write("      14.003074          15.994915             14.003074"+"\n")
        inp.write("      14.003074          15.994915             14.003074"+"\n")                                                                      
        inp.write("      10005000.0"+"\n")
        inp.write("       4.10                0.2                 0.0085"+"\n") #r1 coordinate re, De, we
        inp.write("       0.35                0.2                 0.0305"+"\n") #r2 coordinate re, De, we


if __name__ == '__main__':
    gen_grid2D()