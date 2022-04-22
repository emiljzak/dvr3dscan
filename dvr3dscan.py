import numpy as np
import os
import shutil
import sys
import time
import subprocess


mode        = "analyze" #run
executable  = "run.sh"
inputfile   = "dvr.inp"

NALF        = 10
MAX3D       = 2000

De1         = 0.2
De2         = 0.2
NPNTfixed   = 10
omegafixed  = 0.0305
refixed     = 0.35

rmin        = 3.9
rmax        = 4.2
Nr          = 4

omegamin    = 0.0085
omegamax    = 0.0085
Nomega      = 2

NPNT_max    = 15
NPNT_min    = 10
NPNT_incr   = 1    # increment for NPNT
thr         = 1.0   # convergence threshold in cm^-1
convmode    = "rms" # "band_origin"
nlevels     = 20    # number of lowest J=0 energy levels taken in RMS calculation. Note that levels are printed in 4-columns format in DVR3D.

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
    params['convmode']      = convmode
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
                                NPNTfixed,
                                omegafixed,
                                refixed,  
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
                                convmode,     
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
                gen_input3D(params,r,w,npnt)
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

    return params,rlist,omegalist,npntlist     

def gen_input3D(params,r,w,npnt):
    with open(inputfile,'w') as inp:
        inp.write("&PRT zrot=.true.,ztran=.true.,zlin=.true.,zpfun=.true.,ztheta=.false.,zr2r1=.false.,zembed=.false. /"+"\n") 
        inp.write("    3"+"\n")
        if scan_coord == "1":
            inp.write('{:5d}'.format(params['NPNTfixed']) + "    0   60   50 4000 4000    1    2" +'{:5d}'.format(npnt) + "\n") #kmin = 0 implicitly, see input explanation
        elif scan_coord == "2":
            inp.write('{:5d}'.format(npnt) + "    0   60   50 4000 4000    1    2" +'{:5d}'.format(params['NPNTfixed']) + "\n") #kmin = 0 implicitly, see input explanation
        else:
            raise ValueError("Incorrect name of the radial coordinate")
        #inp.write("   50    0   60   50 4000 4000    1    2   50    "+"\n") #kmin = 0 implicitly, see input explanation
        inp.write(" N2O: JACOBI NON-SYMMETRISED COORDINATES, MASSES 14.003074 15.994915 NUCL 13.995394 15.986138"+"\n")
        inp.write("\n")
        inp.write("      14.003074          15.994915             14.003074"+"\n")
        inp.write("      14.003074          15.994915             14.003074"+"\n")                                                                      
        inp.write("      10005000.0"+"\n")
        if scan_coord == "1":
            inp.write('{:11.2f}'.format(r) +     '{:19.1f}'.format(De1) +  '{:23.4f}'.format(w) +"\n") #r1 coordinate re, De, we
            inp.write('{:11.2f}'.format(params['refixed']) +  '{:19.1f}'.format(De2) +   '{:23.4f}'.format(params['omegafixed']) +"\n") #r2 coordinate re, De, we
        elif scan_coord == "2":
            inp.write('{:11.2f}'.format(params['refixed']) +  '{:19.1f}'.format(De1) +   '{:23.4f}'.format(params['omegafixed']) +"\n") #r1 coordinate re, De, we
            inp.write('{:11.2f}'.format(r) +     '{:19.1f}'.format(De2) +  '{:23.4f}'.format(w) +"\n") #r2 coordinate re, De, we
        else:
            raise ValueError("Incorrect name of the radial coordinate")


def postprocess():
    """This function exctracts appropriate energy levels from DVR3D output files and calculates appropriate RMSDs
    """

    params = gen_params_dict(   NALF,    
                                MAX3D,     
                                De1,
                                De2,
                                NPNTfixed,
                                omegafixed,
                                refixed,  
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
                                convmode,     
                                nlevels )

    rlist       = list(np.linspace(params['rmin'],params['rmax'],params['Nr'],endpoint=True))
    omegalist   = list(np.linspace(params['omegamin'],params['omegamax'],params['Nomega'],endpoint=True))
    npntlist    = list(np.arange(params['NPNT_min'],params['NPNT_max'],params['NPNT_incr'],dtype=int))


    path        = os.getcwd()
    elevels    = np.zeros((int(params['nlevels']/5),5),dtype=float)
    energies    = np.zeros((len(rlist),len(omegalist),len(npntlist),params['nlevels']),dtype = float)
    rmsd        = np.zeros((len(rlist),len(omegalist),params['nlevels']),dtype = float)
    epoint      = np.zeros((4,5))

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
                os.chdir(path+"/runs/"+dirname)

                with open("dvr.out",'r') as outputfile:
                    flag = 0
                    for line in outputfile:
                        words = line.split()
                        #print(words)

                        if flag == 1 and len(words)>0:
                            energylist.append(words)

                        if len(words)>0 and words[0] == "Band" and words[1] == "origins":
                            flag        = 1
                            energylist  = []

                for ielem,elem in enumerate(energylist):
                    if ielem*5 < params['nlevels']:
                        #print(elem)
                        for ienergy,energy in enumerate(elem):
                            elevels[ielem,ienergy] = energy

                elevels_flat = elevels.reshape(-1)
                print(energies.shape)
                print(elevels_flat.shape)
                energies[ir,iw,inpnt,:] = elevels_flat
            
            Ndiff = energies.shape[2]-1 #number of differences taken wrt NPNT parameter
            
            for h in range(params['nlevels']):
                for k in range(Ndiff):
                    rmsd[ir,iw,h] += (energies[ir,iw,k+1,h]-energies[ir,iw,k,h])**2
                rmsd[ir,iw,h] /= Ndiff
                rmsd[ir,iw,h] = np.sqrt(rmsd[ir,iw,h])
                print("rmsd for r= " + str(r) + " omega= " + str(w) + " ID= " +str(h)+ " is: "  + str(rmsd[ir,iw,h]))
            mean_rmsd = np.average(rmsd,axis=2)
    print("mean rmsd for lowest " + str(params['nlevels']) + " is: " + str(mean_rmsd))
    return rmsd,mean_rmsd

if __name__ == '__main__':

    if mode == "run":
        params,rlist,omegalist,npntlist = gen_grid3D()
    elif mode == "analyze":
        postprocess()