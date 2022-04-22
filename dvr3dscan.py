	
    
from concurrent.futures import thread
from socket import RCVALL_MAX
from typing import NoReturn
import numpy as np
import os
import shutil
import sys
import time
import subprocess
import matplotlib.pyplot as plt



NALF        = 80
MAX3D       = 2000
De          = 0.3

rmin        = 0.25
rmax        = 0.6
Nr          = 5

omegamin    = 0.01
omegamax    = 0.2
Nomega      = 4

NPNT_max    = 100
NPNT_min    = 20
NPNT_incr   = 10    # increment for NPNT
thr         = 1.0   # convergence threshold in cm^-1
mode        = "rms" # "band_origin"
nlevels     = 10    # number of lowest J=0 energy levels taken in RMS calculation


def gen_params_dict(**kwargs):

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

    path = os.getcwd()
    print ("The current working directory is %s" % path)
    
    START_FILES = os.listdir(path+"/START")

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
    omegalist   = list(np.linspace(params['omegamin'],params['omegamax'],params['Nomega'],endpoint=True))


    for ir,r in enumerate(rlist):

        for iw,w in enumerate(omegalist):


            dirname = "beta%6.4f"%beta+"e%6.2E"%e0

            try:
                os.mkdir(outputpath+"/"+dirname)
            except OSError:
                print ("Creation of the directory %s failed" % outputpath)
            else:
                print ("Successfully created the directory %s " % outputpath)


            #generate pulse

            if state_to_state == False:

                print("state_to_state = False")
                print("Field type selected: ", field_type)
                print("Field sub-type selected: ", field_subtype)
                print("Envelope type selected: ", envelope_type)
                tres=50
                params['tres'] = tres
                # generate envelope
                fenv= np.zeros(len(params['t_richmol']))
                fenv = gp.calc_envelope(params,param_mode,envelope_type)

                # generate field
                field = np.zeros(shape=(len(params['t_richmol']),4))
                field = gp.calc_field(params,param_mode,field_type,field_subtype,fenv)
                #print(field)

                # print field
                print("saving field to file:",output)
                fl = open(output,'w')
                for elem in field:
                    fl.write("%12.4f"%elem[0]+"  %16.8e"%elem[1]+"  %16.8e"%elem[2]+" %16.8e"%elem[3]+"\n")
                fl.close()

            elif state_to_state == True:
                print("state_to_state = True")

                if path_from_file == True:
                    print("Path read from file",params['pathfile'])
                    tres=gt.gen_tres(field_type,state_to_state,path_from_file,optimize_path,params)
                    params['tres'] = tres
                    print("state_to_state = True")
                    print("Field type selected: ", field_type)
                    print("Field sub-type selected: ", field_subtype)
                    print("Envelope type selected: ", envelope_type)

                    # generate envelope
                    fenv= np.zeros(len(params['t_richmol']))
                    fenv = gp.calc_envelope(params,param_mode,envelope_type)

                    # generate field
                    field = np.zeros(shape=(len(params['t_richmol']),4))
                    field = gp.calc_field(params,param_mode,field_type,field_subtype,fenv)
                    #print(field)

                    # print field
                    print("saving field to file:",output)
                    fl = open(output,'w')
                    for elem in field:
                        fl.write("%12.4f"%elem[0]+"  %16.8e"%elem[1]+"  %16.8e"%elem[2]+" %16.8e"%elem[3]+"\n")
                    fl.close()


                elif path_from_file == False:
                    print("Path read manually")
                else:
                    print("Error: incorrect state_to_state value")
                    exit()

            else:
                print("Error: incorrect state_to_state value")
                exit()



            for f in START_FILES: #copy all files from START to the active directory
                shutil.copy2(path+"/START/"+f, outputpath+"/"+dirname)
            shutil.copy2(path+"/elfield.txt",outputpath+"/"+dirname) #copy elfield.txt to the active directory

            os.chdir(outputpath+"/"+dirname)
            subprocess.call("./master_script.sh oc-run.inp", shell=True) 
            os.chdir(path)
        