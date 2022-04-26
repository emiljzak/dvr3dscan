import numpy as np
import os
import shutil
import sys
import time
import subprocess
import itertools
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib

mode        = "analyze" #run
executable  = "./dvr.n2o.Sch.x<dvr.inp"
inputfile   = "dvr.inp"

NALF        = 10
MAX3D       = 2000

De1         = 0.2
De2         = 0.2
NPNTfixed   = 10
omegafixed  = 0.0305
refixed     = 0.35

rmin        = 3.9
rmax        = 4.4
Nr          = 2

omegamin    = 0.0085
omegamax    = 0.0185
Nomega      = 2

NPNT_max    = 40
NPNT_min    = 10
NPNT_incr   = 10    # increment for NPNT
thr         = 1.0   # convergence threshold in cm^-1
convmode    = "rms" # "origin"
nlevels     = 20    # number of lowest J=0 energy levels taken in RMS calculation. Note that levels are printed in 4-columns format in DVR3D.

scan_coord  = "1" # which of the radial coordinates we take as active in the scan

partitions  = True # use partitioned job grid?
Nbatches    = 1 # number of batches to be executed on different machines
ibatch      = 0 # id of the present batch

Npacks      = 1 # number of packets exectuted serially on a single machine

#Note: we divide the entire job into batches and packets. Batches represent runs on independent machines, while individual packets are collections of jobs executed simulatenously on a single machine. 

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
    params['convmode']  = convmode
    params['nlevels']   = nlevels
    params['Nbatches']  = Nbatches
    params['ibatch']    = ibatch
    params['Npacks']    = Npacks

    return params

def gen_grid3D():
    """This function generates a 3D grid of basis set parameters, produces appropriate input files and submits DVR3D jobs. All jobs are executed simultaneously on a single machine. Suitable for HPCs.
    
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
    npntlist    = list(np.arange(params['NPNT_min'],params['NPNT_max'],params['NPNT_incr'],dtype=int))
    failed_list = [] #list of failed jobs

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
                proc_dvrrun = subprocess.Popen([executable], stdout=outputfile, stderr=errorfile, shell=True)
                proc_dvrrun.wait()
                
                if proc_dvrrun.returncode == 0:
                    print("job finished succesfully!")
                else:
                    print("---###job unsuccesfull###---")
                    failed_list.append([iw,ir,inpnt])
                
                time.sleep(2)
                outputfile.close()
                errorfile.close()
                os.chdir(path)

    return params,rlist,omegalist,npntlist     


def gen_grid3D_partitions():
    """This function generates a 3D grid of basis set parameters, produces appropriate input files and submits DVR3D jobs. Jobs are partitioned into batches and packets. Suitable for small local computer clusters.
    
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
                                nlevels,
                                Nbatches,
                                ibatch,
                                Npacks )

    rlist       = list(np.linspace(params['rmin'],params['rmax'],params['Nr'],endpoint=True))
    omegalist   = list(np.linspace(params['omegamin'],params['omegamax'],params['Nomega'],endpoint=True))
    npntlist    = list(np.arange(params['NPNT_min'],params['NPNT_max'],params['NPNT_incr'],dtype=int))

    G   = np.array(list(itertools.product(*[rlist, omegalist, npntlist]))) 

    failed_list = [] #list of failed jobs


    Ntotal      = len(npntlist) * len(omegalist) * len(rlist)
    print("Total number of jobs = " + str(Ntotal))
    batch_sizes = []
    #for ibatch in range(Nbatches):
    batch_size      = int(Ntotal/Nbatches)
    batch_reminder  = Ntotal%Nbatches
    print("Batch size = " + str(batch_size))
    print("Batch reminder = " + str(batch_reminder))
    print("Reconstructed total number of jobs = " + str(Nbatches*batch_size+batch_reminder))


    pack_size = int(batch_size/Npacks)
    pack_reminder = batch_size%Npacks
    print("pack size = " + str(pack_size))
    print("pack reminder = " + str(pack_reminder))
    
    #global_grid = [[[]] for i in range(Nbatches)]
    global_grid = np.empty((Nbatches,Npacks,pack_size)).tolist()
    print(np.shape(global_grid))

    #exit()
    proc_list   = [None for p in range(pack_size)]
    counter = 0
    for ib in range(Nbatches):

        for ip in range(Npacks):
            
            for k in range(pack_size):    
                global_grid[ib][ip][k] = G[counter,0:3]
                counter +=1
                print(counter)


    #exit()
    
    for jj,ipack in enumerate(global_grid[ibatch]):
        for i,ijob in enumerate(ipack):
            print("running " + str(jj) + "-th pack's job no. " + str(i))
            if ijob[2] >= 100:
                dirname = "r%2.1f"%ijob[0] +"w%5.4f"%ijob[1]+"N%3d"%ijob[2] 
            elif ijob[2]  < 100:
                dirname = "r%2.1f"%ijob[0]+"w%5.4f"%ijob[1]+"N%2d"%ijob[2] 
            else:
                dirname = "r%2.1f"%ijob[0]+"w%5.4f"%ijob[1]+"N%1d"%ijob[2] 

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
            gen_input3D(params,ijob[0],ijob[1],ijob[2])
            outputfile = open('dvr.out','w')
            outputfile.write('Generated with dvr3dscan\n')
            outputfile.flush()  

            errorfile = open('dvr.err','w')
            errorfile.write('Generated with dvr3dscan\n')
            errorfile.flush()  
            print("executing command: "+ executable)
            
            #
            print(i)
            proc_list[i] = subprocess.Popen([executable], stdout=outputfile, stderr=errorfile, shell=True)

            time.sleep(2)
            outputfile.close()
            errorfile.close()
            os.chdir(path)

            
        print(proc_list)
        print("Shape of process list: " + str(np.shape(proc_list)))
        exit_codes = [p.wait() for p in proc_list]

        for _,ijob in enumerate(ipack):

            if ijob[2] >= 100:
                dirname = "r%2.1f"%ijob[0] +"w%5.4f"%ijob[1]+"N%3d"%ijob[2] 
            elif ijob[2]  < 100:
                dirname = "r%2.1f"%ijob[0]+"w%5.4f"%ijob[1]+"N%2d"%ijob[2] 
            else:
                dirname = "r%2.1f"%ijob[0]+"w%5.4f"%ijob[1]+"N%1d"%ijob[2] 

            os.remove(path+"/runs/"+dirname+'/fort.16')
            os.remove(path+"/runs/"+dirname+'/fort.15')   
        
        print("exit_codes: " +str(exit_codes))
        failed_list.append( [i for i, e in enumerate(exit_codes) if e != 0] )
        print(failed_list)

    return params

def gen_input3D(params,r,w,npnt):
    with open(inputfile,'w') as inp:
        inp.write("&PRT zrot=.true.,ztran=.true.,zlin=.true.,zpfun=.true.,ztheta=.false.,zr2r1=.false.,zembed=.false. /"+"\n") 
        inp.write("    3"+"\n")
        if scan_coord == "1":
            inp.write('{:5d}'.format(int(params['NPNTfixed'])) + "    0   60   50 4000 4000    1    2" +'{:5d}'.format(int(npnt)) + "\n") #kmin = 0 implicitly, see input explanation
        elif scan_coord == "2":
            inp.write('{:5d}'.format(int(npnt)) + "    0   60   50 4000 4000    1    2" +'{:5d}'.format(int(params['NPNTfixed'])) + "\n") #kmin = 0 implicitly, see input explanation
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
    elevels     = np.zeros((int(params['nlevels']/5),5),dtype=float)
    energies    = np.zeros((len(rlist),len(omegalist),len(npntlist),params['nlevels']),dtype = float)
    rmsd        = np.zeros((len(rlist),len(omegalist),len(npntlist)-1,params['nlevels']),dtype = float)
    Ndiff       = energies.shape[2]-1 # number of differences taken wrt NPNT parameter
    NPNT_conv   = np.zeros((len(rlist),len(omegalist)),dtype=int)
    NPNT_opt    =np.zeros((len(rlist),len(omegalist)),dtype=int)

    for ir,r in enumerate(rlist):

        for iw,w in enumerate(omegalist):

            for inpnt, npnt in enumerate(npntlist):

                if npnt >= 100:
                    dirname = "r%2.1f"%r+"w%5.4f"%w+"N%3d"%npnt 
                elif npnt  < 100:
                    dirname = "r%2.1f"%r+"w%5.4f"%w+"N%2d"%npnt 
                else:
                    dirname = "r%2.1f"%r+"w%5.4f"%w+"N%1d"%npnt 
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

                elevels_flat            = elevels.reshape(-1)
                energies[ir,iw,inpnt,:] = elevels_flat
            
           
            if params['convmode'] == 'rms':  
                for k in range(Ndiff):
                    for h in range(params['nlevels']):
                    
                        rmsd[ir,iw,k,h] += (energies[ir,iw,k+1,h]-energies[ir,iw,k,h])**2

                    rmsd[ir,iw,k,h] /= params['nlevels']
                    rmsd[ir,iw,k,h] = np.sqrt(rmsd[ir,iw,k,h])
                    print("rmsd for r= " + str(r) + " omega= " + str(w) + " NPNT = " + str(npntlist[k]) + " is: "  + str(rmsd[ir,iw,k,h]))
                
                    mean_rmsd = np.average(rmsd,axis=3)
                    print("mean rmsd for lowest " + str(params['nlevels']) + " for NPNT = " + str(npntlist[k]) + " is: " + str(mean_rmsd[ir,iw,k]))

                    if mean_rmsd[ir,iw,k] < thr:
                        NPNT_conv[ir,iw] = np.argmin(mean_rmsd[ir,iw])
                        print("NPNT for the minimum mean rmsd is " + str(npntlist[NPNT_conv[ir,iw]]))
                        print("NPNT for rmsd below the threshold of " + str(thr) + " is " + str(npntlist[k]) )
                        NPNT_opt[ir,iw] = npntlist[NPNT_conv[ir,iw]]
                    else:
                        print("rmsd for all values of NPNT exceed threshold = " + str(thr)+ " . Assuming maximum NPNT = " + str(npntlist[-1]) + " for the map.")
                        NPNT_opt[ir,iw] = npntlist[-1]

            
            elif params['convmode'] == 'origin':
                for k in range(Ndiff):
                    rmsd[ir,iw,0] += (energies[ir,iw,k+1,0]-energies[ir,iw,k,0])**2

                rmsd[ir,iw,0] /= Ndiff
                rmsd[ir,iw,0] = np.sqrt(rmsd[ir,iw,0])
                print("rmsd for r= " + str(r) + " omega= " + str(w) + " ID= " +str(0)+ " is: "  + str(rmsd[ir,iw,0]))
    
    plot_maps(rlist,omegalist,mean_rmsd,NPNT_opt)      
    
    return rmsd,mean_rmsd

def plot_maps(rlist,omegalist,rmsd,NPNT_opt):
   
    rwgrid = np.meshgrid(rlist,omegalist,indexing="ij")
    #print(np.shape(rwgrid))
    #print(rwgrid)

    figsizex    = 50 #size of the figure on screen
    figsizey    = 50  #size of the figure on screen
    resolution  = 100  #resolution in dpi

    fig         = plt.figure(   figsize = (figsizex, figsizey), 
                                dpi     = resolution,
                                constrained_layout = True)
                                
    grid_fig    = gridspec.GridSpec(ncols=1, nrows=1, figure=fig)
    ax         = fig.add_subplot(grid_fig[0, 0], projection='rectilinear')

    cmap        = matplotlib.cm.jet #jet, cool, etc
    norm        = matplotlib.colors.Normalize(vmin=0, vmax=50)
    plot_W2D  = ax.contour(     rwgrid[0], 
                                rwgrid[1], 
                                NPNT_opt,  
                                10, 
                                cmap = 'jet')
    plt.show()

if __name__ == '__main__':

    if mode == "run":
        if partitions == True:
            params = gen_grid3D_partitions()
        else:
            params,rlist,omegalist,npntlist = gen_grid3D()
    elif mode == "analyze":
        postprocess()