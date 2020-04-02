# mcrun8_26.py
#   Monte Carlo program for electrons (ADP and intervalley scattering)
#          Jan Isberg 2013-2017

import datetime                        #used in filenames
import numpy as np                     #standard math library
#from defineheader import defineheader  #defineheader.py creates a header for ToF-gui readable files
# from MCinloop6 import inloop6          #inloop6 is the fast fortran scattering code using f2py
# from MCinloop6 import gammamaxfind 
# from MCinloop6 import efieldfind  
# from MCinloop6 import emuching   
import multiprocessing as mp           #for paralell processing
import matplotlib.pyplot as plt        #for plotting
import time as tm                      #to record total execution time
from memory_profiler import profile
#from scipy.spatial import ckdtree
#from scipy.spatial import distance 

gmax = np.zeros(1)
TBU_check = np.zeros(3)
def multi_run_wrap(a):    #    
    """
    collects all arguments for inloop6 into one (technicality to make paralell processing routine poolen.map work)

    Parameters
    ----------
    a : list 
    [length m, temp K, Efield V/m, bfield T, Xid, xiu,
    tstop s, maxscats, intevally scatter 1 or 0, gmax=0, which vally 1-6]
        
    Returns
    -------
    gammamax: flaot
    gmax for the parameters i  a
    """
    
    #X_np      = np.frombuffer(Efield) # V2.0
    bins_t    = int(np.frombuffer(Efield)[-2])  #number of temporal bins to record current data in # NEEDS  FIXING, NOW I HAVE TO SPECIFY BINS AT 3 DIFFERENT PLACES
    bins_z    = int(np.frombuffer(Efield)[-3]) 
    bins_xy   = int(np.frombuffer(Efield)[-4])         
    xtotinval = np.zeros(6)
    ttotinval = np.zeros(6)
    Ettot     = np.zeros(6)
    xinttot   = np.zeros(6)
    xvtot     = np.zeros(6)
    jz2       = np.zeros(bins_t)  # FIX THIS
    Pos       = np.zeros([bins_t*4])
    Pos       = np.random.random(bins_t*4)# FIX THIS v2.0
    E_field = np.reshape(np.frombuffer(Efield)[0:bins_t*bins_z*bins_xy*bins_xy*3], (bins_t,bins_z,bins_xy,bins_xy,3),order = 'F') #np.zeros([100,100,29,29,3]) 
    #s_tid = tm.time()
    #inloop6(*a,xtotinval,ttotinval,Ettot,xinttot,xvtot,jz2,Pos,E_field) # Pos in V2.0   
    #b_tid = tm.time()
    #print('Elapsed time scatt_in: ', b_tid-s_tid)
    return([xtotinval,ttotinval,Ettot,xinttot,xvtot,jz2,Pos]) # Pos in V2.0

def gmax(a):
    """
    gmax estimates gmax for the list of paramters in a

    Parameters
    ----------
    a : list of varibles [length m, temp K, Efield V/m, bfield T, Xid, xiu,
    tstop s, maxscats, intevally scatter 1 or 0, gmax=0, which vally 1-6]
        
    Returns
    -------
    gammamax: flaot
    gmax for the parameters i  a

    """
    gammamax  = np.zeros(1)
    #gammamaxfind(*a,gammamax)
    return(gammamax)

def Efind(a):
    '''
    Take a list with the electrons pos and givs back the eletric field from 
    electrons using a fortran script 

    Parameters
    ----------
    a :list
        [t_bin,z_bin,xy_bin,extend_number,electron_pos_matrix,L/W]

    Returns
    -------
    Efield: np.array
    electrifield matrix [z_bin,xy_bin,xy_bin,3] 

    '''    
    
    Efield  = np.zeros([100,29,29,3],order='F') # hårdkodat
    #efieldfind(*a,Efield)
    return(Efield)

def electron_smoothing(bins_z,bins_t,bins_xy,n_n_s,Pos_size,Pos,N,M):
    '''
    takes a list of all postions of the electrons and put them the right
    in a matrix form 
    
    Parameters
    ----------
    bins_z : int
        DESCRIPTION.
    bins_t : int
        DESCRIPTION.
    bins_xy : int
        DESCRIPTION.
    n_n_s : int
        DESCRIPTION.
    Pos_size : int
        number of uesed postions
    Pos : np.array
        [:,4] array with all the postions of electron
    N : np.array
        [:] number of electrons in each position
    M : np.array
        [n_n_s*2+1,n_n_s*2+1,n_n_s*2+1], describes the smoothing

    Returns
    -------
    NN : np.array
    matrix which describes the simulated space with number of electrons 
    in each subspace
    
    '''    
    NN  = np.zeros([bins_t,bins_z+n_n_s*2,bins_xy+n_n_s*2,bins_xy+n_n_s*2],order='F')
    pos  = np.zeros([bins_t*bins_z*bins_xy*bins_xy,4],order='F')
    n = np.zeros([bins_t*bins_z*bins_xy*bins_xy],order='F')
    pos[:Pos_size,:] = Pos
    n[:Pos_size] = N
    M = np.asfortranarray(M)
    #emuching(n_n_s,Pos_size,pos,n,NN,M)
    return(NN)

def initProcess(share):
  global Efield
  Efield = share
@profile
def mcrun():
  intervalley_scattering = 0 #0 = intervalley scattering OFF, 1= intervalley scattering ON
  autogammacorrect = 3000 #1000 #3000 # # scatterings to correct value of gamma, =0 means autocorrect off
                
  L=300e-6              #sample thickness (or large value for "infinite" thickness)
  W= 300e-6
  maxscats=15000000     #maximum number of scattering events 
  tstop=1000.0e-9       #max time for simulation in seconds
  incycles =12       #number of electrons for each TBU value (should be multiple of 6 for equal valley distribution)

  #Natural and material constants
  me=9.10953e-31 #electron mass in kg
  kB=1.3807e-23  #Bolzmann constant J/K
  q=1.602e-19    #elementary charge As
  ml=0.98*me     #effective masses, values from N. Naka cyclotron paper. 
  mt=0.19*me
  m = np.array(( (ml, mt, mt), (ml, mt, mt), (mt, ml, mt), (mt, ml, mt), (mt, mt, ml), (mt, mt, ml) ),float) #array with masses in (x,y,z) direction for each valley
  
  # create bins
  bins_z = 100
  bins_xy = 29
  bins_t=100    #number of temporal bins to record current data in # NEEDS  FIXING, NOW I HAVE TO SPECIFY BINS AT 3 DIFFERENT PLACES
  # create time vector
  mint=-2e-9  # first travel time
  maxt= 180e-9 # maximum travel time # ?set to tstop ?
  extrabaselinezeros = 60 # extend time vector for ToFGUI
  
  deltat=(maxt-mint)/bins_t;
  tmid=np.linspace(mint+deltat/2, maxt-deltat/2, bins_t)  # time vector with "middle time" of each bin
  time = np.concatenate([np.linspace(mint+deltat/2-extrabaselinezeros*deltat, mint-deltat/2, extrabaselinezeros), tmid]) # extend time vector for ToFGUI

  #Building an array called TBU with which T=Temperature (Kelvin), B=Magnetic field (Tesla), U=Bias voltage (Volts) values to use
  count=0     #a counter for the total number of (T,B,U) multiplets
  TBU=[]
  
  for T in [78]: #two temperatures in this example
    for B in [0.0]: #no B-field in this example
      for U in [15]:  #np.logspace(-2.0, 0.0, 7): #a few voltages from 0.01 to 1 V in this example
            count += 1;
            TBU.append([T, B, U]); #In these loops we define a list with all the T,B,U values to use 
  TBU = np.array(TBU)
  # finished building the TBU array

  velmatrix2 = np.zeros((count,20),float)     # initialize velmatrix2 which has in first column the time vector
  datamatrix = np.zeros((time.size,count+1),float);   # initialise datamatrix 
  datamatrix[:,0] = time                    #write time vector to first column of datamatrix
  
  # open two files  with unique filenames including date and time
  # now = str(datetime.datetime.now())
  # now2 = now.replace(":",".") #because we cannot have colons in filenames
  # exportfilename="MC_ToF_Export_" + now2 + ".txt"             #datamatrix will be written to this file
  # velocityfilename="MC_ToF_Export_Velocity_" + now2 + ".txt"  #velmatrix will be written to this file
  # outputformat_velocity = '%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g'
  # outputformat, headerstring = defineheader(TBU)  #uses defineheader.py

  #make a place to store the electric field vector V2.0
  #E_field = np.zeros([100,100,29,29,3])
  #Siz = 7*1000000 + bins_t + bins_t*2*bins_z + np.array([bins_xy,bins_z,bins_t,B]).size
  E_field = np.zeros(bins_t*bins_z*bins_xy**2*3+4)
  E_field[bins_t*bins_z*bins_xy**2*3:] = np.array([bins_xy,bins_z,bins_t,W])
  X = mp.RawArray('d', int(bins_t*bins_z*bins_xy**2*3+4))
  X_np = np.frombuffer(X, dtype=np.float64) #.reshape(X_shape)
  np.copyto(X_np, E_field)

  
  #define one parallel process for each cpu
  workers = 2 #mp.cpu_count()    
  poolen = mp.Pool(workers, initializer= initProcess, initargs = (X,)) # V2.0


  diffrentcycles = TBU.shape[0]*6
  #Xid = 8.6*q  #Dilatation deformation potential (Joule) #this value gives good fit
  #Xiu = -0.1*q   #Uniaxial deformation potential (Joule) #this value gives good fit
  Xid = -6.57*q  #Dilatation deformation potential (Joule) #this value gives good fit
  Xiu = 16.5*q   #Uniaxial deformation potential (Joule) #this value gives good fit
  Ls = [L for x in range(diffrentcycles)]     
  TLs = [TBU[int(np.floor(x/6)),0] for x in range(diffrentcycles)]  
  Es = [TBU[int(np.floor(x/6)),2]/L for x in range(diffrentcycles)]
  Bs = [TBU[int(np.floor(x/6)),1] for x in range(diffrentcycles)]
  Xids = [Xid for x in range(diffrentcycles)]
  Xius = [Xiu for x in range(diffrentcycles)]
  tstops = [tstop for x in range(diffrentcycles)]  
  scatss = [maxscats for x in range(diffrentcycles)]
  ints = [intervalley_scattering for x in range(diffrentcycles)]
  autos = [autogammacorrect for x in range(diffrentcycles)]
  valleys = [1+(x%6) for x in range(diffrentcycles)] # repeating 0,1,2,3,4,5 list
  gmax_list = poolen.map(gmax,list(zip(*[Ls,TLs,Es,Bs,Xids,Xius,tstops,scatss,ints,autos,valleys])))

  starttime = tm.time()     #start counting cpu time
  for outloop in range(TBU.shape[0]):   #outloop does one run for each TBU value
    TL = TBU[outloop,0]
    B = TBU[outloop,1]
    U = TBU[outloop,2]
    #Xid = -6.57*q  #Dilatation deformation potential (Joule) #this value gives good fit
   # Xiu = 16.5*q   #Uniaxial deformation potential (Joule) #this value gives good fit
    print('TL= ', TL, '   B= ', B,'   U= ', U,'   Xid= ', Xid/q,'   Xiu= ', Xiu/q) #to be able to follow the progress in the command window
    jz2 = np.zeros(bins_t) # zero the current vector etc.
    xtotinval = np.zeros(6); ttotinval = np.zeros(6); Ettot = np.zeros(6); xinttot = np.zeros(6); xvtot = np.zeros(6)

    #jz2x = np.zeros(bins)
    #xtotinvalx = np.zeros(6); ttotinvalx = np.zeros(6); Ettotx = np.zeros(6)
#   Parallel call via poolen.map, first all the indata is merged into an array for all electrons (Nr=incycles)
    Ls = [L for x in range(incycles)]     
    TLs = [TL for x in range(incycles)]  
    Es = [U/L for x in range(incycles)]
    Bs = [B for x in range(incycles)]
    Xids = [Xid for x in range(incycles)]
    Xius = [Xiu for x in range(incycles)]
    tstops = [tstop for x in range(incycles)]  
    scatss = [maxscats for x in range(incycles)]
    ints = [intervalley_scattering for x in range(incycles)]
    autos = [autogammacorrect for x in range(incycles)]
    valleys = [1+(x%6) for x in range(incycles)] # repeating 0,1,2,3,4,5 list
    Gammamax = [gmax_list[int(outloop*6 +(x%6))] for x in range(incycles)] 
    #E_field2 = [np.zeros(1000000) for x in range(incycles)]
    #size2 = [L for x in range(incycles)]


#############################################################################################################
    def weird_division(n, d):
        return n / d if d else 0
    n_n_s  = int(2)
    M2 = np.zeros([n_n_s*2+1,n_n_s*2+1,n_n_s*2+1])
    M2[n_n_s,n_n_s,n_n_s] = 1
 
    EE = np.zeros([100,100,29,29,3])
    NN = np.zeros([100,100+n_n_s*2,29+n_n_s*2,29+n_n_s*2])
#########################################################################################
    #pos_raw = np.asarray(np.loadtxt('Pos_end'+'.txt'))
    #pos,N = np.unique(np.floor(np.divide(pos_raw,[L/bins_z,L/bins_xy,L/bins_xy,22e-9/bins_t])), return_index=False, return_inverse=False, return_counts=True, axis=0) # kan gå snabbare om jag gör den för varje tidsteg

    #III = (N > 0) # ta bort alla små väden
    #pos = pos[III ,:]
    #N   = N[III]
    #pos_end = pos[(pos[:,0] > bins_z-1),:]
    #N_end = N[(pos[:,0] > bins_z-1)]
    #III = (pos[:,0] < bins_z) & (0<pos[:,0]) & (pos[:,1]<15) & (pos[:,1]>-15) & (pos[:,2]<15) & (pos[:,2]>-15) # tar bort alla som är utanför 
    #pos = pos[III ,:]
    #N   = N[III]
    #ind = np.lexsort((pos[:,0],pos[:,3]), axis=0) # borde vara på alla
    #pos = pos[ind,:]
    #N   = N[ind]
    #NN = electron_smoothing(bins_z,bins_t,bins_xy,n_n_s,pos[:,0].size,pos,N,M2)

    ##E_field = np.zeros(bins_t*bins_z*bins_xy**2*3+4)
    ##E_field[bins_t*bins_z*bins_xy**2*3:] = np.array([bins_xy,bins_z,bins_t,W])
    #T_m = [bins_t for x in range(bins_t)]     
    #Z_m = [bins_z for x in range(bins_t)]  
    #XY_m = [bins_xy for x in range(bins_t)]
    #n_n_m = [n_n_s for x in range(bins_t)]
    #NN_m = [NN[x,:,:,:] for x in range(bins_t)]
    #ratio_m = [L/W for x in range(bins_t)]

    #mapper = poolen.map(Efind,list(zip(*[T_m,Z_m,XY_m,n_n_m,NN_m,ratio_m])))  
    #RRR = list(mapper)   #this is the outdata from inloop6 (fortran)
    #EE = np.stack(RRR, axis=0)

    ##X = mp.RawArray('d', int(bins_t*bins_z*bins_xy**2*3+4))
    #X_np = np.frombuffer(X, dtype=np.float64) #.reshape(X_shape)
    ##np.copyto(X_np, E_field)
    #np.copyto(X_np, np.append( np.reshape(EE*10000000/incycles, (1,-1),order='F' ) , np.array([bins_xy,bins_z,bins_t,W]) ) ) #50000000/incycles
############################################################################################################################################
    
    for run in range(3): #V2.0
        print(run)
        s_tid = tm.time()
        mapper = poolen.map(multi_run_wrap,list(zip(*[Ls,TLs,Es,Bs,Xids,Xius,tstops,scatss,ints,autos,valleys,Gammamax,])))  
        resultlist = list(mapper)   #this is the outdata from inloop6 (fortran)
        summed = [sum(x) for x in zip(*resultlist)] #outdata summed over the electrons
        xtotinval = summed[0]      #xtotinval(n)= total distance travelled in valley n (n=1..6)
        ttotinval = summed[1]      #ttotinval(n)= total time spent in valley n (n=1..6)
        Ettot     = summed[2]      #Ettot(n)= energy integrated over time in valley n (n=1..6)  ???
        xinttot   = summed[3]      #xinttot(n)= distance integrated over time in valley n (n=1..6)
        xvtot     = summed[4]      #xvtot(n)= distance times velocity integrated over time in valley n (n=1..6) ???
        jz2       = summed[5]      #jz2(n)= current  in valley n (n=1..6)  ???
        b_tid = tm.time()
        print('Elapsed time scatt: ', b_tid-s_tid) 

    

        #V2.0  sorterar ut alla positionenr för elektronerna
        Pos = np.zeros(bins_t*4) # zero the position vector etc.
        for i in range(incycles): # allt under är V2.0 
            Pos = np.array(resultlist[i][6]).reshape((bins_t,4),order='F')
            #print(np.nonzero(Pos[1:100,:]))
            #Pos = Pos[np.nonzero(Pos[:,3]),:]
            #print(Pos)
            resultlist[i][6] = Pos[np.nonzero(Pos[:,3]),:]
        pos_raw = np.array(np.hstack( (np.reshape(resultlist,incycles*7)[6::7])))
        pos,N = np.unique(np.floor(np.divide(pos_raw,[L/bins_z,L/bins_xy,L/bins_xy,22e-9/bins_t])), return_index=False, return_inverse=False, return_counts=True, axis=1) # kan gå snabbare om jag gör den för varje tidsteg
        pos = np.reshape(pos,(N.size,4))
        pos_raw = np.reshape(pos_raw,(-1,4))



        III = (N > 0) # ta bort alla små väden
        pos = pos[III ,:]
        N   = N[III]
        #print(N.size) 
        pos_end = pos[(pos[:,0] > bins_z-1),:]
        N_end = N[(pos[:,0] > bins_z-1)]
        #print(N_end.size) Är inte alltid lika med incykle?? vet inte varför

#######################################################################
        k = 0
        for i in range(pos_end[:,0].size):
            for j in range(int(pos_end[i,3]),int(bins_t)):
                k = k+1
        ## lägger till att electronerna att the end 
        pos_e = np.zeros([k,4])
        N_e = np.zeros(k)
        part_e_Nv = 0.1 # part of electrons that hits Nv center
        f = 0
        for i in range(pos_end[:,0].size):
            for j in range(int(pos_end[i,3]),int(bins_t)):
                pos_e[f,:] = [bins_z-1,pos_end[i,1],pos_end[i,2],j]
                N_e[f] = part_e_Nv*N_end[i]
                f = f+1
        pos = np.vstack((pos,pos_e))
        N = np.append(N,N_e)
#######################################################################################

        III = (pos[:,0] < bins_z) & (0<pos[:,0]) & (pos[:,1]<=int((bins_xy-1)/2)) & (pos[:,1]>=-int((bins_xy-1)/2)) & (pos[:,2]<=int((bins_xy-1)/2)) & (pos[:,2]>=-int((bins_xy-1)/2)) # tar bort alla som är utanför 
        pos = pos[III ,:]
        N   = N[III]
        ind = np.lexsort((pos[:,0],pos[:,3]), axis=0) # borde vara på alla
        pos = pos[ind,:]
        N   = N[ind]

        b_tid = tm.time()

        ratio_add = 0.02
        NN  = (1-ratio_add)*NN + ratio_add*electron_smoothing(bins_z,bins_t,bins_xy,n_n_s,pos[:,0].size,pos,N,M2)

        s_tid = tm.time()
        print('Elapsed time NN: ', s_tid-b_tid) 
################################################################################33
        T_m = [bins_t for x in range(bins_t)]     
        Z_m = [bins_z for x in range(bins_t)]  
        XY_m = [bins_xy for x in range(bins_t)]
        n_n_m = [n_n_s for x in range(bins_t)]
        NN_m = [NN[x,:,:,:] for x in range(bins_t)]
        ratio_m = [L/W for x in range(bins_t)]

        mapper = poolen.map(Efind,list(zip(*[T_m,Z_m,XY_m,n_n_m,NN_m,ratio_m])))  
        RRR = list(mapper)   #this is the outdata from inloop6 (fortran)
        EE = np.stack(RRR, axis=0)
        #print(np.mean(EE))
        b_tid = tm.time()
        print('Elapsed time EE: ', b_tid-s_tid) 
########################################################################


        plt.figure(20)
        plt.imshow(NN[40,:,:,:].sum(axis=1))
        plt.figure(21)
        plt.imshow(EE[40,:,:,:,2].sum(axis=1))

        #plt.figure(23)
        #plt.imshow(EE[:,81,:,:,2].sum(axis=0))

        
        
        X_np = np.frombuffer(X, dtype=np.float64) #.reshape(X_shape)
        np.copyto(X_np, np.append( np.reshape(EE*1e8/incycles, (1,-1),order='F' ) , np.array([bins_xy,bins_z,bins_t,W]) ) ) 
        np.savetxt('Pos_end'+str(run)+'.txt',pos_raw)

####################################################################################################################    
    # nearest = np.argsort(np.sum((Y[:,np.newaxis, :]-Y[np.newaxis,:, :]) ** 2, axis=-1), axis=1) https://jakevdp.github.io/PythonDataScienceHandbook/02.08-sorting.html 
    np.savetxt('Pos_end.txt',pos_raw)
    #np.savetxt('NN.txt',np.reshape(NN, (-1,1)))
    #np.savetxt('pos_all.txt',np.reshape(pos, (-1,1)))
    #np.savetxt('Pos_end'+'.txt',pos_raw)

########################################################################################################################
    #fancy way of summing squares (which I no longer understand, but it works (!)
    sq = lambda x:[y**2 for y in x]
    xsqtotinval  = [sum(y) for y in zip(*[sq(x) for x in list(zip(*resultlist))[0]])]
        
    Emean = Ettot/ttotinval    #average energy for each valley 
    Edrift = 1/2*m[:,2]*(xtotinval/ttotinval)**2 #average drift energy for each valley (only for x-drift)
    TC=2/(3*kB)*(Emean-Edrift) #energy in disordered DOF = 3/2*kT (equipartition theorem for ideal gas)
    
    # velmatrix2 (written to "velocityfilename") contains the following:
    #0 TL (K), 1 B (T), 2 U (V), 3 E (V/cm), 4 vt (cm/s), 5 vl (cm/s), 6 ttott (s), 7 ttotl (s),
    #8 xtott (m), 9 xtotl (m), 10 x2t (m2), 11 x2l (m2), 12 TCt (K), 13 TCl (K), 14 Xid (eV), 15 Xiu (eV)
    #16 xintl (ms), 17 xintt (ms) #18 xvl (m2), 19 xvt (m2)
    velmatrix2[outloop,:] = np.array([TL, B, U, U/L/100, 100*np.average(xtotinval[0:4]/ttotinval[0:4]),\
                          100*np.average(xtotinval[4:6]/ttotinval[4:6]), sum(ttotinval[0:4]),\
                          sum(ttotinval[4:6]), sum(xtotinval[0:4]), sum(xtotinval[4:6]),\
                          sum(xsqtotinval[0:4]), sum(xsqtotinval[4:6]),\
                          np.average(TC[0:4]), np.average(TC[4:6]), Xid/q, Xiu/q,\
                          sum(xinttot[0:4]), sum(xinttot[4:6]),sum(xvtot[0:4]), sum(xvtot[4:6])  ]) 
    
    #save current vector to datamatrix
    # datamatrix[:,outloop+1] = np.array(np.concatenate((np.zeros(extrabaselinezeros), jz2)))
    # with open(exportfilename,"w") as fid:           # open file for overwriting data
    #      fid.write(headerstring)         # write header to exportfile
    #      for row in range(datamatrix.shape[0]):
    #          print(outputformat % tuple(datamatrix[row,:]), file=fid)


    # # save velocity data
    # with open(velocityfilename,"w") as fv:           # open file for overwriting data
    #     fv.write('T \t B \t U \t E(V/cm) \t hot_v(cm/s) \t cool_v(cm/s) \t t_hot \t t_cool \t x_hot \t x_cool \t x_hot^2 \t x_cool^2 \t hot_TC \t cool_TC \t Xid (eV) \t Xiu (eV) \n')
    #     for row in range(velmatrix2.shape[0]):
    #         print(outputformat_velocity % tuple(velmatrix2[row,:]), file=fv);  # write velocitymatrix to exportfile

  #end of outloop
  
  #hooray, we are finished !
  print('Simulation finished ');
  endtime =tm.time()
  print('Elapsed time: ', endtime-starttime)
  
  # plotting some example figures
  #plt.figure(1) #plotting veloctity vs. E-field (transversal valleys blue, longitudinal valleys green)
  #plt.loglog(velmatrix2[:,3], velmatrix2[:,4], 'b^', velmatrix2[:,3], velmatrix2[:,5], 'g^')
  #plt.xlabel('Electric field (V/cm)')
  #plt.ylabel('Velocity (cm/s)')
  #plt.title('Velocity vs. E-field')
  #plt.grid(True,'both')
  
  plt.figure(23)
  plt.plot(jz2) 
  #plt.figure(2) #plotting electron temeperature vs. lattice temperature (transversal valleys blue, longitudinal valleys green)
  #plt.loglog(velmatrix2[:,0], velmatrix2[:,12], 'b^', velmatrix2[:,0], velmatrix2[:,13], 'g^')
  #plt.xlabel('Temperature (K)')
  #plt.ylabel('Carrier temperature (K)')
  #plt.title('Carrier temperature vs. Lattice temperature')
  #plt.grid(True,'both')

  plt.figure(3)
  #np.meshgrid()
  plt.hist2d(pos_raw[:,1],pos_raw[:,2], bins = 30)
  #for i in range(run+1):
  #  plt.figure(i)
  #  A = np.loadtxt('Pos_end'+str(i)+'.txt') 
  #  plt.hist2d(A[:,1],A[:,2], bins = 30)
  #plt.figure(3) #plotting electron temeperature vs. lattice temperature (transversal valleys blue, longitudinal valleys green)
  #plt.loglog(velmatrix2[:,3], velmatrix2[:,14], 'b^', velmatrix2[:,3], velmatrix2[:,15], 'g^')
  #plt.xlabel('Electric field (V/cm)')
  #plt.ylabel('Xi')
  #plt.title('Carrier temperature vs. Lattice temperature')
  #plt.grid(True,'both')


  #plt.figure(4) #plotting veloctity vs. E-field (transversal valleys blue, longitudinal valleys green)
  #plt.loglog(velmatrix2[:,3], velmatrix2[:,4], 'b^', velmatrix2[:,3], velmatrix2[:,5], 'g^')
  #plt.xlabel('Electric field (V/cm)')
  #plt.ylabel('Velocity (cm/s)')
  #plt.title('Velocity vs. E-field')
  #plt.grid(True,'both')
  
  #plt.figure(5)
  ##plt.subplot(1, 2, 2)
  #plt.loglog(velmatrix2[:,0], velmatrix2[:,12], 'b^', velmatrix2[:,0], velmatrix2[:,13], 'g^')
  ##plt.xlabel('E-field (V/cm)')
  #plt.xlabel('Lattice temperature (K)')
  #plt.ylabel('Carrier temperature (K)')
  ##plt.title('Carrier temperature vs. Xi ratio')
  #plt.grid(True,'both')

  #plt.figure(6)
  ##plt.subplot(1, 2, 1)
  ##plt.loglog(velmatrix2[:,0], velmatrix2[:,4], 'b^', velmatrix2[:,0], velmatrix2[:,5], 'g^')
  #plt.loglog(velmatrix2[:,0], 10000*q/(kB*velmatrix2[:,0])*(velmatrix2[:,10]-3/(2*incycles)*velmatrix2[:,8]**2)/(2*velmatrix2[:,6]), 'b^',\
  #           velmatrix2[:,0], 10000*q/(kB*velmatrix2[:,0])*(velmatrix2[:,11]-3/incycles*velmatrix2[:,9]**2)/(2*velmatrix2[:,7]), 'g^')
  ##plt.xlabel('E-field (V/cm)')
  #plt.xlabel('T (K)')
  #plt.ylabel('Mobility (cm2/Vs)')
  #plt.title('Mobility vs. T-lattice')
  #plt.grid(True,'both')



  plt.show()  
  poolen.close() #close the paralell workers
  poolen.join()
#end mcrun

def main():
    mcrun()
    
if __name__ == "__main__":
    main()
