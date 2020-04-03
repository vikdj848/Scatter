# mcrun8_26.py
#   Monte Carlo program for electrons (ADP and intervalley scattering)
#          Jan Isberg 2013-2017

import datetime                        #used in filenames
import numpy as np                     #standard math library
from defineheader import defineheader  #defineheader.py creates a header for ToF-gui readable files
from MCinloop6 import inloop6          #inloop6 is the fast fortran scattering code using f2py
from MCinloop6 import gammamaxfind 
from MCinloop6 import efieldfind  
from MCinloop6 import emuching   
import multiprocessing as mp           #for paralell processing
import matplotlib.pyplot as plt        #for plotting
import time as tm                      #to record total execution time
from scipy.spatial import ckdtree
from scipy.spatial import distance 
from memory_profiler import profile
import scattervariables as sv

gmax = np.zeros(1)

@profile
def mcrun():
  '''
  Monte Carlo simulation for electrons with ADP and intervalley scattering.
  Set parameters between line 155 and 190.
  Returns
  -------
  None.

  '''

  # create bins
  bins_z  = sv.bins_z
  bins_xy = sv.bins_xy
  bins_t  = sv.bins_t    #number of temporal bins to record current data in # NEEDS  FIXING, NOW I HAVE TO SPECIFY BINS AT 3 DIFFERENT PLACES

  #make a place to store the electric field vector V2.0
  e_field = np.zeros(bins_t*bins_z*bins_xy**2*3)
  raw_array = mp.RawArray('d', int(bins_t*bins_z*bins_xy**2*3))
  raw_array_np = np.frombuffer(raw_array, dtype=np.float64) #.reshape(X_shape)
  np.copyto(raw_array_np, e_field)

  #define one parallel process for each cpu
  workers = mp.cpu_count()    
  poolen = mp.Pool(workers, initializer= initProcess, initargs = (raw_array,)) # V2.0

  # calculates gmax for each vally
  diffrent_cycles = 6
  valleys = [1+(x%6) for x in range(diffrent_cycles)] # repeating 0,1,2,3,4,5 list
  gmax_list = poolen.map(gmax,list(zip(*[valleys])))
###################################################
  start_time = tm.time()     #start counting cpu time
  print('TL= ', sv.T_L, '   B= ', sv.B,'   U= ', sv.U,'   Xid= ', sv.Xid/sv.q,'   Xiu= ', sv.Xiu/sv.q) #to be able to follow the progress in the command window
  jz2 = np.zeros(bins_t) # zero the current vector etc.
  xtotinval = np.zeros(6); 
  ttotinval = np.zeros(6); 
  Ettot = np.zeros(6);
  xinttot = np.zeros(6); 
  xvtot = np.zeros(6)

############################################################################################################
  E_field = np.zeros([100,100,29,29,3])
  number_e_matrix = np.zeros([100,100,29,29])

  for run in range(2): #V2.0
        print(run)
        run_start_time = tm.time()
        # Parallel call via poolen.map, first all the indata is merged into an array for all electrons (Nr=incycles)
        valleys = [1+(x%6) for x in range(sv.incycles)] # repeating 0,1,2,3,4,5 list
        Gammamax = [gmax_list[int(x%6)] for x in range(sv.incycles)] 
        mapper = poolen.map(multi_run_wrap,list(zip(*[valleys,Gammamax])))  
        resultlist = list(mapper)   #this is the outdata from inloop6 (fortran)
        summed = [sum(x) for x in zip(*resultlist)] #outdata summed over the electrons
        xtotinval = summed[0]      #xtotinval(n)= total distance travelled in valley n (n=1..6)
        ttotinval = summed[1]      #ttotinval(n)= total time spent in valley n (n=1..6)
        Ettot     = summed[2]      #Ettot(n)= energy integrated over time in valley n (n=1..6)  ???
        xinttot   = summed[3]      #xinttot(n)= distance integrated over time in valley n (n=1..6)
        xvtot     = summed[4]      #xvtot(n)= distance times velocity integrated over time in valley n (n=1..6) ???
        jz2       = summed[5]      #jz2(n)= current  in valley n (n=1..6)  ???
        run_end_time = tm.time()
        print('Elapsed time scatt: ', run_end_time-run_start_time) 
        #V2.0  sorterar ut alla positionenr för elektronerna

        pos_i = np.zeros(bins_t*4) # zero the position vector etc. make sure that the format of the vectors is correct and all zero values are gone
        for i in range(sv.incycles): # allt under är V2.0 
            pos_i = np.array(resultlist[i][6]).reshape((bins_t,4),order='F')
            resultlist[i][6] = pos_i[np.nonzero(pos_i[:,3]),:]
        
        pos_raw = np.array(np.hstack( (np.reshape(resultlist,sv.incycles*7)[6::7]))) 
        pos_bins = [] 
        pos,number_e = np.unique(np.floor(np.divide(pos_raw,[sv.length/bins_z,sv.length/bins_xy,sv.length/bins_xy,22e-9/bins_t])),
             return_index=False, return_inverse=False, return_counts=True, axis=1) # kan gå snabbare om jag gör den för varje tidsteg
        pos = np.reshape(pos,(number_e.size,4))
        pos_raw = np.reshape(pos_raw,(-1,4))

        ind_greater_zero = (number_e > 0) # ta bort alla små väden
        pos = pos[ind_greater_zero ,:]
        number_e   = number_e[ind_greater_zero]
        pos_end = pos[(pos[:,0] > bins_z-1),:]
        number_e_end = number_e[(pos[:,0] > bins_z-1)]
##simulats that the electrons get stuck in Nv centers
        part_e_Nv = 0.1 # part of electrons that hits Nv center
        count = 0
        k = 0
        for i in range(pos_end[:,0].size):
            for j in range(int(pos_end[i,3]),int(bins_t)):
                count += 1
        ## lägger till electroner till slutet the end 
        pos_e_Nv = np.zeros([count,4])
        number_e_Nv = np.zeros(count)
        for i in range(pos_end[:,0].size):
            for j in range(int(pos_end[i,3]),int(bins_t)):
                pos_e_Nv[k,:] = [bins_z-1,pos_end[i,1],pos_end[i,2],j]
                number_e_Nv[k] = part_e_Nv*number_e_end[i]
                k += 1
        pos = np.vstack((pos,pos_e_Nv))
        number_e = np.append(number_e,number_e_Nv)
###############################################################
        ind_inside = ((pos[:,0] < bins_z) & (0<pos[:,0]) &
            (pos[:,1]<=int((bins_xy-1)/2)) & 
            (pos[:,1]>=-int((bins_xy-1)/2)) &  
            (pos[:,2]<=int((bins_xy-1)/2)) & 
            (pos[:,2]>=-int((bins_xy-1)/2))) 
        pos = pos[ind_inside ,:]
        number_e   = number_e[ind_inside]
        ind_sort = np.lexsort((pos[:,0],pos[:,3]), axis=0)
        pos = pos[ind_sort,:]
        number_e   = number_e[ind_sort]
        b_tid = tm.time()
        ratio_add = 0.02
        number_e_matrix  = ((1-ratio_add)*number_e_matrix + 
            ratio_add*electron_smoothing(pos[:,0].size,pos,number_e))
        run_start_time = tm.time()
        print('Elapsed time number_e_matrix: ', run_start_time-b_tid) 
####################################################################
        e_pos_matrix_i = [number_e_matrix[x,:,:,:] for x in range(bins_t)]
        mapper = poolen.map(Efind,list(zip(*[e_pos_matrix_i])))  
        result_efield_pool = list(mapper)   #this is the outdata from inloop6 (fortran)
        E_field = np.stack(result_efield_pool, axis=0)
        #print(np.mean(E_field))
        b_tid = tm.time()
        print('Elapsed time E_field: ', b_tid-run_start_time) 
########################################################################
        raw_array_np = np.frombuffer(raw_array, dtype=np.float64) #.reshape(X_shape)
        np.copyto(raw_array_np, np.reshape(E_field*1e8/sv.incycles, (1,-1),order='F' ) )
        np.savetxt('Pos_end'+str(run)+'.txt',pos_raw)
############################################################################    
  np.savetxt('Pos_end.txt',pos_raw)
  #np.savetxt('number_e_matrix.txt',np.reshape(number_e_matrix, (-1,1)))
  #np.savetxt('pos_all.txt',np.reshape(pos, (-1,1)))
  #np.savetxt('Pos_end'+'.txt',pos_raw)
  plt.figure(20)
  plt.imshow(number_e_matrix[40,:,:,:].sum(axis=1))
  plt.figure(21)
  plt.imshow(E_field[40,:,:,:,2].sum(axis=1))
########################################################################################################################
  end_time =tm.time()
  print('Elapsed time: ', end_time-start_time)


  plt.show()  
  poolen.close() #close the paralell workers
  poolen.join()
#end mcrun
def multi_run_wrap(a):    #collects all arguments for inloop6 into one (technicality to make paralell processing routine poolen.map work)
    """
    collects all arguments for inloop6 into one (technicality to make paralell 
    processing routine poolen.map work)

    Parameters
    ----------
    a : list 
    [gmax, which vally 1-6]
        
    Returns
    -------
    : list
    [xtotinval,ttotinval,Ettot,xinttot,xvtot,jz2,Pos]
    """
    bins_z    = sv.bins_z
    bins_xy   = sv.bins_xy
    bins_t    = sv.bins_t          
    xtotinval = np.zeros(6)
    ttotinval = np.zeros(6)
    Ettot     = np.zeros(6)
    xinttot   = np.zeros(6)
    xvtot     = np.zeros(6)
    jz2       = np.zeros(bins_t)  # FIX THIS
    Pos       = np.zeros([bins_t*4])  # FIX THIS v2.0
    E_field = np.reshape(np.frombuffer(Efield), (bins_t,bins_z,bins_xy,bins_xy,3),order = 'F') #np.zeros([100,100,29,29,3]) #
    #s_tid = tm.time()
    inloop6(sv.length,sv.T_L,sv.U/sv.length,sv.B,sv.Xid,sv.Xiu,sv.stop_t,
        sv.max_scats,sv.intervalley_scattering,sv.auto_gamma_correct,
        *a,xtotinval,ttotinval,Ettot,xinttot,xvtot,jz2,Pos,E_field) # Pos in V2.0   
    #b_tid = tm.time()
    #print('Elapsed time scatt_in: ', b_tid-s_tid)
    return([xtotinval,ttotinval,Ettot,xinttot,xvtot,jz2,Pos]) # Pos in V2.0

def gmax(a):
    """
    gmax estimates gmax for the list of paramters in a

    Parameters
    ----------
    a : int
    which vally, 1-6
        
    Returns
    -------
    gammamax: flaot

    """
    gammamax  = np.zeros(1)
    gammamaxfind(sv.length,sv.T_L,sv.U/sv.length,sv.B,sv.Xid,sv.Xiu,sv.stop_t,
        sv.max_scats,sv.intervalley_scattering,sv.auto_gamma_correct,*a,gammamax)
    return(gammamax)

def Efind(a):
    '''
    Take a list with the electrons pos and givs back the eletric field from 
    electrons using a fortran script 

    Parameters
    ----------
    a :list
        [electron_pos_matrix]

    Returns
    -------
    Efield: np.array
    electrifield matrix [z_bin,xy_bin,xy_bin,3] 

    '''    

    Efield  = np.zeros([sv.bins_z,sv.bins_xy,sv.bins_xy,3],order='F')
    efieldfind(sv.bins_z,sv.bins_xy,sv.length/sv.width,*a,Efield)
    return(Efield)

def electron_smoothing(Pos_size,Pos,N):

    '''
    takes a list of all postions of the electrons and put them the right
    in a matrix form 
    
    Parameters
    ----------
    Pos_size : int
        number of uesed postions
    Pos : np.array
        [:,4] array with all the postions of electron
    N : np.array
        [:] number of electrons in each position


    Returns
    -------
    e_number_matrix : np.array
    matrix which describes the simulated space with number of electrons 
    in each subspace
    
    '''    

    pos  = np.zeros([sv.bins_t*sv.bins_z*sv.bins_xy*sv.bins_xy,4],order='F')
    pos[:Pos_size,:] = Pos # fortran need to know the length beforehand
    n = np.zeros([sv.bins_t*sv.bins_z*sv.bins_xy*sv.bins_xy],order='F')
    n[:Pos_size] = N
    e_number_matrix  = np.zeros([sv.bins_t,sv.bins_z,sv.bins_xy,sv.bins_xy],order='F')
    emuching(Pos_size,pos,n,e_number_matrix)
    return(e_number_matrix)

def initProcess(share):
  global Efield
  Efield = share
def main():
    mcrun()
    
if __name__ == "__main__":
    main()


#########################################################################################
    #pos_raw = np.asarray(np.loadtxt('Pos_end'+'.txt'))
    #pos,number_e = np.unique(np.floor(np.divide(pos_raw,[sv.length/bins_z,sv.length/bins_xy,sv.length/bins_xy,22e-9/bins_t])), return_index=False, return_inverse=False, return_counts=True, axis=0) # kan gå snabbare om jag gör den för varje tidsteg

    #III = (number_e > 0) # ta bort alla små väden
    #pos = pos[III ,:]
    #number_e   = number_e[III]
    #pos_end = pos[(pos[:,0] > bins_z-1),:]
    #number_e_end = number_e[(pos[:,0] > bins_z-1)]
    #III = (pos[:,0] < bins_z) & (0<pos[:,0]) & (pos[:,1]<15) & (pos[:,1]>-15) & (pos[:,2]<15) & (pos[:,2]>-15) # tar bort alla som är utanför 
    #pos = pos[III ,:]
    #number_e   = number_e[III]
    #ind = np.lexsort((pos[:,0],pos[:,3]), axis=0) # borde vara på alla
    #pos = pos[ind,:]
    #number_e   = number_e[ind]
    #number_e_matrix = electron_smoothing(bins_z,bins_t,bins_xy,extend_matrix,pos[:,0].size,pos,N,smooth_matrix)

    #T_m = [bins_t for x in range(bins_t)]     
    #Z_m = [bins_z for x in range(bins_t)]  
    #XY_m = [bins_xy for x in range(bins_t)]
    #n_n_m = [extend_matrix for x in range(bins_t)]
    #NN_m = [number_e_matrix[x,:,:,:] for x in range(bins_t)]
    #ratio_m = [length/width for x in range(bins_t)]

    #mapper = poolen.map(Efind,list(zip(*[T_m,Z_m,XY_m,n_n_m,NN_m,ratio_m])))  
    #RRR = list(mapper)   #this is the outdata from inloop6 (fortran)
    #E_field = np.stack(RRR, axis=0)

    ##X = mp.RawArray('d', int(bins_t*bins_z*bins_xy**2*3+4))
    #raw_array_np = np.frombuffer(X, dtype=np.float64) #.reshape(X_shape)
    ##np.copyto(raw_array_np, E_field)
    #np.copyto(raw_array_np, np.append( np.reshape(E_field*10000000/sv.incycles, (1,-1),order='F' ) , np.array([bins_xy,bins_z,bins_t,width]) ) ) #50000000/incycles
############################################################################################################################################
      #  #fancy way of summing squares (which I no longer understand, but it works (!)
  #  sq = lambda x:[y**2 for y in x]
  #  xsqtotinval  = [sum(y) for y in zip(*[sq(x) for x in list(zip(*resultlist))[0]])]
        
  #  Emean = Ettot/ttotinval    #average energy for each valley 
  #  Edrift = 1/2*m[:,2]*(xtotinval/ttotinval)**2 #average drift energy for each valley (only for x-drift)
  #  TC=2/(3*kB)*(Emean-Edrift) #energy in disordered DOF = 3/2*kT (equipartition theorem for ideal gas)
    
  #  # velmatrix2 (written to "velocityfilename") contains the following:
  #  #0 T_L (K), 1 B (T), 2 U (V), 3 E (V/cm), 4 vt (cm/s), 5 vl (cm/s), 6 ttott (s), 7 ttotl (s),
  #  #8 xtott (m), 9 xtotl (m), 10 x2t (m2), 11 x2l (m2), 12 TCt (K), 13 TCl (K), 14 Xid (eV), 15 Xiu (eV)
  #  #16 xintl (ms), 17 xintt (ms) #18 xvl (m2), 19 xvt (m2)
  #  velmatrix2[outloop,:] = np.array([T_L, B, U, U/length/100, 100*np.average(xtotinval[0:4]/ttotinval[0:4]),\
  #                        100*np.average(xtotinval[4:6]/ttotinval[4:6]), sum(ttotinval[0:4]),\
  #                        sum(ttotinval[4:6]), sum(xtotinval[0:4]), sum(xtotinval[4:6]),\
  #                        sum(xsqtotinval[0:4]), sum(xsqtotinval[4:6]),\
  #                        np.average(TC[0:4]), np.average(TC[4:6]), Xid/q, Xiu/q,\
  #                        sum(xinttot[0:4]), sum(xinttot[4:6]),sum(xvtot[0:4]), sum(xvtot[4:6])  ]) 
    
  #  #save current vector to datamatrix
  #  datamatrix[:,outloop+1] = np.array(np.concatenate((np.zeros(extra_baseline_zeros), jz2)))
  #  with open(exportfilename,"w") as fid:           # open file for overwriting data
  #       fid.write(headerstring)         # write header to exportfile
  #       for row in range(datamatrix.shape[0]):
  #           print(outputformat % tuple(datamatrix[row,:]), file=fid)


  #  # save velocity data
  #  with open(velocityfilename,"w") as fv:           # open file for overwriting data
  #      fv.write('T \t B \t U \t E(V/cm) \t hot_v(cm/s) \t cool_v(cm/s) \t t_hot \t t_cool \t x_hot \t x_cool \t x_hot^2 \t x_cool^2 \t hot_TC \t cool_TC \t Xid (eV) \t Xiu (eV) \n')
  #      for row in range(velmatrix2.shape[0]):
  #          print(outputformat_velocity % tuple(velmatrix2[row,:]), file=fv);  # write velocitymatrix to exportfile

  ##end of outloop
  
  ##hooray, we are finished !
  #print('Simulation finished ');
  #endtime =tm.time()
  #print('Elapsed time: ', endtime-start_time)
  
  ## plotting some example figures
  ##plt.figure(1) #plotting veloctity vs. E-field (transversal valleys blue, longitudinal valleys green)
  ##plt.loglog(velmatrix2[:,3], velmatrix2[:,4], 'b^', velmatrix2[:,3], velmatrix2[:,5], 'g^')
  ##plt.xlabel('Electric field (V/cm)')
  ##plt.ylabel('Velocity (cm/s)')
  ##plt.title('Velocity vs. E-field')
  ##plt.grid(True,'both')
  
  #plt.figure(23)
  #plt.plot(jz2) 
  ##plt.figure(2) #plotting electron temeperature vs. lattice temperature (transversal valleys blue, longitudinal valleys green)
  ##plt.loglog(velmatrix2[:,0], velmatrix2[:,12], 'b^', velmatrix2[:,0], velmatrix2[:,13], 'g^')
  ##plt.xlabel('Temperature (K)')
  ##plt.ylabel('Carrier temperature (K)')
  ##plt.title('Carrier temperature vs. Lattice temperature')
  ##plt.grid(True,'both')

  #plt.figure(3)
  ##np.meshgrid()
  #plt.hist2d(pos_raw[:,1],pos_raw[:,2], bins = 30)
  ##for i in range(run+1):
  ##  plt.figure(i)
  ##  A = np.loadtxt('Pos_end'+str(i)+'.txt') 
  ##  plt.hist2d(A[:,1],A[:,2], bins = 30)
###################################################

  #velmatrix2 = np.zeros((count,20),float)     # initialize velmatrix2 which has in first column the time vector
  #datamatrix = np.zeros((time.size,count+1),float);   # initialise datamatrix 
  #datamatrix[:,0] = time                    #write time vector to first column of datamatrix
  
  # open two files  with unique filenames including date and time
  #now = str(datetime.datetime.now())
  #now2 = now.replace(":",".") #because we cannot have colons in filenames
  #exportfilename="MC_ToF_Export_" + now2 + ".txt"             #datamatrix will be written to this file
  #velocityfilename="MC_ToF_Export_Velocity_" + now2 + ".txt"  #velmatrix will be written to this file
  #outputformat_velocity = '%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g\t%5.5g'
  #outputformat, headerstring = defineheader(T_B_U)  #uses defineheader.py

###########################################################################

  #extra_baseline_zeros = 60 # extend time vector for ToFGUI
  
  #delta_t=(max_t-min_t)/bins_t;
  #mid_t=np.linspace(min_t+delta_t/2, max_t-delta_t/2, bins_t)  # time vector with "middle time" of each bin
  #time = np.concatenate([np.linspace(min_t+delta_t/2-extra_baseline_zeros*delta_t, min_t-delta_t/2, extra_baseline_zeros), mid_t]) # extend time vector for ToFGUI
###############################################333