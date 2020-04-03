subroutine inloop6(L,TL,Ez,B,Xid,Xiu,tmax,scats,intervalley_scattering,autogammacorrect,val,gmax,xtotinval,ttotinval,Ettot,xinttot,xvtot,jz2,Pos,E_field) BIND(C,NAME="inloop6_")
	implicit none   !All variables need to be properly declared
	integer, parameter :: bins=100 !1000 Mus be same as in python code
    integer, parameter  :: dbl=selected_real_kind(p=15)
	real (kind=dbl) :: L, TL, Ez, B, Xid, Xiu, tmax  ! L sample thickness, TL lattice temperature, Ez electric field along z, B-field along y,Xid & Xiu deformation poteneials, tmax maximum time
!f2py	intent(in) :: L, TL, Ez, B, Xid, Xiu,  tmax
    real (kind=dbl),INTENT(INOUT) :: gmax
    integer :: scats,intervalley_scattering, autogammacorrect  !scats max number of scatterings, intervalley_scattering 0 off 1 on, autogammacorrect number of loops to correct gamma
!f2py		intent(in) :: scats,intervalley_scattering, autogammacorrect
    integer  :: val !the occupied valley 1-6
!f2py		intent(inout) :: val
    real (kind=dbl)  :: t1,t2
    
    real (kind=dbl), dimension(6) :: xtotinval, ttotinval, Ettot, xinttot, xvtot 
    real (kind=4),INTENT(IN), dimension(100,100,-14:14,-14:14,3) :: E_field
    !real (kind=dbl),INTENT(IN), dimension(200,100) :: z_ind
    !real (kind=dbl), INTENT(IN),dimension(100) :: t_ind
    !integer,INTENT(IN), dimension(1000000,4) :: position

	!        xtotinval: z-distance travelled in valley n (n=1-6)
	!        ttotinval: time spent in valley n
    !        Ettot: energy integrated over time for valley n 
    !        xinttot: z-coordinate integrated over time for valley n
    !        xvtot: z-coordinate times z-velocity integrated over time for valley n
!f2py		    intent(inout) :: xtotinval, ttotinval, Ettot, xinttot, xvtot
    real (kind=dbl), dimension(bins) ::  jz2 !current in z-direction (binned in "bins" bins)
    real (kind=dbl), dimension(bins*4) ::  Pos !position of e at time of the end of bin V2.0
!f2py			intent(inout) ::  jz2

    
    !Natural and material constants
    real (kind=dbl),parameter :: pi=3.1415926535
	real (kind=dbl),parameter :: q=1.602e-19  	!electron charge in C
	real (kind=dbl),parameter :: me=9.10953e-31 	!electron mass in kg
	real (kind=dbl),parameter :: kB=1.3807e-23 	!Bolzmann constant J/K
	real (kind=dbl),parameter :: hbar=1.05459e-34 	!Plancks constant Js
    
	!real (kind=dbl),parameter :: rho=3.51e3 		!diamond density in kg/m3
	!real (kind=dbl),parameter :: vt=12.82e3 		!TA sound velocity in m/s 
	!real (kind=dbl),parameter :: vl=17.52e3 		!LA sound velocity in m/s 
	!real (kind=dbl),parameter :: ml=1.56*me 		!Values from N. Naka stability paper. 
	!real (kind=dbl),parameter :: mt=0.28*me
 !   real (kind=dbl),parameter,dimension(4) :: maxg  = (/ 0.735, 0.1708, 9.6806, 1.6398  /)		       !maximum in g (for diamond, they depend on ml & mt) see AGMplot.m
    real (kind=dbl),parameter :: W=300E-6 		!diamond density in kg/m3
	real (kind=dbl),parameter :: rho=3.51e3 		!diamond density in kg/m3
	real (kind=dbl),parameter :: vt=12.82e3 		!TA sound velocity in m/s 
	real (kind=dbl),parameter :: vl=17.52e3 		!LA sound velocity in m/s 
	real (kind=dbl),parameter :: ml=1.56*me 		!Values from N. Naka stability paper. 
	real (kind=dbl),parameter :: mt=0.28*me
    real (kind=dbl),parameter,dimension(4) :: maxg  = (/ 0.735, 0.1708, 9.6806, 1.6398  /)		       !maximum in g (for diamond, they depend on ml & mt) see AGMplot.m
    
    
	real (kind=dbl),parameter,dimension(6,3) :: m = (/ml,ml,mt,mt,mt,mt,mt,mt,ml,ml,mt,mt, mt,mt,mt,mt,ml,ml/) !  array with mass values for each valley 
    real (kind=dbl),parameter ::  mdos =(ml*mt*mt)**(1/3.0) !DOS effective mass in kg
	real (kind=dbl),parameter ::  mavg = sqrt(ml*mt)/mdos   !NO longer used, remove
    integer,parameter ::  Zf = 4 !valleys to f-scatter into
    
    real (kind=dbl), parameter :: DOevcm=4.0e8   ! Fitting parameter optical scattering - order of magnitude 1e9-1e10 (eV/cm) (intervally)
    real (kind=dbl),parameter :: DI = DOevcm*q*100 !used in intervalley (eV/m)
    real (kind=dbl),parameter :: hbaromegaf = 110E-3*q !Threshold for f-scattering LA phonons (TO has higher threshold)
    real (kind=dbl),parameter :: hbaromegag = 165E-3*q !Threshold for g-scattering (only LO allowed, see Nava "Electron Effective masses..."
    !real (kind=dbl),parameter :: DA = DAev*q
	real (kind=dbl),parameter :: emmax = 0.6476102      !maximum for emission Q^2 *Bose-Einstein function
	real (kind=dbl),parameter :: C6    = mdos/(2*hbar**2*rho) !a useful parameter
	real (kind=dbl),dimension(4) :: g

    
    real (kind=dbl), dimension(6,3,3) :: Tm, Tminv  !Heering-Voigt matrix and its inverse (ISSUE: Tm does not seem to be used, Tminv is used in intervalley scattering but the definition has disappeared in this verson)
	integer, dimension(6,3,3) :: Pm  !six 3x3 permutation (rotation) matrices going from (x,y,z) to "(t,t,l)" 
    integer :: i,vali,errors,bin,valdum,j,ind, binz, binx, biny,P1,P2,P3, binprev,k !various looping variables
    integer, dimension(8) :: scatstat !collects scattering statistics (for diagnostics only)
    real (kind=dbl),parameter :: mint = -2e-9           ! v2.0   before -30.0e-9      ! first travel time
    real (kind=dbl),parameter :: maxt = 20.0e-9      ! maximum travel time
    real (kind=dbl),parameter :: deltat =(maxt-mint)/bins; !width of the bins
        
	real (kind=dbl):: kT !kB*TL
    real (kind=dbl)  :: tmp,t,maxXi2,C1,C3,DA,tprev, t_at_bin
    real (kind=dbl),dimension(3) :: x,p,xprev,E, x_at_bin,pprev
    integer, dimension(3) :: P23
    real (kind=dbl):: maxen, gamma_guess, gammamax, gammamax_1, gammamax_2, toterrmargin, toterr, tau0, tau1, tau2, tc, kmax_low, kmax_mid, test
    real (kind=dbl):: En, errmargin
    integer,dimension(bins) :: N
	real (kind=dbl),dimension(bins) :: Nx,Nt,Nxt,Ntt
    real (kind=dbl),dimension(6) :: xinval,tinval,Et,xint,xv
    
	!call cpu_time ( t1 )
    
    Tm = 0.0   !unused ??
	Tminv = 0.0  ! Tminv definition missing
	Pm = 0
    do vali = 1,6
    	do i = 1,3
        	Pm(vali,i,1+mod(i-1+(vali+1)/2,3)  )  = 1 !Axis Permutation Matrices, vali=1,2:120deg rotation around [1,1,1], vali=3,4: -120deg, vali=5,6: identity matrix
		end do
	end do
    
    !print *, E_field(2,1,-14,-13,1)
    !print *, E_field(2,1,-14,-13,2)
    !print *, E_field(2,1,-14,-13,3)
	!can be used kept for diagnostics, Warning: gives lots of output
	!print *,'indata',L,TL,Ez,B,Xid/q,Xiu/q,intervalley_scattering,autogammacorrect,val,xtotinval,ttotinval,Ettot,jz2

    kT = kB*TL
	DA=max(abs(Xid),abs(Xid+Xiu)) ! DA now only used in C1 and for initial gamma guess
	maxXi2=max(Xid**2,(Xid+Xiu)**2)
	g(1)=kT**2/(hbar**2*vl**2) *emmax*maxXi2*maxg(1)  !SUSPICIOUS factor 2 REMOVED
	g(2)=kT**2/(hbar**2*vt**2) *emmax*Xiu**2*maxg(2)  !SUSPICIOUS factor 2 REMOVED
    C1 = mdos**1.5*sqrt(2.0)*DA**2*kT /(pi*hbar**4*rho*((vl+2*vt)/3)**2)   !(J^-0.5 s^-1) !C1 is the constant in Jacoboni 2.4.19, now only used for initial gamma guess
    !C2 = mdos**0.5*DA**2*kT**3 /(2**2.5*pi*hbar**4*rho*vavg**4) !(J^+0.5 s^-1) !no longer used
    C3 = mdos**1.5*DI**2/(sqrt(2.0)*pi*hbar**2*rho) !intervalley Jacoboni (3.74) (J^-0.5 s^-1)
	!C5 = mdos*DA**2/(rho*vavg*hbar**2) !(m^2 s^-1) !no longer used
         
    scatstat = 0

    call random_seed
    
    gammamax = gmax
    kmax_low = sqrt( (gammamax - 2*C6/vl*g(1) - 2*C6/vt*g(2) )/( C6/vl*maxXi2*maxg(3) + C6/vt*Xiu**2*maxg(4) ) )*0.2
    kmax_mid = kmax_low/0.2*0.4
    gammamax_1 = 2*C6/vl*g(1) + 2*C6/vt*g(2) + (C6/vl*maxXi2*maxg(3) + C6/vt*Xiu**2*maxg(4))*kmax_mid**2
    gammamax_2 = 2*C6/vl*g(1) + 2*C6/vt*g(2) + (C6/vl*maxXi2*maxg(3) + C6/vt*Xiu**2*maxg(4))*kmax_low**2
    tau0 = 1/gammamax !(re)define tau0, Jacoboni ch 3.4
    tau1 = 1/gammamax_1
    tau2 = 1/gammamax_2
    
    !print *, tau0, tau1
    
    toterrmargin = 0; errors = 0;
    x = 0.0
    x(3) = 3e-6
    p = 0.0 !initial position and momentum
    N = 0; Nx = 0.0; Nt = 0.0; Nxt = 0.0; Ntt = 0.0; jz2 = 0.0
    xinval = 0.0; tinval = 0.0; Et = 0.0; xint=0.0; xv = 0.0
    xtotinval = 0.0; ttotinval = 0.0; Ettot = 0.0; xinttot = 0.0; xvtot = 0.0
    t=0.0;		!alternatively t=3/2.4*1e-9*np.random.normal() gives Gaussian distributed start time (emulating 3 ns laser pulse) ***Python remnant revrite to Fortran before use
				!x(3)=-3e-6*log(np.random.random()) !Exponentially distributed z-start emulating 3um penetration depth ***Python remnant revrite to Fortran before use
    i=0			!counter incremented for each scattering 
    bin = 0
    binprev = max(1,min(1+int( (t-mint)/deltat ), bins))
    E(1) = 0
    E(2) = 0
    E(3) = Ez
    
    do
      	if ((x(3) > L) .or. (i > scats) .or. (t > tmax) ) exit ! .or. (x(3)<0) .or. (x(1)**2+x(2)**2>L**2) ) exit
    	i = i + 1
        xprev = x		! keep current x coordinate
		tprev = t		! keep current time
        pprev = p
        !print *, E
        call Scatter6(TL,E,B,t,x,p,val,En,errmargin) !perform one scattering event, updating t,x,p,val,En,errmargin
        toterrmargin = max(toterrmargin,errmargin) !a max function for the error margin, for diagnostics only
        if (errmargin > 1) errors=errors+1			!records the number of error events (gamma > gammamax), for diagnostics only
		
        
		bin = max(1,min(1+int( (t-mint)/deltat ), bins))  ! find the right temporal bin. Results outside container are put in first or last bin.
        !if ( bin .NE. max(1,min(1+int( (tprev-mint)/deltat ), bins)) ) then
        !    do i = 2, 5
                !print *, position(i,2)+1
         !       print *,position(1,1)
                !P23 = position(i+1,(/1,2,3/))
              !  E_field_matrix(1,i,1,2) = E_field(i,2)
                !E_field_matrix(position(i,1)+1,position(i,2)+1,position(i,3)+1,3) = E_field(i,3)
           ! end do
        !endif
        
        
		N(bin) = N(bin) + 1         ! N = number of events in bin
		Nx(bin) = Nx(bin) + x(3)	! Nx = total z-distance traversed in bin
		Nt(bin) = Nt(bin) + t		 ! used in current computation
		Nxt(bin) = Nxt(bin) + x(3)*t ! used in current computation
		Ntt(bin) = Ntt(bin) + t**2   ! used in current computation
        
        if( ABS(x(1))<W/2 .or. ABS(x(2))<W/2 ) then
            binz = int(min(max(x(3)/L*100,1.0),100.0))
            binx = int(min(max(x(1)/W*28,-14.0),14.0))
            biny = int(min(max(x(2)/W*28,-14.0),14.0))
        endif 
        !print *, max(1,min(1+int( (x(3)/L*1000) ), bins))
        !ind =  int(t_ind(bin+1)+1) +int(1+z_ind(bin*2+1,min(int(1+x(3)/L*1000),1000)) ) !varför +1 kolla att detta ger samma e_fält från båda.
        !E = Ez + E_field(int(t_ind(bin)+z_ind(bin,int(min(max(x(3)/L*200,1.0),200.0)) )+1)) 
        
        
        t_at_bin = (bin*deltat+mint - t)/deltat ! del av tiden som är i förgånde bin
        E(1) = (1-t_at_bin)*E_field(bin,binz,binx,biny,2) + t_at_bin*E_field(bin-1,binz,binx,biny,2)
        E(2) = (1-t_at_bin)*E_field(bin,binz,binx,biny,3) +t_at_bin*E_field(bin-1,binz,binx,biny,3)
        E(3) = (1-t_at_bin)*E_field(bin,binz,binx,biny,1) + t_at_bin*E_field(bin-1,binz,binx,biny,1) + Ez
        
        
        if (binprev .NE. bin) then
            do k = 0, bin-binprev-1
                t_at_bin = (binprev+k)*deltat+mint
                x_at_bin = xprev + pprev(:)/m(val,:)*(t_at_bin-tprev)+ E*q/(2*m(val,:))*(t_at_bin-tprev)**2
                !print *, x(3), x_at_bin(3)!, t_at_bin, t
                !print *, bin
                Pos(bin) = x_at_bin(3)           ! take the position of the at the end of t bin. borde fixas till (V2.0)
                Pos(bin+bins) = x_at_bin(1)           ! take the time of the at the end of t bin. borde fixas till (V2.0)
                Pos(bin+bins*2) = x_at_bin(2)           ! take the time of the at the end of t bin. borde fixas till (V2.0)
                Pos(bin+bins*3) = t_at_bin           ! take the time of the at the end of t bin. borde fixas till (V2.0)
                binprev = bin
            enddo
        endif
        !    binprev = bin
        !end if
        !if (x(3) == 0) print *, x(3) 
        !print *, Pos
        
		xinval(val) = xinval(val) + x(3) - xprev(3)		!traversed distance in a given valley
		tinval(val) = tinval(val) + t - tprev			!time spent in a given valley
		Et(val) =  Et(val) + En*(t-tprev)				!energy integrated over time in a given valley
		xint(val) = xint(val) + x(3)*(t-tprev)			!z-distance energy integrated over time in a given valley
		xv(val) =  xv(val) + x(3)*p(3)/m(val,3)*(t-tprev)   !z-distance * velocity integrated over time in a given valley
    end do
  !  print *, x
 !      print *,'Test error : ', toterr, '  Run margin: ', toterrmargin,' run errors: ', errors,'  Valley: ', val !for diagnostics
    xtotinval = xtotinval + xinval
    ttotinval = ttotinval + tinval
    Ettot = Ettot + Et
	xinttot = xinttot + xint
	xvtot = xvtot + xv
    !save last pos
    !print *, x(3)
    Pos(bin) = x(3)           ! take the position of the at the end of t bin. borde fixas till (V2.0)
    Pos(bin+bins) = x(1)           ! take the time of the at the end of t bin. borde fixas till (V2.0)
    Pos(bin+bins*2) = x(2)           ! take the time of the at the end of t bin. borde fixas till (V2.0)
    Pos(bin+bins*3) = t           ! take the time of the at the end of t bin. borde fixas till (V2.0)
 !     print *,scats,i,tau0   !for diagnostics
 !	   print *, scatstat(1:4) !for diagnostics
     
	! this loop computes the total current (summed for all electrons) as a function of time (in bins)
    do i = 1,bins  !do for each time interval (bin)
          if (N(i) > 1) then  !make sure that there actually has been an event in this bin, to avoid 0/0 error in next line for empty bins
              	jz2(i) = jz2(i) + q/L*(Nxt(i)-Nx(i)*Nt(i)/N(i))/(Ntt(i)-Nt(i)**2/N(i))  !Current (including displacement current) in Synchronous ensemble !NEEDS CITATION!
          end if
	end do
    
	! remainder of the file is the definition of subroutine Scatter6, which drifts one electron in the field and performs one scattering event
    !call cpu_time ( t2 )
    !write ( *, * ) 'Elapsed CPU time = ', t2 - t1
    contains
    include 'scatter6sub.f90'
end subroutine inloop6

subroutine gammamaxfind(L,TL,Ez,B,Xid,Xiu,tmax,scats,intervalley_scattering,autogammacorrect,val,gmax) BIND(C,NAME="gammamaxfind_")
	implicit none   !All variables need to be properly declared
	integer, parameter :: bins=100 !1000 Mus be same as in python code
    integer, parameter  :: dbl=selected_real_kind(p=15)
	real (kind=dbl) :: L, TL, Ez, B, Xid, Xiu, tmax  ! L sample thickness, TL lattice temperature, Ez electric field along z, B-field along y,Xid & Xiu deformation poteneials, tmax maximum time
!f2py	intent(in) :: L, TL, Ez, B, Xid, Xiu,  tmax
    real (kind=dbl),INTENT(INOUT) :: gmax
    integer :: scats,intervalley_scattering, autogammacorrect  !scats max number of scatterings, intervalley_scattering 0 off 1 on, autogammacorrect number of loops to correct gamma
!f2py		intent(in) :: scats,intervalley_scattering, autogammacorrect
    integer  :: val !the occupied valley 1-6
!f2py		intent(inout) :: val

	

    !Natural and material constants
    real (kind=dbl),parameter :: pi=3.1415926535
	real (kind=dbl),parameter :: q=1.602e-19  	!electron charge in C
	real (kind=dbl),parameter :: me=9.10953e-31 	!electron mass in kg
	real (kind=dbl),parameter :: kB=1.3807e-23 	!Bolzmann constant J/K
	real (kind=dbl),parameter :: hbar=1.05459e-34 	!Plancks constant Js
    
    
    !real (kind=dbl),parameter :: rho=3.51e3 		!diamond density in kg/m3
	!real (kind=dbl),parameter :: vt=12.82e3 		!TA sound velocity in m/s 
	!real (kind=dbl),parameter :: vl=17.52e3 		!LA sound velocity in m/s 
	!real (kind=dbl),parameter :: ml=1.56*me 		!Values from N. Naka stability paper. 
	!real (kind=dbl),parameter :: mt=0.28*me
 !   real (kind=dbl),parameter,dimension(4) :: maxg  = (/ 0.735, 0.1708, 9.6806, 1.6398  /)		       !maximum in g (for diamond, they depend on ml & mt) see AGMplot.m
    
    
	real (kind=dbl),parameter :: rho=2.33e3 		!silicon density in kg/m3
	real (kind=dbl),parameter :: vt=5.843e3 		!TA sound velocity in m/s 
	real (kind=dbl),parameter :: vl=8.433e3 		!LA sound velocity in m/s 
	real (kind=dbl),parameter :: ml=0.98*me 		!Values from N. Naka stability paper. 
	real (kind=dbl),parameter :: mt=0.19*me
    real (kind=dbl),parameter,dimension(4) :: maxg  = (/ 0.732, 0.172, 9.214, 1.626  /)		       !maximum in g (for diamond, they depend on ml & mt) see AGMplot.m
    
    
	real (kind=dbl),parameter,dimension(6,3) :: m = (/ml,ml,mt,mt,mt,mt,mt,mt,ml,ml,mt,mt, mt,mt,mt,mt,ml,ml/) !  array with mass values for each valley 
    integer,parameter ::  Zf = 4 !valleys to f-scatter into
    real (kind=dbl), parameter :: DOevcm=4.0e8   ! Fitting parameter optical scattering - order of magnitude 1e9-1e10 (eV/cm)
    real (kind=dbl),parameter :: DI = DOevcm*q*100 !used in intervalley (eV/m)
    real (kind=dbl),parameter :: hbaromegaf = 110E-3*q !Threshold for f-scattering LA phonons (TO has higher threshold)
    real (kind=dbl),parameter :: hbaromegag = 165E-3*q !Threshold for g-scattering (only LO allowed, see Nava "Electron Effective masses..."
    !real (kind=dbl),parameter :: DA = DAev*q
    
    
	real (kind=dbl),parameter :: emmax = 0.6476102      !maximum for emission Q^2 *Bose-Einstein function    
    real (kind=dbl),parameter ::  mdos =(ml*mt*mt)**(1/3.0) !DOS effective mass in kg
	real (kind=dbl),parameter ::  mavg = sqrt(ml*mt)/mdos   !NO longer used, remove
	real (kind=dbl),parameter :: C6    = mdos/(2*hbar**2*rho) !a useful parameter
	real (kind=dbl),dimension(4) :: g

    
    real (kind=dbl), dimension(6,3,3) :: Tm, Tminv  !Heering-Voigt matrix and its inverse (ISSUE: Tm does not seem to be used, Tminv is used in intervalley scattering but the definition has disappeared in this verson)
	integer, dimension(6,3,3) :: Pm  !six 3x3 permutation (rotation) matrices going from (x,y,z) to "(t,t,l)" 
    integer :: i,vali,errors,bin,valdum !various looping variables
    integer, dimension(8) :: scatstat !collects scattering statistics (for diagnostics only)
    real (kind=dbl),parameter :: mint = -30.0e-9      ! first travel time
    real (kind=dbl),parameter :: maxt = 100.0e-9      ! maximum travel time
    real (kind=dbl),parameter :: deltat =(maxt-mint)/bins; !width of the bins
        
	real (kind=dbl):: kT !kB*TL
    real (kind=dbl)  :: tmp,t,maxXi2,C1,C3,DA,tprev
    real (kind=dbl),dimension(3) :: x,p,xprev,E

    real (kind=dbl):: maxen, gamma_guess, gammamax, gammamax_1, gammamax_2, toterrmargin, toterr, tau0, tau1, tau2, tc, kmax_low, kmax_mid
    real (kind=dbl):: En, errmargin
    integer,dimension(bins) :: N
	real (kind=dbl),dimension(bins) :: Nx,Nt,Nxt,Ntt
    real (kind=dbl),dimension(6) :: xinval,tinval,Et,xint,xv
	
    Tm = 0.0   !unused ??
	Tminv = 0.0  ! Tminv definition missing
	Pm = 0
    do vali = 1,6
    	do i = 1,3
        	Pm(vali,i,1+mod(i-1+(vali+1)/2,3)  )  = 1 !Axis Permutation Matrices, vali=1,2:120deg rotation around [1,1,1], vali=3,4: -120deg, vali=5,6: identity matrix
		end do
	end do
	
	!can be used kept for diagnostics, Warning: gives lots of output
	!print *,'indata',L,TL,Ez,B,Xid/q,Xiu/q,intervalley_scattering,autogammacorrect,val,xtotinval,ttotinval,Ettot,jz2
	
    kT = kB*TL
	DA=max(abs(Xid),abs(Xid+Xiu)) ! DA now only used in C1 and for initial gamma guess
	maxXi2=max(Xid**2,(Xid+Xiu)**2)
	g(1)=kT**2/(hbar**2*vl**2) *emmax*maxXi2*maxg(1)  !SUSPICIOUS factor 2 REMOVED
	g(2)=kT**2/(hbar**2*vt**2) *emmax*Xiu**2*maxg(2)  !SUSPICIOUS factor 2 REMOVED
    C1 = mdos**1.5*sqrt(2.0)*DA**2*kT /(pi*hbar**4*rho*((vl+2*vt)/3)**2)   !(J^-0.5 s^-1) !C1 is the constant in Jacoboni 2.4.19, now only used for initial gamma guess
    !C2 = mdos**0.5*DA**2*kT**3 /(2**2.5*pi*hbar**4*rho*vavg**4) !(J^+0.5 s^-1) !no longer used
    C3 = mdos**1.5*DI**2/(sqrt(2.0)*pi*hbar**2*rho) !intervalley Jacoboni (3.74) (J^-0.5 s^-1)
	!C5 = mdos*DA**2/(rho*vavg*hbar**2) !(m^2 s^-1) !no longer used
         
    scatstat = 0
    E(1) = 0
    E(2) = 0
    E(3) = Ez
    call random_seed
    ! autocorrection of gamma0, the scattering rate ceiling ========
    !if (gmax == 0) then
    maxen = 0.5*(10*kB*TL + q*Ez/3.0e4) !Estimated cutoff energy in eV
    gamma_guess = 3*C1*sqrt(maxen)
    gammamax=gamma_guess
    kmax_low = 0.0; kmax_mid = 0.0
    toterrmargin=0.0; toterr = 0.0
    if (autogammacorrect > 0) then
      	    x = 0.0; p = 0.0
            tau0 = 1/gammamax  !Jacoboni ch 3.4
            valdum = val; En=0.0; errmargin=0.0 !errmargin is the ratio between calculated scattering rate and the scattering rate ceiling
            do i = 1, autogammacorrect
              	call Scatter6(TL,E,B,t,x,p,valdum,En,errmargin)
              	toterrmargin = toterrmargin + errmargin**6  !a softmax function for the error
            end do
            toterr = toterrmargin**(1/6.0) !a softmax function for the error
            gammamax = gammamax * sqrt(toterr)*1.2 !redefine gamma0 with extra safety factor 20%
    end if
    gmax = gammamax
    contains
    include 'scatter6sub.f90'
end subroutine gammamaxfind

subroutine efieldfind(z_bin,xy_bin,L_over_W,N,E) BIND(C,NAME="efieldfind_")
    implicit none   !All variables need to be properly declared E,N,M,
    integer,INTENT(IN) :: z_bin,xy_bin,L_over_W
    integer :: t,z,x,y,i,j,k
    real (kind=8) :: dist, ratio
    real (kind=8),INTENT(INOUT) :: E(100,29,29,3)
    real (kind=8),INTENT(IN), dimension(100,29,29) :: N
    !real (kind=8),INTENT(IN), dimension(-3:3,-3:3,-3:3,3) :: M
    !real (kind=4),INTENT(IN), dimension(z_bin+n_n_s*2,xy_bin+n_n_s*2,xy_bin+n_n_s*2,t_bin) :: N
    !real (kind=4),INTENT(IN), dimension(n_n_s*2+1,n_n_s*2+1,n_n_s*2+1) :: M
    
     ratio = z_bin/xy_bin/L_over_W
     do z = 1, z_bin
        do x = 1, xy_bin
            do y = 1, xy_bin
                if (N(z,x,y) .NE. 0.0) then
                    dist = E(z,x,y,3)
                    do i = 1, z_bin
                        do j = 1, xy_bin
                            do k = 1, xy_bin
                               ! dist = sqrt(real( (z-i)**2 + (x-j)**2 +(y-k)**2 )) # testaatt utgå ifrån en punkt till alla andra punkter ha ett test för att ta bort on den är noll 
                                E(i,j,k,3) = E(i,j,k,3) - ratio*N(z,x,y)*(y-k)/( (z-i)**2 + (ratio*(x-j))**2 +(ratio*(y-k))**2 )**1.5 !/sqrt(real( (z-i)**2 + (x-j)**2 +(y-k)**2 ))
                            end do                         
                        end do                                    !print *, E(z,x,y,t,3)
                    end do
                    E(z,x,y,3) = dist
                endif
            end do
        end do
     end do
     do z = 1, z_bin
        do x = 1, xy_bin
            do y = 1, xy_bin
                if (N(z,x,y) .NE. 0.0) then
                    dist = E(z,x,y,2)
                    do i = 1, z_bin
                        do j = 1, xy_bin
                            do k = 1, xy_bin
                               ! dist = sqrt(real( (z-i)**2 + (x-j)**2 +(y-k)**2 )) # testaatt utgå ifrån en punkt till alla andra punkter ha ett test för att ta bort on den är noll 
                                E(i,j,k,2) = E(i,j,k,2) - ratio*N(z,x,y)*(x-j)/( (z-i)**2 + (ratio*(x-j))**2 +(ratio*(y-k))**2 )**1.5 !/sqrt(real( (z-i)**2 + (x-j)**2 +(y-k)**2 ))
                            end do                         
                        end do                                    !print *, E(z,x,y,t,3)
                    end do
                    E(z,x,y,2) = dist
                endif
            end do
        end do
     end do
     do z = 1, z_bin
        do x = 1, xy_bin
            do y = 1, xy_bin
                if (N(z,x,y) .NE. 0.0) then
                    dist = E(z,x,y,1)
                    do i = 1, z_bin
                        do j = 1, xy_bin
                            do k = 1, xy_bin
                               ! dist = sqrt(real( (z-i)**2 + (x-j)**2 +(y-k)**2 )) # testaatt utgå ifrån en punkt till alla andra punkter ha ett test för att ta bort on den är noll 
                                E(i,j,k,1) = E(i,j,k,1) - N(z,x,y)*(z-i)/( (z-i)**2 + (ratio*(x-j))**2 +(ratio*(y-k))**2 )**1.5 !/sqrt(real( (z-i)**2 + (x-j)**2 +(y-k)**2 ))
                            end do                         
                        end do                                    !print *, E(z,x,y,t,3)
                    end do
                    E(z,x,y,1) = dist
                endif
            end do
        end do
     end do
    
    !do t =1, t_bin
    !    do z = 1, z_bin
    !        do x = 1, xy_bin
    !            do y = 1, xy_bin
    !                do i = -n_n_s,n_n_s
    !                    do j = -n_n_s,n_n_s
    !                        do k = -n_n_s,n_n_s
    !                            E(t,z,x,y,1) = E(t,z,x,y,1) - N(t,z+n_n_s+i,x+n_n_s+j,y+n_n_s+k)*M(i,j,k,1)
    !                            E(t,z,x,y,2) = E(t,z,x,y,2) - N(t,z+n_n_s+i,x+n_n_s+j,y+n_n_s+k)*M(i,j,k,2)
    !                            E(t,z,x,y,3) = E(t,z,x,y,3) - N(t,z+n_n_s+i,x+n_n_s+j,y+n_n_s+k)*M(i,j,k,3)
    !                            !print *, E(z,x,y,t,3)
    !                        end do 
    !                    end do    
    !                end do
    !            end do
    !        end do
    !    end do
    !end do
end subroutine efieldfind

subroutine electronsmooth(Pos_size,Pos,N,NN) BIND(C,NAME="electronsmooth_")
    implicit none   !All variables need to be properly declared E,N,M,
    integer,INTENT(IN) ::Pos_size
    integer :: t,i,j,k,f
    integer ,INTENT(IN) :: Pos(100*100*29*29,4)
    real (kind=8) ,INTENT(IN) :: N(100*100*29*29)
    real (kind=8),INTENT(INOUT), dimension(100,100,29,29) :: NN !dimension(100,104,-16:16,-16:16)

    
    do t =1, Pos_size
        NN(Pos(t,4)+1,Pos(t,1)+1,Pos(t,2)+14+1,Pos(t,3)+14+1) = NN(Pos(t,4)+1,Pos(t,1)+1,Pos(t,2)+14+1,Pos(t,3)+14+1) + N(t) ! fungerar inte längre NN är fel
        !print *, NN(Pos(t,4)+1,Pos(t,1)+1+n_n_s,Pos(t,2)+14+1+n_n_s,Pos(t,3)+14+1+n_n_s)
        !print *, t
        !print *, Pos(t,2)+14+1+n_n_s
        !print *, Pos(t,3)+14+1+n_n_s
        !NN(Pos(t,4)+1,Pos(t,1)+n_n_s,Pos(t,2)+14+1+n_n_s,Pos(t,3)+14+1+n_n_s) = N(t)*0.3
        !NN(Pos(t,4)+1,Pos(t,1)+2+n_n_s,Pos(t,2)+14+1+n_n_s,Pos(t,3)+14+1+n_n_s) = N(t)*0.3
        !NN(Pos(t,4)+1,Pos(t,1)+1+n_n_s,Pos(t,2)+14+n_n_s,Pos(t,3)+14+1+n_n_s) = N(t)*0.3
        !NN(Pos(t,4)+1,Pos(t,1)+1+n_n_s,Pos(t,2)+14+2+n_n_s,Pos(t,3)+14+1+n_n_s) = N(t)*0.3
        !NN(Pos(t,4)+1,Pos(t,1)+1+n_n_s,Pos(t,2)+14+1+n_n_s,Pos(t,3)+14+n_n_s) = N(t)*0.3
        !NN(Pos(t,4)+1,Pos(t,1)+1+n_n_s,Pos(t,2)+14+1+n_n_s,Pos(t,3)+14+2+n_n_s) = N(t)*0.3
        !
        !NN(Pos(t,4)+1,Pos(t,1)+n_n_s,Pos(t,2)+14+n_n_s,Pos(t,3)+14+1+n_n_s) = N(t)*0.2
        !NN(Pos(t,4)+1,Pos(t,1)+2+n_n_s,Pos(t,2)+14+n_n_s,Pos(t,3)+14+1+n_n_s) = N(t)*0.2
        !NN(Pos(t,4)+1,Pos(t,1)+n_n_s,Pos(t,2)+14+2+n_n_s,Pos(t,3)+14+1+n_n_s) = N(t)*0.2
        !NN(Pos(t,4)+1,Pos(t,1)+2+n_n_s,Pos(t,2)+14+2+n_n_s,Pos(t,3)+14+1+n_n_s) = N(t)*0.2
        !NN(Pos(t,4)+1,Pos(t,1)+n_n_s,Pos(t,2)+14+1+n_n_s,Pos(t,3)+14+n_n_s) = N(t)*0.2
        !NN(Pos(t,4)+1,Pos(t,1)+2+n_n_s,Pos(t,2)+14+1+n_n_s,Pos(t,3)+14+n_n_s) = N(t)*0.2
        !NN(Pos(t,4)+1,Pos(t,1)+n_n_s,Pos(t,2)+14+1+n_n_s,Pos(t,3)+14+2+n_n_s) = N(t)*0.2
        !NN(Pos(t,4)+1,Pos(t,1)+2+n_n_s,Pos(t,2)+14+1+n_n_s,Pos(t,3)+14+2+n_n_s) = N(t)*0.2
        !NN(Pos(t,4)+1,Pos(t,1)+1+n_n_s,Pos(t,2)+14+n_n_s,Pos(t,3)+14+2+n_n_s) = N(t)*0.2
        !NN(Pos(t,4)+1,Pos(t,1)+1+n_n_s,Pos(t,2)+14+n_n_s,Pos(t,3)+14+n_n_s) = N(t)*0.2
        !NN(Pos(t,4)+1,Pos(t,1)+1+n_n_s,Pos(t,2)+14+2+n_n_s,Pos(t,3)+14+n_n_s) = N(t)*0.2
        !NN(Pos(t,4)+1,Pos(t,1)+1+n_n_s,Pos(t,2)+14+2+n_n_s,Pos(t,3)+14+2+n_n_s) = N(t)*0.2
        !
        !NN(Pos(t,4)+1,Pos(t,1)+n_n_s,Pos(t,2)+14+n_n_s,Pos(t,3)+14+2+n_n_s) = N(t)*0.1
        !NN(Pos(t,4)+1,Pos(t,1)+n_n_s,Pos(t,2)+14+n_n_s,Pos(t,3)+14+n_n_s) = N(t)*0.1
        !NN(Pos(t,4)+1,Pos(t,1)+n_n_s,Pos(t,2)+14+2+n_n_s,Pos(t,3)+14+n_n_s) = N(t)*0.1
        !NN(Pos(t,4)+1,Pos(t,1)+n_n_s,Pos(t,2)+14+2+n_n_s,Pos(t,3)+14+2+n_n_s) = N(t)*0.1
        !NN(Pos(t,4)+1,Pos(t,1)+2+n_n_s,Pos(t,2)+14+n_n_s,Pos(t,3)+14+2+n_n_s) = N(t)*0.1
        !NN(Pos(t,4)+1,Pos(t,1)+2+n_n_s,Pos(t,2)+14+n_n_s,Pos(t,3)+14+n_n_s) = N(t)*0.1
        !NN(Pos(t,4)+1,Pos(t,1)+2+n_n_s,Pos(t,2)+14+2+n_n_s,Pos(t,3)+14+n_n_s) = N(t)*0.1
        !NN(Pos(t,4)+1,Pos(t,1)+2+n_n_s,Pos(t,2)+14+2+n_n_s,Pos(t,3)+14+2+n_n_s) = N(t)*0.1
        
        !print *, N(t)*M(0,0,1)
        !do i = -n_n_s, n_n_s
        !    do j = -n_n_s, n_n_s
        !         do k = -n_n_s, n_n_s
        !            f = M(i,j,k)
        !            NN(Pos(t,4)+1,Pos(t,1)+1+n_n_s+i,Pos(t,2)+14+1+n_n_s+j,Pos(t,3)+14+1+n_n_s+k) = N(t)*f
        !
        !         end do 
        !    end do    
        !end do
    end do
end subroutine electronsmooth