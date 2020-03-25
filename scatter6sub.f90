subroutine Scatter6(TL,E3D1,B,t,x,p,val,En,errmargin) !out: x,p,val,En,errmargin  
    !  Evolves the electron state by the E- and B- fields
    !  and by scattering
    implicit none
    real (kind=dbl),intent(in) :: TL,B  ! Lattice temperature, Efield in z-direction, magnetic field in y-direction
    real (kind=dbl),intent(inout) :: t    ! time 
	integer, intent(inout) :: val			! valley
    real (kind=dbl),dimension(3),intent(inout) :: x,p,E3D1 !(x,y,z) position and (px,py,pz) momentum
    real (kind=dbl),dimension(3) :: E3D
    real (kind=dbl),intent(out) :: En,errmargin		! energy and errormargin
    real (kind=dbl) :: qE, qB, cc, sqm, sqEn, r2, klen, beta, agm, tr, ts
    real (kind=dbl) :: Qph, qth, EQs, ksprim, Eprim, cosa, sina, cosb, sinb, phi, cosphi, sinphi, xsp, ysp, zsp
    real (kind=dbl),dimension(3) :: k, ks, qs, randsphere, xd, pd, tmp3
    real (kind=dbl),dimension(2) :: ksproj, scd, tmp2, tmp_2
    real (kind=dbl),dimension(8) :: gam
    integer :: sgn, scatinit
    qB = q*B   !defined for compactness
    E3D = E3D1*q ![0,0,33333]*q   !E3D*q
    xd = 0.0    !xd is the change in position by drift in the field
    pd = 0.0
	scatinit =sum(scatstat(1:4))
	gam = 0.0
	gam(1) = C6/vl * g(1)
	gam(2) = C6/vt * g(2)
            
						
	do 
		if (sum(scatstat(1:4)) > scatinit) exit !only perform one real scattering event
        		
		!update time
		call random_number(tmp_2)
		!tc = -tau0*log(tmp) !Exponentially distributed random timestep, Jacoboni (3.4.5)
        if (abs(B)<1e-6) then	
                
            if (sqrt(sum(p*p))/hbar < kmax_low ) then !starting with low energy
		        tc = -tau2*log(tmp_2(1))
                xd(:) = p(:)/m(val,:)*tc+ E3D/(2*m(val,:))*tc**2   ! moved distance= velocity* time (in x&y directions)+ extra term for the z-direction = due to acceleration in the E-field along z
        		pd(:) =  E3D*tc ! z-momentum increases also due to Ez field (???, why not just write: p(3)= p(3) + qE*tc, slightly faster ???)                    
                k = (p+pd)/hbar                  !k wave vector
                klen = sqrt(sum(k*k))		!length of wave vector k (same as length of ks)                   
                if (klen > kmax_mid) then 
                    tr = -sum(p*E3D)/sum(E3D*E3D) + sqrt( (sum(p*E3D)/sum(E3D*E3D))**2 +(kmax_low**2*hbar**2-sum(p*p))/sum(E3D*E3D) )
                    ts = -sum(p*E3D)/sum(E3D*E3D) + sqrt( (sum(p*E3D)/sum(E3D*E3D))**2 +(kmax_mid**2*hbar**2-sum(p*p))/sum(E3D*E3D) )
                    tc =-tau0*log(tmp_2(1)) + ts - tr*tau0/tau2-(ts-tr)*tau0/tau1
                    xd(:) = p(:)/m(val,:)*tc+ E3D/(2*m(val,:))*tc**2   ! moved distance= velocity* time (in x&y directions)+ extra term for the z-direction = due to acceleration in the E-field along z
        		    pd(:) =  E3D*tc ! z-momentum increases also due to Ez field (???, why not just write: p(3)= p(3) + qE*tc, slightly faster ???)  
                    r2 = tmp_2(2)*gammamax
                else if (klen > kmax_low) then 
                    tr = -sum(p*E3D)/sum(E3D*E3D) + sqrt( (sum(p*E3D)/sum(E3D*E3D))**2 +(kmax_low**2*hbar**2-sum(p*p))/sum(E3D*E3D) )
                    tc =-tau1*log(tmp_2(1))+tr*(1-tau1/tau2)
                    xd(:) = p(:)/m(val,:)*tc+ E3D/(2*m(val,:))*tc**2   ! moved distance= velocity* time (in x&y directions)+ extra term for the z-direction = due to acceleration in the E-field along z
        		    pd(:) =  E3D*tc ! z-momentum increases also due to Ez field (???, why not just write: p(3)= p(3) + qE*tc, slightly faster ???)  
                    r2 = tmp_2(2)*gammamax_1
                else 
                    r2 = tmp_2(2)*gammamax_2
                end if
                    
            else if (sqrt(sum(p*p))/hbar < kmax_mid ) then !starting with midel energy
                tc = -tau1*log(tmp_2(1))
                xd(:) = p(:)/m(val,:)*tc+ E3D/(2*m(val,:))*tc**2   ! moved distance= velocity* time (in x&y directions)+ extra term for the z-direction = due to acceleration in the E-field along z
                pd(:) =  E3D*tc ! z-momentum increases also due to Ez field (???, why not just write: p(3)= p(3) + qE*tc, slightly faster ???)  
                k = (p+pd)/hbar                  !k wave vector
                klen = sqrt(sum(k*k))		!length of wave vector k (same as length of ks)                      
                if (klen > kmax_mid) then  
                    tr = -sum(p*E3D)/sum(E3D*E3D) + sqrt( (sum(p*E3D)/sum(E3D*E3D))**2 +(kmax_mid**2*hbar**2-sum(p*p))/sum(E3D*E3D) )
                    tc =-tau0*log(tmp_2(1))+tr*(1-tau0/tau1)
                    xd(:) = p(:)/m(val,:)*tc+ E3D/(2*m(val,:))*tc**2   ! moved distance= velocity* time (in x&y directions)+ extra term for the z-direction = due to acceleration in the E-field along z
                    pd(:) =  E3D*tc ! z-momentum increases also due to Ez field (???, why not just write: p(3)= p(3) + qE*tc, slightly faster ???)  
                    r2 = tmp_2(2)*gammamax
                else 
                    r2 = tmp_2(2)*gammamax_1
                end if 
                        
            else !starting with high energy
                tc = -tau0*log(tmp_2(1))  !tau most relate to the lowe tc
                xd(:) = p(:)/m(val,:)*tc+ E3D/(2*m(val,:))*tc**2   ! moved distance= velocity* time (in x&y directions)+ extra term for the z-direction = due to acceleration in the E-field along z
                pd(:) =  E3D*tc ! z-momentum increases also due to Ez field (???, why not just write: p(3)= p(3) + qE*tc, slightly faster ???)  
                r2 = tmp_2(2)*gammamax
		    end if
		else !in case there is magnetic field
   !         tc = -tau0*log(tmp_2(1))
			!cc = sqrt(m(val,3)/m(val,1))
   !     	sqm = sqrt(m(val,3)*m(val,1))
   !     	xd(1) = cc*(p(1)/qB + qE*m(val,1)/qB**2)*sin(qB/sqm*tc) - p(3)/qB*(1-cos(qB/sqm*tc)) - qE/qB*tc !as above but with more complicated expressions in the prescence of B-field
   !     	xd(2) = p(2)/m(val,2)*tc																		!must FIND these expressions in my notes somewhere
   !     	xd(3) = 1/cc*p(3)/qB*sin(qB/sqm*tc) + (qE*m(val,1)/qB**2 + p(1)/qB)*(1-cos(qB/sqm*tc))
   !     	pd(1) = -qB*xd(3)
   !         pd(3) = qE*tc + qB*xd(1)
            !call random_number(tmp)
            !r2 = tmp*gammamax
		end if
        

        
        En = 0.5*sum(p**2/m(val,:)) !kinetic energy 1/2*mv^2
        sqEn = sqrt(En)
		k = (p+pd)/hbar                  !k wave vector
		klen = sqrt(sum(k*k))		!length of wave vector k (same as length of ks)
		ks = matmul(Pm(val,:,:),k)	!k in (X,Y,Z) coordinate system (Z parallel to longitudinal direction)
		cosa = ks(3)/klen			!cosine of alpha (alpha = angle between k and Z-axis, [0,pi] )
		sina = sqrt(1-cosa**2)		!sine of alpha  
		g(3) = g(1) + maxXi2*klen**2*maxg(3) !SUSPICIOUS factor 2 REMOVED
		g(4) = g(2) + (Xiu*klen)**2 *maxg(4) !SUSPICIOUS factor 2 REMOVED
		gam(3) = C6/vl * g(3)
		gam(4) = C6/vt * g(4)
			
		!print *, 'gamma1', gam(1), 'gamma2', gam(2) !for diagnostics
		!print * ,'En ',En,' klen ',klen             !for diagnostics
                
        if (intervalley_scattering == 1) then   ! ****  This needs to be updated in .v6 
            if (En > hbaromegaf) then
                gam(3) = Zf*C3/hbaromegaf*(1.0+1.0/(exp(hbaromegaf/(kB*TL))-1))*sqrt(En-hbaromegaf) 	!f-scattering by emission
			end if
            gam(4) = Zf*C3/hbaromegaf/(exp(hbaromegaf/(kB*TL))-1)*sqrt(En+hbaromegaf)  				!f-scattering by absorbtion
			if (En > hbaromegag) then
				gam(5) = C3/hbaromegag*(1.0+1.0/(exp(hbaromegag/(kB*TL))-1))*sqrt(En-hbaromegag)  	!g-scattering by emission
			end if
            gam(6)=C3/hbaromegag/(exp(hbaromegag/(kB*TL))-1)*sqrt(En+hbaromegag)     !g-scattering by absorbtion
		end if !(intervalley_scattering =1)

                
        ! detta är dumt borde flyttas upp till det som är ovanför
		!call random_number(tmp)
!            if (klen < kmax_low) then !for two-level self-scattering method, Jacoboni ch 3.4
!                if (sqrt(sum(p*p))/hbar < kmax_low ) then !starting with low energy
!                    r2 = tmp*gammamax_2
!                    errmargin = sum(gam)/gammamax_2 !does the total scattering rate exceed the maximukm scattering rate
!                else if (sqrt(sum(p*p))/hbar < kmax_mid ) then
!                    r2 = tmp*gammamax_1
!                    errmargin = sum(gam)/gammamax_1 !does the total scattering rate exceed the maximukm scattering rate 
!                else 
!                    r2 = tmp*gammamax
!                    errmargin = sum(gam)/gammamax !does the total scattering rate exceed the maximukm scattering rate
!                end if
!            else if (klen < kmax_mid) then 
!                if (sqrt(sum(p*p))/hbar < kmax_mid ) then
!                    r2 = tmp*gammamax_1
!                    errmargin = sum(gam)/gammamax_1 !does the total scattering rate exceed the maximukm scattering rate
!                else
!	                r2 = tmp*gammamax
!                    errmargin = sum(gam)/(gammamax) !does the total scattering rate exceed the maximukm scattering rate
!                end if
!            else
!	            r2 = tmp*gammamax
!                errmargin = sum(gam)/(gammamax) !does the total scattering rate exceed the maximukm scattering rate
!            end if
                
        errmargin = sum(gam)/(gammamax)
        t = t + tc !update time
        x = x + xd !update x
		p = p + pd !update momentum
				
		acoustic: if (r2 < sum(gam(1:4))) then           !acoustic phonon scattering
        
            call random_number(tmp3) !three (0-1) random numbers
			beta = pi * tmp3(1)      !beta = angle between phonon wavevector q and Z-axis   **********   CHECK this does not give isotropic q distribution ???? , should be "cosb=2*rand-1" 
			cosb = cos(beta)
			sinb = sin(beta)
			phi = 2*pi*tmp3(2)		  !phi azimuthal angle between ks and q projections in X-Y plane
			cosphi = cos(phi)
			sinphi = sin(phi)
			agm = mdos/ml*cosb**2 + mdos/mt*sinb**2  !useful quantity below
				
			sgn = 0 !(sgn indicate type of scattering event: -1 phonon emission, +1 phonon absorption)
			if (r2 < sum(gam(1:2))) then   !absorption
				if (r2 < gam(1)) then      !longitudinal phonon
					Qph = 2.0*( -klen*(mdos/ml*cosa*cosb + mdos/mt*sina*sinb*cosphi) + mdos*vl/hbar )/agm  !length of q, Notes w/o H-W Max 1 (Q1)
					if (Qph > 0) then
						if  (tmp3(3)*g(1) < sinb/agm*Qph**2*(Xid+Xiu*cosb**2)**2 /(exp(hbar*vl*Qph/(kB*TL))-1.0) ) then !SUSPICIOUS factor 2 REMOVED
						scatstat(1) = scatstat(1) + 1   !absorption of longitudinal acoustic phonon
						sgn = 1
						end if
					end if
				else						!transveral phonon
					Qph = 2.0*( -klen*(mdos/ml*cosa*cosb + mdos/mt*sina*sinb*cosphi) + mdos*vt/hbar )/agm  !length of q, Notes w/o H-W Max 1 (Q1)
					if (Qph > 0) then
						if (tmp3(3)*g(2) < sinb/agm*Qph**2*(Xiu*sinb*cosb)**2 /(exp(hbar*vt*Qph/(kB*TL))-1.0) ) then   !SUSPICIOUS factor 2 REMOVED
						scatstat(2) = scatstat(2) + 1   !absorption of transversal acoustic phonon
						sgn = 1
						end if
					end if
				end if                      
			else  !emission
				if (r2 < sum(gam(1:3))) then      !longitudinal phonon
					Qph = 2.0*( klen*(mdos/ml*cosa*cosb + mdos/mt*sina*sinb*cosphi) - mdos*vl/hbar )/agm  !length of q, Notes w/o H-W Max 2 (Q2)
					if (Qph > 0) then
						if  (tmp3(3)*g(3) < sinb/agm*Qph**2*(Xid+Xiu*cosb**2)**2 *(1.0+1.0/(exp(hbar*vl*Qph/(kB*TL))-1.0) ) )then  !SUSPICIOUS factor 2 REMOVED
						scatstat(3) = scatstat(3) + 1
						sgn = -1
						end if
					end if
				else						!transveral phonon
					Qph = 2.0*( klen*(mdos/ml*cosa*cosb + mdos/mt*sina*sinb*cosphi) - mdos*vt/hbar )/agm  !length of q, Notes w/o H-W Max 2 (Q2)
					if (Qph > 0) then
						if (tmp3(3)*g(4) < sinb/agm*Qph**2*(Xiu*sinb*cosb)**2 *(1.0+1.0/(exp(hbar*vt*Qph/(kB*TL))-1.0) ) ) then   !SUSPICIOUS factor 2 REMOVED
						scatstat(4) = scatstat(4) + 1
						sgn = -1
						end if
					end if
				end if  
			end if
					
			if (sgn /= 0) then                          !if there was a real scattring event
				ksproj = (/ks(1), ks(2)+1e-50_dbl/)     ! ks projected on X-Y plane ks(2)+eps ?)
				scd = ksproj/sqrt(sum(ksproj*ksproj))   ! (1) cos delta (2) sin delta   (delta=angle between ks projection on X-Y plane and X axis)
				qs = Qph * (/sinb*(scd(1)*cosphi-scd(2)*sinphi), sinb*(scd(2)*cosphi+scd(1)*sinphi), cosb/)  !phonon wavevector in XYZ coordinates "(sin(beta)*cos(delta+phi),sin(beta)*sin(delta+phi),cos(beta))"
				p = hbar * (k + sgn*matmul(qs, Pm(val,:,:))) !updating momentum after scattering by adding or subtracting phonon wave vector. The transformation works since Pm are orthogonal matrices (transpose=inverse)
			end if
				                  
        else  acoustic !not acoustic scattering   		! ***** This case needs to be fixed to work in .v6  ****
            f_or_g: if (r2 < sum(gam(1:4))) then    		!  phonon f-scattering
         		call random_number(tmp2); !call random_number(tmp2);
                val = 1 + nint(tmp2(1)) + 2*mod((val-1)/2 + 1 + nint(tmp2(2)), 3)  ! function for random f-scattering
                EQs = hbaromegaf
                if (r2 < sum(gam(1:3))) then     	 ! f-scattering through phonon emission
                    scatstat(3)= scatstat(3) + 1
                    Eprim = En-EQs
                else		                         ! f-scattering through phonon absorption
                    scatstat(4)= scatstat(4) + 1
                    Eprim = En+EQs
				end if
                
                ksprim = sqrt(2*mdos*Eprim)/hbar  
                call random_number(tmp2)           
                phi = 2*pi*tmp2(1)							  ! generating isotropically distributed unit vector
                !call random_number(tmp)
                zsp = 2*tmp2(2)-1   						  !random z value between -1 and 1
                xsp = sqrt(1-zsp**2)*cos(phi) 
                ysp = sqrt(1-zsp**2)*sin(phi)
                randsphere = (/xsp,ysp,zsp/) 				   !random unit vector
                p=hbar*ksprim*matmul(randsphere,Tminv(val,:,:))
           
            else f_or_g           !  phonon g-scattering
                if (r2 < sum(gam)) then    
                    val = val + 2*mod(val,2) - 1     ! function for g-scattering (jump to opposite valley)
                    EQs = hbaromegag
                    if (r2 < sum(gam(1:5))) then     ! g-scattering through phonon emission
                        scatstat(5) = scatstat(5) + 1
                        Eprim = En-EQs
                    else                      !g-scattering through phonon absorption
                        scatstat(6) = scatstat(6) + 1
                        Eprim = En+EQs
					end if
                      
                	ksprim = sqrt(2*mdos*Eprim)/hbar  
                	call random_number(tmp2)           
                	phi = 2*pi*tmp2(1)							  ! generating isotropically distributed unit vector
                	!call random_number(tmp)
                	zsp = 2*tmp2(2)-1   						  !random z value between -1 and 1
                	xsp = sqrt(1-zsp**2)*cos(phi) 
                	ysp = sqrt(1-zsp**2)*sin(phi)
                	randsphere = (/xsp,ysp,zsp/) 				   !random unit vector
                	p=hbar*ksprim*matmul(randsphere,Tminv(val,:,:))
				end if
			end if f_or_g
		end if acoustic
		!exit ! to just go through once even if it was just a self-scattering event
	end do	
end subroutine Scatter6  !Scatter ends here  ===================================================================================   