
!------------------------------------------------------------------
	subroutine getheating0(teqx,teqy,teqs,ionemittx,ionemitty,iondponp,ionsigs,ionnum)
	parameter(pi=3.141592,clight=2.998e8,qe=1.602e-19)
	parameter(aatomi=1.,zatomi=1.)
	common/parms/pnum,zatom,aatom,epsx,epsy,gamma,sigs,dponp,coulomb
	common/betafunc/ betax0,betay0,scool,circ,betaxI0,betayI0,dispex,dispix
	real ionemittx,ionemitty,iondponp,ionsigs,ionnum

    real(8):: me,mi,kb
    	me = 9.10938356e-31
    	mi = 1836.0*me*zatomi
    	kb = 1.38064852e-23
!average density of ion beam (estimate)
    	denion = ionnum/(1.0*ionsigs)/(1.0*sqrt(ionemittx*betaxI0))/(1.0*sqrt(ionemitty*betayI0))/2/pi/sqrt(2*pi)
!fractor of cooling length on circumference
    	frcool = scool/circ
!get the tempreture of ion and electron beam
    	call temrest(temex,temey,temes,epsx,epsy,dponp,betaxI0,betay0,me,gamma)
    	call temrest(temIx,temIy,temIs,ionemittx,ionemitty,iondponp,betaxI0,betay0,mi,gamma)
	!write(*,*) "tem_e",temex,temey,temes
	!write(*,*) "tem_i",temIx,temIy,temIs

!calculation in rest frame (spitzer formula)
    	teq1 = 3.*me*mi/qe**2
    	teq2 = (4.*pi*8.854e-12)**2/(8.*sqrt(2*pi)*denion/gamma*zatom**2*zatomi**2*qe**2*coulomb)

!the results below are transfered to the lab frame by gamma (dt = gamma*dt')
!there teqx=d(emittx)/emittx/dt and teqs=d(sigmap**2)/sigmap**2/dt
		teqx = (temIx/temex-1.)/(teq1*teq2*(temex/me+temIx/mi)**(3./2.))/gamma
		teqy = (temIy/temey-1.)/(teq1*teq2*(temey/me+temIy/mi)**(3./2.))/gamma
		teqs = (temIs/temes-1.)/(teq1*teq2*(temes/me+temIs/mi)**(3./2.))/gamma
		!write(*,*) temex,temey,temes

	return
	end

!new subroutine
!get the beam tempreture in rest frame
    
    
	subroutine temrest(temx,temy,tems,ex,ey,dpp,betax,betay,mass,gamma)
	real(8):: mass
		clight=2.998e8
		beta=sqrt(1-1/gamma**2)
		sigmavx = sqrt(ex/betax)*beta*clight*gamma
		sigmavy = sqrt(ey/betay)*beta*clight*gamma
		sigmadpp = dpp*beta*clight

		temx = mass*(sigmavx**2)
		temy = mass*(sigmavy**2)
		tems = mass*sigmadpp**2
	return
	end
!------------------------------------------------------------------

	
!------------------------------------------------------------------
	
! Beam heating based on the Boltzmann transport equation and coulomb collision operator Cab
	
	subroutine getheating11(alfxC,alfyC,alfpC,ionemittx,ionemitty,iondponp,ionsigs,ionnum,steps)
!alfxC,alfyC,alfpC are velocity square growth rates in per sec
!In this subroutine, the vx, vy ,vz of ion beam and electron beam are not real velocity.They are just vx/beta/clight...
!Based on the formulae, the alfxC,alfyC and alfpC should should devided by (beta*clight)^3 in that frame
! only aplly uspan in the integration, should be more accurate
	parameter(pi=3.141592,clight=2.998e8,qe=1.602e-19,dcool=200.)
	parameter(aatomi=1.,zatomi=1.)
	common/parms/pnum,zatom,aatom,epsx,epsy,gamma,sigs,dponp,coulomb
	common/betafunc/ betax0,betay0,scool,circ,betaxI0,betayI0,dispex,dispix
	real ionemittx,ionemitty,iondponp,ionsigs,ionnum
	real Ix0,Ix1,Ix2,Iy0,Iy1,Iy2,Ip0,Ip1,Ip2,sumx,sumy,sump
	real iv2x,iv2y,iv2p,ev2x,ev2y,ev2p
	real sumx1,sumx2,sumx3,kvalue,kvalue1,kvalue2,gab,me
	real alfxC,alfyC,alfpC
	integer steps

	beta=sqrt(1.-1./gamma**2)
	me = 9.10938356e-31
	epsilon = 8.854187e-12
	qi = qe * zatomi
	gab = (qi*qe/me)**2 * coulomb /4/pi/epsilon**2
	!gab = real(gab)
	!fractor of cooling length on circumference
    	frcool = scool/circ
	
	sigix = sqrt(ionemittx*betaxI0)
	sigiy = sqrt(ionemitty*betayI0)
	sigis = ionsigs*gamma
	
! <vx2> <vy2> and <vp2> in rest frame
	iv2x = ionemittx / betaxI0 * gamma**2  
	iv2y = ionemitty / betayI0 * gamma**2 
	iv2p = iondponp**2 
	iv2max = iv2x
	if(iv2y.gt.iv2max) iv2max = iv2y
	if(iv2p.gt.iv2max) iv2max = iv2p
	
	ev2x = epsx / betax0 * gamma**2	 
	ev2y = epsy / betay0 * gamma**2	
	ev2p = dponp**2	
	ev2max = ev2x
	if(ev2y.gt.ev2max) ev2max = ev2y
	if(ev2p.gt.ev2max) ev2max = ev2p
	uspan = 5.*sqrt(iv2max+ev2max)
	sumx=0.0
	sumy=0.0
	sump=0.0	

	do i=1,steps
		do j=1,steps
			do k=1,steps
				
				ux = (2.*real(i+0.3)/real(steps) -1.) * uspan
				uy = (2.*real(j+0.3)/real(steps) -1.) * uspan
				up = (2.*real(k+0.3)/real(steps) -1.) * uspan
				uu = sqrt(ux**2 + uy**2 + up**2)
				
				call getI(Ix0,Ix1,Ix2,0.5/ev2x,0.5/iv2x,ux)
				call getI(Iy0,Iy1,Iy2,0.5/ev2y,0.5/iv2y,uy)
				call getI(Ip0,Ip1,Ip2,0.5/ev2p,0.5/iv2p,up)
				sumx1 = (uu**2-ux**2)/uu**3
				sumx2 = ux*uy/uu**3
				sumx3 = ux*up/uu**3
					
				sumx1 = sumx1 * Ix2 * Iy0 * Ip0/ev2x
				sumx2 = sumx2 *	Ix1 * Iy1 * Ip0/ev2y
				sumx3 = sumx3 * Ix1 * Iy0 * Ip1/ev2p
				sumx = sumx + sumx1 - sumx2 - sumx3
				
				sumy1 =	(uu**2-uy**2)/uu**3
				sumy2 = uy*ux/uu**3
				sumy3 = uy*up/uu**3
				sumy1 = sumy1 * Ix0 * Iy2 * Ip0/ev2y
				sumy2 = sumy2 *	Ix1 * Iy1 * Ip0/ev2x
				sumy3 = sumy3 * Ix0 * Iy1 * Ip1/ev2p
				sumy = sumy + sumy1 - sumy2 - sumy3
				
				sump1 =	(uu**2-up**2)/uu**3
				sump2 = up*ux/uu**3
				sump3 = up*uy/uu**3
				sump1 = sump1 * Ix0 * Iy0 * Ip2/ev2p
				sump2 = sump2 *	Ix1 * Iy0 * Ip1/ev2x
				sump3 = sump3 * Ix0 * Iy1 * Ip1/ev2y
				sump = sump + sump1 - sump2 - sump3
			   !write(*,*) ux,uy,up,sumx
			end do
		end do
	end do
	
	kvalue1 = gab * ionnum /(2.0*pi)/sqrt(2.0*pi)/sigis/sigix/sigiy/2.0/sqrt(2.0)
	kvalue2 = ((2.0*pi)**3*(sqrt(iv2x*iv2y)*sqrt(iv2p*ev2x)*sqrt(ev2y*ev2p)))

	kvalue =kvalue1/kvalue2
	sumindex = 8.0 * uspan * uspan * uspan / steps /steps /steps
	! taux=1/emitt*d(emitt)/dt = 2/sigmax*d(sigmax)/dt
	!<vx2>/2=sigmavx2	gaussian distribution
	alfxC = sumx * kvalue * sumindex/ev2x /gamma/(beta*clight)**3/2.
	alfyC = sumy * kvalue * sumindex/ev2y /gamma/(beta*clight)**3/2.
	alfpC = sump * kvalue * sumindex/ev2p /gamma/(beta*clight)**3/2.

	!write(*,*) sumx,sumy,sump
	return
	end
	

!------------------------------------------------------------------
	subroutine getI(I0,I1,I2,aa,bb,uuk)
	parameter(pi=3.141592)
	real I0,I1,I2,aa,bb,uuk,kk1,kk2

	kk1 = aa+bb
	kk2 = aa*bb
	I0 = sqrt(pi/kk1)*exp(-kk2/kk1*uuk*uuk)
	I1 = -bb/kk1*uuk*I0
	I2 = (0.5/kk1+bb*bb/kk1/kk1*uuk*uuk)*I0

	end

!------------------------------------------------------------------
!------------------------------------------------------------------	
	subroutine getheating12(alfxC,alfyC,alfpC,ionemittx,ionemitty,iondponp,ionsigs,ionnum,steps)
!alfx,alfy,alfp are emittance growth rates in per sec
!the vx, vy, vp are the real beam velocity, which are timed beta*clight already.
! aplly uspanx,uspany,uspanp in the integration, should be more accurate
	parameter(pi=3.141592,clight=2.998e8,qe=1.602e-19,dcool=200.)
	parameter(aatomi=1.,zatomi=1.)
	!parameter(aatomi=197.,zatomi=79.)	!gold beam
	common/parms/pnum,zatom,aatom,epsx,epsy,gamma,sigs,dponp,coulomb
	common/betafunc/ betax0,betay0,scool,circ,betaxI0,betayI0,dispex,dispix
	real ionemittx,ionemitty,iondponp,ionsigs,ionnum
	real Ix0,Ix1,Ix2,Iy0,Iy1,Iy2,Ip0,Ip1,Ip2,sumx,sumy,sump
	real iv2x,iv2y,iv2p,ev2x,ev2y,ev2p
	real sumx1,sumx2,sumx3,kvalue,gab,me
	integer steps
	beta=sqrt(1.-1./gamma**2)
	me = 9.10938356e-31
	epsilon = 8.854187e-12
	qi = qe * zatomi
	gab = (qi*qe/me)**2 * coulomb /4./pi/epsilon**2
	!gab = real(gab)
	!fractor of cooling length on circumference
    	frcool = scool/circ
	!write(*,*) scool,circ,frcool
	sigix = sqrt(ionemittx*betaxI0)
	sigiy = sqrt(ionemitty*betayI0)
	sigis = ionsigs*gamma

! <vx2> <vy2> and <vp2> in rest frame
	iv2x = ionemittx / betaxI0 * gamma**2  *(beta*clight)**2	*1.
	iv2y = ionemitty / betayI0 * gamma**2 *(beta*clight)**2 *1.
	iv2p = iondponp**2  *(beta*clight)**2*1.
	
	!iv2max = iv2x
	!if(iv2y.gt.iv2max) iv2max = iv2y
	!if(iv2p.gt.iv2max) iv2max = iv2p
	
	ev2x = epsx / betax0 * gamma**2	*(beta*clight)**2  *1.
	ev2y = epsy / betay0 * gamma**2	*(beta*clight)**2 *1.
	ev2p = dponp**2	*(beta*clight)**2*1.
	!ev2max = ev2x
	!if(ev2y.gt.ev2max) ev2max = ev2y
	!if(ev2p.gt.ev2max) ev2max = ev2p

	uspanx = 4.*sqrt(iv2x+ev2x)
	uspany = 4.*sqrt(iv2y+ev2y)
	uspanp = 4.*sqrt(iv2p+ev2p)
	!write(*,*) uspanx,uspany,uspanp

	sumx=0.0
	sumy=0.0
	sump=0.0
	
	do i=1,steps
		ux = (2.*real(i+0.13131)/real(steps) -1.) * uspanx
		do j=1,steps
			uy = (2.*real(j-0.23133)/real(steps) -1.) * uspany
			do k=1,steps	
				up = (2.*real(k+0.13133)/real(steps) -1.) * uspanp
				uu = sqrt(ux**2 + uy**2 + up**2)
				call getI(Ix0,Ix1,Ix2,0.5/ev2x,0.5/iv2x,ux)				
				call getI(Iy0,Iy1,Iy2,0.5/ev2y,0.5/iv2y,uy)
				call getI(Ip0,Ip1,Ip2,0.5/ev2p,0.5/iv2p,up)

				sumx1 = (uu**2-ux**2)/uu**3
				sumx2 = ux*uy/uu**3
				sumx3 = ux*up/uu**3
					
				sumx1 = sumx1 * Ix2 * Iy0 * Ip0/ev2x
				sumx2 = sumx2 *	Ix1 * Iy1 * Ip0/ev2y
				sumx3 = sumx3 * Ix1 * Iy0 * Ip1/ev2p
				sumx = sumx + sumx1 - sumx2 - sumx3
				!write(*,*) sumx,sumx1,sumx2,sumx3,sumx1 - sumx2 - sumx3
				
				sumy1 =	(uu**2-uy**2)/uu**3
				sumy2 = uy*ux/uu**3
				sumy3 = uy*up/uu**3
				sumy1 = sumy1 * Ix0 * Iy2 * Ip0/ev2y
				sumy2 = sumy2 *	Ix1 * Iy1 * Ip0/ev2x
				sumy3 = sumy3 * Ix0 * Iy1 * Ip1/ev2p
				sumy = sumy + sumy1 - sumy2 - sumy3
				
				sump1 =	(uu**2-up**2)/uu**3
				sump2 = up*ux/uu**3
				sump3 = up*uy/uu**3
				sump1 = sump1 * Ix0 * Iy0 * Ip2/ev2p
				sump2 = sump2 *	Ix1 * Iy0 * Ip1/ev2x
				sump3 = sump3 * Ix0 * Iy1 * Ip1/ev2y
				sump = sump + sump1 - sump2 - sump3
			    
			end do
		end do
	end do
	
	kvalue = gab * ionnum/(2.0*pi)/sqrt(2.0*pi)/sigis/sigix/sigiy/2.0/sqrt(2.0)
	kvalue = kvalue/(2.0*pi)**3/sqrt(iv2x*iv2y)/sqrt(iv2p*ev2x)/sqrt(ev2y*ev2p)
	
	!the number 8 is 2*2*2
	sumindex = 8.0 * uspanx * uspany * uspanp / steps /steps /steps
	!write(*,*) sumx1*sumindex,sumx2,sumx3,sumx
	! taux=1/emitt*d(emitt)/dt = 2/sigmax*d(sigmax)/dt
	!<vx2>/2=sigmavx2	gaussian distribution
	alfxC = sumx * kvalue * sumindex /gamma /ev2x/2.
	alfyC = sumy * kvalue * sumindex /gamma	/ev2y/2.
	alfpC = sump * kvalue * sumindex /gamma	/ev2p/2.
	
	return
	end



!------------------------------------------------------------------
	subroutine getIp(I0,I1,I2,aa,bb,Ma,Mb,uuk)
!new getI for condition with dispersion.
!special for p with Ma, Mb
	parameter(pi=3.141592)
	real I0,I1,I2,aa,bb,uuk,kk1,kk2,Ma,Mb

	kk1 = aa+bb
	kk2 = aa*bb
	I0 = sqrt(pi/kk1)*exp(-kk2/kk1*(uuk+Ma-Mb)**2)
	I1 = -bb/kk1*(uuk+Ma-Mb)*I0

	kk3 = kk1 - 2.*kk2*Ma*(uuk+Ma-Mb)
	kk3 = kk3 - 2*bb**2*(Mb-uuk)*(uuk+Ma-Mb)
	I2 = kk3*I0/2.0/kk1**2

	end

!------------------------------------------------------------------
!	function gauss_random(average,sigma)
!		real gauss_random1,gauss_random2
!		real average,sigma,mm,nn,w
!2		call random(xsl)
!		mm = 2.0*xsl-1.0
!		call random(xsl)
!		nn = 2.0*xsl-1.0
!		w = mm*mm + nn*nn
!		if(w.gt.1.0) goto 2
!		if(w.lt.3.0e-7) goto 2
!		gauss_random = mm*sqrt((-2.0*log(w))/w)*sigma+average
!	end

!------------------------------------------------------------------
!------------------------------------------------------------------		
	subroutine getheating13(alfxC,alfyC,alfpC,ionemittx,ionemitty,iondponp,ionsigs,ionnum,steps)
!alfxC,alfyC,alfpC are velocity square growth rates in per sec
!In this subroutine, the vx, vy ,vz of ion beam and electron beam are not real velocity.They are just vx/beta/clight...
!Based on the formulae, the alfxC,alfyC and alfpC should should devided by (beta*clight)^3 in that frame
! only aplly uspan in the integration, should be more accurate
! add dispersion (De,Di) 4/14/2021
	parameter(pi=3.141592,clight=2.998e8,qe=1.602e-19,dcool=200.)
	parameter(aatomi=1.,zatomi=1.)
	common/parms/pnum,zatom,aatom,epsx,epsy,gamma,sigs,dponp,coulomb
	common/betafunc/ betax0,betay0,scool,circ,betaxI0,betayI0,dispex,dispix
	real ionemittx,ionemitty,iondponp,ionsigs,ionnum
	real Ix0,Ix1,Ix2,Iy0,Iy1,Iy2,Ip0,Ip1,Ip2,sumx,sumy,sump
	real iv2x,iv2y,iv2p,ev2x,ev2y,ev2p,ev2p0,sumxf,sumyf,sumpf
	real sumx1,sumx2,sumx3,kvalue,kvalue1,kvalue2,gab,me
	real alfxC,alfyC,alfpC,Mii,Mei,Nx,pxi,xi
	integer steps

	beta=sqrt(1.-1./gamma**2)
	me = 9.10938356e-31
	epsilon = 8.854187e-12
	qi = qe * zatomi
	gab = (qi*qe/me)**2 * coulomb /4/pi/epsilon**2
	!gab = real(gab)
	!fractor of cooling length on circumference
    frcool = scool/circ
	
	sigix0 = sqrt(ionemittx*betaxI0)
	sigix = sqrt(ionemittx*betaxI0+(dispix*iondponp)**2)
	sigiy = sqrt(ionemitty*betayI0)
	sigis = ionsigs*gamma
	
! <vx2> <vy2> and <vp2> in rest frame
	iv2x = ionemittx / betaxI0 * gamma**2  
	iv2y = ionemitty / betayI0 * gamma**2 
	iv2p = iondponp**2*sigix0**2/sigix**2
	ivpk = dispix*iondponp**2/sigix**2

	ev2x = epsx / betax0 * gamma**2	 
	ev2y = epsy / betay0 * gamma**2	
	ev2p0 = dponp**2
	ev2p = dponp**2	* (epsx*betax0)/(epsx*betax0+(dispex*dponp)**2)
	evpk = dispex*dponp**2/(epsx*betax0+(dispex*dponp)**2)

	uspanx = 4.*sqrt(iv2x+ev2x)
	uspany = 4.*sqrt(iv2y+ev2y)
	uspanp = 4.*sqrt(iv2p+ev2p)

	sumxf=0.0
	sumyf=0.0
	sumpf=0.0
	Nx = 100.

	do l=1,Nx
		xi = gauss_random(0.0,sigix)

		Mii = ivpk*xi
		Mei = evpk*xi
		pxi = 1./sqrt(2.*pi)/sigix*exp(-xi**2/sigix**2/2.0)

		sumx=0.0
		sumy=0.0
		sump=0.0
		do i=1,steps
			do j=1,steps
				do k=1,steps
					
					ux = (2.*real(i+0.3)/real(steps) -1.) * uspanx
					uy = (2.*real(j+0.3)/real(steps) -1.) * uspany
					up = (2.*real(k+0.3)/real(steps) -1.) * uspanp
					uu = sqrt(ux**2 + uy**2 + up**2)
					
					call getI(Ix0,Ix1,Ix2,0.5/ev2x,0.5/iv2x,ux)
					call getI(Iy0,Iy1,Iy2,0.5/ev2y,0.5/iv2y,uy)
					call getIp(Ip0,Ip1,Ip2,0.5/ev2p,0.5/iv2p,Mei,Mii,up)

					sumx1 = (uu**2-ux**2)/uu**3
					sumx2 = ux*uy/uu**3
					sumx3 = ux*up/uu**3
						
					sumx1 = sumx1 * Ix2 * Iy0 * Ip0/ev2x
					sumx2 = sumx2 *	Ix1 * Iy1 * Ip0/ev2y
					sumx3 = sumx3 * Ix1 * Iy0 * Ip1/ev2p
					sumx = sumx + sumx1 - sumx2 - sumx3
					
					sumy1 =	(uu**2-uy**2)/uu**3
					sumy2 = uy*ux/uu**3
					sumy3 = uy*up/uu**3
					sumy1 = sumy1 * Ix0 * Iy2 * Ip0/ev2y
					sumy2 = sumy2 *	Ix1 * Iy1 * Ip0/ev2x
					sumy3 = sumy3 * Ix0 * Iy1 * Ip1/ev2p
					sumy = sumy + sumy1 - sumy2 - sumy3
					
					sump1 =	(uu**2-up**2)/uu**3
					sump2 = up*ux/uu**3
					sump3 = up*uy/uu**3
					sump1 = sump1 * Ix0 * Iy0 * Ip2/ev2p
					sump2 = sump2 *	Ix1 * Iy0 * Ip1/ev2x
					sump3 = sump3 * Ix0 * Iy1 * Ip1/ev2y
					sump = sump + sump1 - sump2 - sump3
				   !write(*,*) ux,uy,up,sumx
				end do
			end do
		end do
		sumxf = sumxf + sumx*pxi/Nx
		sumyf = sumyf + sumy*pxi/Nx
		sumpf = sumpf + sump*pxi/Nx
		!sumxf = sumxf + pxi/Nx
		!sumyf = sumyf + pxi/Nx
		!sumpf = sumpf + pxi/Nx
	end do

	kvalue1 = gab * ionnum /(2.0*pi)/sigis/sigiy/2.0
	kvalue2 = ((2.0*pi)**3*(sqrt(iv2x*iv2y)*sqrt(iv2p*ev2x)*sqrt(ev2y*ev2p)))

	kvalue =kvalue1/kvalue2
	sumindex = 8.0 * uspanx * uspany * uspanp / steps /steps /steps
	! taux=1/emitt*d(emitt)/dt = 2/sigmax*d(sigmax)/dt
	!<vx2>/2=sigmavx2	gaussian distribution
	alfxC = sumxf * kvalue * sumindex/ev2x /gamma/(beta*clight)**3/2.
	alfyC = sumyf * kvalue * sumindex/ev2y /gamma/(beta*clight)**3/2.
	alfpC = sumpf * kvalue * sumindex/ev2p0 /gamma/(beta*clight)**3/2.
	alfxC = alfxC + alfpC/betax0*dispex**2/epsx*dponp**2 

	return
	end
	



!------------------------------------------------------------------
!------------------------------------------------------------------	
	subroutine heatingIBS(alfxC,alfyC,alfpC,avebeta,steps)
!alfx,alfy,alfp are emittance growth rates in per sec
!the vx, vy, vp are the real beam velocity, which are timed beta*clight already.
! calculation on ion beam for benchmark with IBS heating
	parameter(pi=3.141592,clight=2.998e8,qe=1.602e-19,dcool=200.)
	parameter(aatomi=1.,zatomi=1.)
	common/parms/pnum,zatom,aatom,epsx,epsy,gamma,sigs,dponp,coulomb

	real Ix0,Ix1,Ix2,Iy0,Iy1,Iy2,Ip0,Ip1,Ip2,sumx,sumy,sump
	real ev2x,ev2y,ev2p
	real sumx1,sumx2,sumx3,sumy1,sumy2,sumy3,sump1,sump2,sump3
	real kvalue,gab,me,avebeta
	integer steps

	beta=sqrt(1.-1./gamma**2)
	me = 9.10938356e-31
	epsilon = 8.854187e-12
	
	gab = (qe*qe/me)**2 * coulomb /4/pi/epsilon**2
	!gab = real(gab)	
	sigex = sqrt(epsx*avebeta)
	sigey = sqrt(epsy*avebeta)
	siges = sigs*gamma
	 
	ev2x = epsx / avebeta * gamma**2 *(beta*clight)**2 *1.
	ev2y = epsy / avebeta * gamma**2 *(beta*clight)**2 *1.
	ev2p = dponp**2	*(beta*clight)**2*1.
	!ev2max = ev2x
	!if(ev2y.gt.ev2max) ev2max = ev2y
	!if(ev2p.gt.ev2max) ev2max = ev2p
	!uspan = 6.*sqrt(ev2max)

	uspanx = 4.*sqrt(ev2x+ev2x)
	uspany = 4.*sqrt(ev2y+ev2y)
	uspanp = 4.*sqrt(ev2p+ev2p)
	sumx=0.0
	sumy=0.0
	sump=0.0

	do i=1,steps
		do j=1,steps
			do k=1,steps
				
				ux = (2.*real(i+0.324)/real(steps) -1.) * uspanx
				uy = (2.*real(j+0.134)/real(steps) -1.) * uspany
				up = (2.*real(k+0.305)/real(steps) -1.) * uspanp

				uu = sqrt(ux**2 + uy**2 + up**2)

				call getI(Ix0,Ix1,Ix2,0.5/ev2x,0.5/ev2x,ux)				
				call getI(Iy0,Iy1,Iy2,0.5/ev2y,0.5/ev2y,uy)
				call getI(Ip0,Ip1,Ip2,0.5/ev2p,0.5/ev2p,up)
				!write(*,*) i,j,sumx
				
				sumx1 = (uu**2-ux**2)/uu**3*ux
				sumx2 = ux*uy/uu**3*uy
				sumx3 = ux*up/uu**3*up
					
				sumx1 = sumx1 * Ix1 * Iy0 * Ip0 / ev2x
				sumx2 = sumx2 *	Ix1 * Iy0 * Ip0 / ev2y
				sumx3 = sumx3 * Ix1 * Iy0 * Ip0 / ev2p
				sumx = sumx + sumx1 - sumx2 - sumx3
				
				sumy1 =	(uu**2-uy**2)/uu**3*uy
				sumy2 = uy*ux/uu**3*ux
				sumy3 = uy*up/uu**3*up

				sumy1 = sumy1 * Ix0 * Iy1 * Ip0/ev2y
				sumy2 = sumy2 *	Ix0 * Iy1 * Ip0/ev2x
				sumy3 = sumy3 * Ix0 * Iy1 * Ip0/ev2p
				sumy = sumy + sumy1 - sumy2 - sumy3
				
				sump1 =	(uu**2-up**2)/uu**3*up
				sump2 = up*ux/uu**3*ux
				sump3 = up*uy/uu**3*uy
				sump1 = sump1 * Ix0 * Iy0 * Ip1/ev2p
				sump2 = sump2 *	Ix0 * Iy0 * Ip1/ev2x
				sump3 = sump3 * Ix0 * Iy0 * Ip1/ev2y
				sump = sump + sump1 - sump2 - sump3
			    	
			end do
		end do
	end do
	
	!averate beam intensity n=integrate(f*f)/integrate(f)
	!<n>=n0/2/sqrt(2)  3-d gaussian
	! the ibs rate is always 2 times larger than B-M model. why?
	kvalue = gab * pnum/(2.0*pi)/sqrt(2.0*pi)/siges/sigex/sigey/2.0/sqrt(2.0)
	kvalue = kvalue/(2.0*pi)**3/ev2x/ev2y/ev2p
	!write(*,*) pnum,gamma
	sumindex = 8.0 * uspanx * uspany * uspanp / steps /steps /steps
	! taux=1/emitt*d(emitt)/dt = 2/sigmax*d(sigmax)/dt
	!<vx2>/2=sigmavx2	gaussian distribution
	alfxC = -sumx * kvalue * sumindex /gamma /ev2x/2.
	alfyC = -sumy * kvalue * sumindex /gamma /ev2y/2.
	alfpC = -sump * kvalue * sumindex /gamma /ev2p/2.

	return
	end



!------------------------------------------------------------------
	subroutine getInew(I0,I1,I2,aa,bb,cc,uux,uuy)
!with coupling
	parameter(pi=3.141592)
	real I0,I1,I2,aa,bb,uux,uuy,kk1,kk2

	kk1 = sqrt(4.*aa*bb-cc*cc)
	I0 = pi/kk1*exp(0.5*(-aa*uux**2-bb*uuy**2+cc*uux*uuy))

	I1 = -I0*uux/2.0
	I2 = -I0*uuy/2.0

	end


!------------------------------------------------------------------
!------------------------------------------------------------------	
	subroutine heatingIBSnew(alfxC,alfyC,alfpC,betax,dispx,steps)
!alfxC,alfyC,alfpC are emittance growth rates in per sec
!the vx, vy, vp are the real beam velocity, which are timed beta*clight already.
! calculation on ion beam for benchmark with IBS heating
! including dispersion and betafunction 4/11/2021
! assume betax=betay, alfx=alfy=0

	parameter(pi=3.141592,clight=2.998e8,qe=1.602e-19,dcool=200.)
	parameter(aatomi=1.,zatomi=1.)
	common/parms/pnum,zatom,aatom,epsx,epsy,gamma,sigs,dponp,coulomb

	real Ix0,Ix1,Ix2,Iy0,Iy1,Iy2,Ip0,Ip1,Ip2,sumx,sumy,sump
	real I0new,Ixnew,Ipnew
	real ev2x,ev2y,ev2p,ev2x0,ev2y0,ev2p0
	real sumx1,sumx2,sumx3,sumy1,sumy2,sumy3,sump1,sump2,sump3
	real kvalue,gab,me,betax,gammax,alfx,dispx,M1,M2,M3,rho
	integer steps

	beta=sqrt(1.-1./gamma**2)
	me = 9.10938356e-31
	epsilon = 8.854187e-12
	
	alfx = 0.0
	gammax = (1.+alfx**2)/betax	!correlation
	gab = (qe*qe/me)**2 * coulomb /4/pi/epsilon**2
	!gab = real(gab)	
	sigex0 = sqrt(epsx*betax)
	sigex = sqrt(epsx*betax + (dispx*dponp)**2)
	sigey = sqrt(epsy*betax)
	siges = sigs*gamma
	 
! sigma_vxo
	ev2x0 = epsx *gammax * gamma**2 *(beta*clight)**2 
	ev2y0 = epsy *gammax * gamma**2 *(beta*clight)**2 
	ev2p0 = dponp**2	*(beta*clight)**2*1.
!sigma_vx_new
	ev2x = epsx /betax * (1.+(alfx*dispx*dponp)**2/sigex**2) 
	ev2x = ev2x * gamma**2 *(beta*clight)**2
	ev2y = epsy /betax * gamma**2 *(beta*clight)**2  !noting:same beta x/y
	ev2p = dponp**2	*sigex0**2/sigex**2
	ev2p = ev2p *(beta*clight)**2

	rho = alfx*dispx*dponp/sqrt(sigex0**2+dispx**2*dponp**2*(1.+alfx**2))
	!ev2max = ev2x
	!if(ev2y.gt.ev2max) ev2max = ev2y
	!if(ev2p.gt.ev2max) ev2max = ev2p
	!uspan = 6.*sqrt(ev2max)

	M1 = 1.0/(1.-rho**2)/ev2x 
	M13 = rho/(1.-rho**2)/sqrt(ev2x*ev2p)
	M2 = 1.0/ev2y
	M3 = 1.0/(1.-rho**2)/ev2p
	
	uspanx = 4.*sqrt(ev2x+ev2x)
	uspany = 4.*sqrt(ev2y+ev2y)
	uspanp = 4.*sqrt(ev2p+ev2p)
	sumx=0.0
	sumy=0.0
	sump=0.0
	
	do i=1,steps
		do j=1,steps
			do k=1,steps
				
				ux = (2.*real(i+0.324)/real(steps) -1.) * uspanx
				uy = (2.*real(j+0.134)/real(steps) -1.) * uspany
				up = (2.*real(k+0.305)/real(steps) -1.) * uspanp

				uu = sqrt(ux**2 + uy**2 + up**2)

				
				call getI(Iy0,Iy1,Iy2,0.5/ev2y,0.5/ev2y,uy)
				call getInew(I0new,Ixnew,Ipnew,0.5/ev2x/(1.-rho**2),0.5/ev2p/(1.-rho**2),rho/(1.-rho**2)/sqrt(ev2x*ev2p),ux,up)

				sumx1 = (uu**2-ux**2)/uu**3
				sumx2 = ux*uy/uu**3
				sumx3 = ux*up/uu**3
					
				sumx1 = sumx1 * Ixnew * Iy0 *(M1*ux-M13*up)
				sumx2 = sumx2 *	Ixnew * Iy0 *M2*uy
				sumx3 = sumx3 * Ixnew * Iy0  *(M3*up-M13*ux)
				sumx = sumx + sumx1 - sumx2 - sumx3
				

				sumy1 =	(uu**2-uy**2)/uu**3
				sumy2 = uy*ux/uu**3
				sumy3 = uy*up/uu**3

				sumy1 = sumy1 * I0new * Iy1 *M2*uy
				sumy2 = sumy2 * I0new * Iy1 *(M1*ux-M13*up)
				sumy3 = sumy3 * I0new * Iy1 *(M3*up-M13*ux)
				sumy = sumy + sumy1 - sumy2 - sumy3
				
				sump1 =	(uu**2-up**2)/uu**3
				sump2 = up*ux/uu**3
				sump3 = up*uy/uu**3
				sump1 = sump1 * Iy0 * Ipnew *(M3*up-M13*ux)
				sump2 = sump2 *	Iy0 * Ipnew *(M1*ux-M13*up)
				sump3 = sump3 * Iy0 * Ipnew *M2*uy
				sump = sump + sump1 - sump2 - sump3
			    	
			end do
		end do
	end do
	
	!averate beam intensity n=integrate(f*f)/integrate(f)
	!<n>=n0/2/sqrt(2)  3-d gaussian
	! the ibs rate is always 2 times larger than B-M model. why?
	kvalue = gab * pnum/(2.0*pi)/sqrt(2.0*pi)/siges/sigex/sigey/2.0/sqrt(2.0)
	kvalue = kvalue/(2.0*pi)**3/ev2x/ev2y/ev2p/(1.-rho**2)
	sumindex = 8.0 * uspanx * uspany * uspanp / steps /steps /steps
	! taux=1/emitt*d(emitt)/dt = 2/sigmax*d(sigmax)/dt
	!<vx2>/2=sigmavx2	gaussian distribution
	alfxC = -sumx * kvalue * sumindex /gamma /ev2x0/2.*gammax*betax
	alfyC = -sumy * kvalue * sumindex /gamma /ev2y0/2.*gammax*betax
	alfpC = -sump * kvalue * sumindex /gamma /ev2p0/2.
	alfxC = alfxC + alfpC*gammax*dispx**2/epsx*dponp**2 

	return
	end
!------------------------------------------------------------------
!------------------------------------------------------------------	
	subroutine getheating111(alfxC,alfyC,alfpC,ionemittx,ionemitty,iondponp,ionsigs,ionnum,steps)
!alfx,alfy,alfp are emittance growth rates in per sec
!the vx, vy, vp are the real beam velocity, which are timed beta*clight already.
! calculation on ion beam for benchmark with IBS heating
	parameter(pi=3.141592,clight=2.998e8,qe=1.602e-19,dcool=200.)
	parameter(aatomi=1.,zatomi=1.)
	common/parms/pnum,zatom,aatom,epsx,epsy,gamma,sigs,dponp,coulomb
	common/betafunc/ betax0,betay0,scool,circ,betaxI0,betayI0
	real ionemittx,ionemitty,iondponp,ionsigs,ionnum
	real Ix0,Ix1,Ix2,Iy0,Iy1,Iy2,Ip0,Ip1,Ip2,sumx,sumy,sump
	real iv2x,iv2y,iv2p,ev2x,ev2y,ev2p
	real sumx1,sumx2,sumx3,kvalue,gab,me
	integer steps

	beta=sqrt(1.-1./gamma**2)
	me = 9.10938356e-31
	epsilon = 8.854187e-12
	qi = qe * zatomi
	gab = (qi*qe/me)**2 * coulomb /4/pi/epsilon**2
	!gab = real(gab)	
	sigix = sqrt(ionemittx*betaxI0)
	sigiy = sqrt(ionemittx*betayI0)
	sigis = ionsigs

! <vx2> <vy2> and <vp2> in rest frame
	iv2x = ionemittx / betaxI0 * gamma**2  *(beta*clight)**2	*1.
	iv2y = ionemitty / betayI0 * gamma**2 *(beta*clight)**2 *1.
	iv2p = iondponp**2  *(beta*clight)**2*1.
	
	!iv2max = iv2x
	!if(iv2y.gt.iv2max) iv2max = iv2y
	!if(iv2p.gt.iv2max) iv2max = iv2p
	
	ev2x = epsx / betax0 * gamma**2	*(beta*clight)**2  *1.
	ev2y = epsy / betay0 * gamma**2	*(beta*clight)**2 *1.
	ev2p = dponp**2	*(beta*clight)**2*1.
	!ev2max = ev2x
	!if(ev2y.gt.ev2max) ev2max = ev2y
	!if(ev2p.gt.ev2max) ev2max = ev2p

	uspanx = 6.*sqrt(iv2x+ev2x)
	uspany = 6.*sqrt(iv2y+ev2y)
	uspanp = 6.*sqrt(iv2p+ev2p)
	sumx=0.0
	sumy=0.0
	sump=0.0
	
	do i=1,steps
		ux = (2.*real(i+0.13131)/real(steps) -1.) * uspanx
		do j=1,steps
			uy = (2.*real(j-0.23133)/real(steps) -1.) * uspany
			do k=1,steps
				
				up = (2.*real(k+0.13133)/real(steps) -1.) * uspanp
				uu = sqrt(ux**2 + uy**2 + up**2)

				call getI(Ix0,Ix1,Ix2,0.5/ev2x,0.5/iv2x,ux)				
				call getI(Iy0,Iy1,Iy2,0.5/ev2y,0.5/iv2y,uy)
				call getI(Ip0,Ip1,Ip2,0.5/ev2p,0.5/iv2p,up)

				sumx1 = (uu**2-ux**2)/uu**3
				sumx2 = ux*uy/uu**3
				sumx3 = ux*up/uu**3
					
				sumx1 = sumx1 * Ix2 * Iy0 * Ip0/iv2x
				sumx2 = sumx2 *	Ix1 * Iy1 * Ip0/iv2y
				sumx3 = sumx3 * Ix1 * Iy0 * Ip1/iv2p
				sumx = sumx + sumx1 - sumx2 - sumx3
				!write(*,*) sumx,sumx1,sumx2,sumx3,sumx1 - sumx2 - sumx3
				
				sumy1 =	(uu**2-uy**2)/uu**3
				sumy2 = uy*ux/uu**3
				sumy3 = uy*up/uu**3
				sumy1 = sumy1 * Ix0 * Iy2 * Ip0/iv2y
				sumy2 = sumy2 *	Ix1 * Iy1 * Ip0/iv2x
				sumy3 = sumy3 * Ix0 * Iy1 * Ip1/iv2p
				sumy = sumy + sumy1 - sumy2 - sumy3
				
				sump1 =	(uu**2-up**2)/uu**3
				sump2 = up*ux/uu**3
				sump3 = up*uy/uu**3
				sump1 = sump1 * Ix0 * Iy0 * Ip2/iv2p
				sump2 = sump2 *	Ix1 * Iy0 * Ip1/iv2x
				sump3 = sump3 * Ix0 * Iy1 * Ip1/iv2y
				sump = sump + sump1 - sump2 - sump3
			    
			end do
		end do
	end do
	
	kvalue = gab * ionnum/(2.0*pi)/sqrt(2.0*pi)/sigis/sigix/sigiy/gamma
	kvalue = kvalue/(2.0*pi)**3/sqrt(iv2x*iv2y*iv2p)/sqrt(ev2x*ev2y*ev2p)
	
	sumindex = 8.0 * uspanx * uspany * uspanp / steps /steps /steps
	! taux=1/emitt*d(emitt)/dt = 2/sigmax*d(sigmax)/dt
	!<vx2>/2=sigmavx2	gaussian distribution
	alfxC = -sumx * kvalue * sumindex /gamma /ev2x/1836./2.
	alfyC = -sumy * kvalue * sumindex /gamma /ev2y/1836./2.
	alfpC = -sump * kvalue * sumindex /gamma /ev2p/1836./2.
	
	return
	end




!------------------------------------------------------------------
!------------------------------------------------------------------	
	subroutine getheating112(alfxC,alfyC,alfpC,ionemittx,ionemitty,iondponp,ionsigs,ionnum,steps)
!alfx,alfy,alfp are emittance growth rates in per sec
!the vx, vy, vp are the real beam velocity, which are timed beta*clight already.

	parameter(pi=3.141592,clight=2.998e8,qe=1.602e-19,dcool=200.)
	parameter(aatomi=1.,zatomi=1.)
	common/parms/pnum,zatom,aatom,epsx,epsy,gamma,sigs,dponp,coulomb
	common/betafunc/ betax0,betay0,scool,betaxI0,betayI0
	real ionemittx,ionemitty,iondponp,ionsigs,ionnum
	real Ix0,Ix1,Ix2,Iy0,Iy1,Iy2,Ip0,Ip1,Ip2,sumx,sumy,sump
	real iv2x,iv2y,iv2p,ev2x,ev2y,ev2p
	real sumx1,sumx2,sumx3,kvalue,gab,me
	integer steps

	beta=sqrt(1.-1./gamma**2)
	me = 9.10938356e-31
	epsilon = 8.854187e-12
	qi = qe * zatomi
	gab = (qi*qe/me)**2 * coulomb /4/pi/epsilon**2
	!gab = real(gab)	
	sigix = sqrt(ionemittx*betaxI0)
	sigiy = sqrt(ionemittx*betayI0)
	sigis = ionsigs

! <vx2> <vy2> and <vp2> in rest frame
	iv2x = ionemittx / betaxI0 * gamma**2  *(beta*clight)**2	*1.
	iv2y = ionemitty / betayI0 * gamma**2 *(beta*clight)**2 *1.
	iv2p = iondponp**2  *(beta*clight)**2*1.
	
	!iv2max = iv2x
	!if(iv2y.gt.iv2max) iv2max = iv2y
	!if(iv2p.gt.iv2max) iv2max = iv2p
	
	ev2x = epsx / betax0 * gamma**2	*(beta*clight)**2  *1.
	ev2y = epsy / betay0 * gamma**2	*(beta*clight)**2 *1.
	ev2p = dponp**2	*(beta*clight)**2*1.
	!ev2max = ev2x
	!if(ev2y.gt.ev2max) ev2max = ev2y
	!if(ev2p.gt.ev2max) ev2max = ev2p

	uspanx = 6.*sqrt(iv2x+ev2x)
	uspany = 6.*sqrt(iv2y+ev2y)
	uspanp = 6.*sqrt(iv2p+ev2p)
	sumx=0.0
	sumy=0.0
	sump=0.0
	
	do i=1,steps
		ux = (2.*real(i+0.13131)/real(steps) -1.) * uspanx
		do j=1,steps
			uy = (2.*real(j-0.23133)/real(steps) -1.) * uspany
			do k=1,steps
				
				up = (2.*real(k+0.13133)/real(steps) -1.) * uspanp
				uu = sqrt(ux**2 + uy**2 + up**2)

				call getI(Ix0,Ix1,Ix2,0.5/ev2x,0.5/iv2x,ux)				
				call getI(Iy0,Iy1,Iy2,0.5/ev2y,0.5/iv2y,uy)
				call getI(Ip0,Ip1,Ip2,0.5/ev2p,0.5/iv2p,up)

				sumx1 = (uu**2-ux**2)/uu**3*ux
				sumx2 = ux*uy/uu**3*uy
				sumx3 = ux*up/uu**3*up
					
				sumx1 = sumx1 * Ix1 * Iy0 * Ip0/iv2x
				sumx2 = sumx2 *	Ix1 * Iy0 * Ip0/iv2y
				sumx3 = sumx3 * Ix1 * Iy0 * Ip0/iv2p
				sumx = sumx + sumx1 - sumx2 - sumx3
				!write(*,*) sumx,sumx1,sumx2,sumx3,sumx1 - sumx2 - sumx3
				
				sumy1 =	(uu**2-uy**2)/uu**3*uy
				sumy2 = uy*ux/uu**3*ux
				sumy3 = uy*up/uu**3*up
				sumy1 = sumy1 * Ix0 * Iy1 * Ip0/iv2y
				sumy2 = sumy2 *	Ix0 * Iy1 * Ip0/iv2x
				sumy3 = sumy3 * Ix0 * Iy1 * Ip0/iv2p
				sumy = sumy + sumy1 - sumy2 - sumy3
				
				sump1 =	(uu**2-up**2)/uu**3*up
				sump2 = up*ux/uu**3*ux
				sump3 = up*uy/uu**3*uy
				sump1 = sump1 * Ix0 * Iy0 * Ip1/iv2p
				sump2 = sump2 *	Ix0 * Iy0 * Ip1/iv2x
				sump3 = sump3 * Ix0 * Iy0 * Ip1/iv2y
				sump = sump + sump1 - sump2 - sump3
			    
			end do
		end do
	end do
	
	kvalue = gab * ionnum/(2.0*pi)/sqrt(2.0*pi)/sigis/sigix/sigiy/gamma
	kvalue = kvalue/(2.0*pi)**3/sqrt(iv2x*iv2y*iv2p)/sqrt(ev2x*ev2y*ev2p)
	
	sumindex = 8.0 * uspanx * uspany * uspanp / steps /steps /steps
	! taux=1/emitt*d(emitt)/dt = 2/sigmax*d(sigmax)/dt
	!<vx2>/2=sigmavx2	gaussian distribution
	alfxC = -sumx * kvalue * sumindex /gamma /ev2x/1836./2.
	alfyC = -sumy * kvalue * sumindex /gamma /ev2y/1836./2.
	alfpC = -sump * kvalue * sumindex /gamma /ev2p/1836./2.
	
	return
	end

