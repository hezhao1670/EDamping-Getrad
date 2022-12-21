!234567 ifort -o getrad7 getrad7.f90 getIBS.f90 heating.f90 
!2020/7/17 revise the heating rate from ions. <vx^2>/2=sigma_vx^2
	program getrad6
	!include 'getIBS.f90'
	parameter(clight=2.998e8,erest = 0.511e6,pi=3.1415926536,qe=1.602e-19)
	parameter(nb=100000)
	common/ibslatt/bxa(nb),bxap(nb),byap(nb),bya(nb),dxa(nb),dxpa(nb),fraca(nb),dya(nb),dypa(nb)
	common/iparms/nslice,kloop,nsim,ibstype
	common/beam/dtsim,betas,betabend,frev
	common/parms/pnum,zatom,aatom,epsx,epsy,gamma,sigs,dponp,coulomb,r0
	common/rad/growp_rad,growx_rad,growy_rad,alphax_rad,alphay_rad,alphap_rad
	common/betafunc/ betax0,betay0,scool,circ,betaxI0,betayI0,dispex,dispix
	common/ionbeam/epsxI,epsyI,dponpI,sigsI,pnumI
	common/switch/ibson,heaton,dampon,intmeth

! given mad output this code calculates electron bunch
! parameters

	open(unit=10,file='para.in',status='unknown')
	!open(unit=10,file='getrad5_ering.in',status='unknown')
! beta functions in middle of cooling section, length of cooling section, momentum compaction
	read(10,*)betax0,betay0,scool,ibstype,slipfactor
! initial values
! lorentz factor,particles per bunch,  rms geometrical emittances, rms length, rms dponp
! simulations step size and number of steps
	read(10,*)gamma,pnum,epsx,epsy,sigs,dponp,dtsim,nsim

	read(10,*)betaxI0,betayI0,dispex,dispix
! ion beam parameters
! emittx,emity,dp/p,sigs,nion,steps
	read(10,*) epsxI,epsyI,dponpI,sigsI,pnumI
	zatom = 1		!charge
	aatom = 1/1836.109	!electron
	!aatom = 1		!ion
! synchrotron rad integrals from mad8 using definitiosn in helm et al pac73 p 900

!bunch length should be recalculated. sigRMS=6.6e-4--> sigs=0.07m (benchmark with Nagaitsev[MathCad])
!sigrms = 6.6E-4
!sigs = sigs * sigrms / dponp    
	write(*,*) "dponp & length:",dponp,sigs

	read(10,*)fsynch1,fsynch2,fsynch3,fsynch4x,fsynch4y,fsynch5x,fsynch5y
! scaling for synchrotron radiation to mimic a stronger undulator
	write(*,*)" estt"
	read(10,*) bratio,intmeth
	
	read(10,*) ibson,heaton,dampon
close(10)
	fsynch2 = fsynch2*bratio**2
	fsynch4x = fsynch4x*bratio**4
	fsynch4y = fsynch4y*bratio**4
	!open(unit=10,file='weasel35.dat',status='unknown')
	!open(unit=10,file='optics_ering.dat',status='unknown')
	open(unit=10,file='optics.dat',status='unknown')
! number of lines of optics file, number of lines to skip after first (15 nominally)
	read(10,*)ntotal,nskip
	write(*,*) ntotal,nskip
	do k=1,nskip
		read(10,*)
	enddo
	nslice = ntotal-1-nskip
	sold=0

	do k=1,nslice
		!read(10,*)s,betx,alfx,bety,dx,dpx
		read(10,*)s,betx,alfx,bety,alfy,dx,dpx,dy,dpy
		!read(10,*)s,betx,alfx,bety,alfy,dx,dpx
		bxa(k)=betx
		! bxap is d beta/ds
		bxap(k)=-2*alfx
		byap(k)=-2*alfy*1.0
		bya(k)=bety
		dxa(k)=dx
		dxpa(k)=dpx
		dya(k)=dy*1.0
		dypa(k)=dpy*1.0
		fraca(k)=s-sold
		sold=s
		!write(*,*) betx,alfx,bety,dx
	end do
	circ = sold
	write(*,*) 'Finish lattice read!'
	close(10)
	beta = sqrt(1-1/gamma**2)
	frev = beta*clight/circ
	re = 2.818e-15
	write(6,*)re,gamma,frev,fsynch2
	alpha0_rad = 0.3333*re*frev*fsynch2*gamma**3
	curlydx = fsynch4x/fsynch2
	curlydy = fsynch4y/fsynch2


	write(6,*)'alpha0_rad = ',alpha0_rad
! energy loss per turn in eV
eloss = 0.66666*erest*gamma**4*re*fsynch2

	write(6,*)' energy loss per turn = ',eloss
! get the emittance damping rates assuming nominal partition
	write(6,*)' curly D = ',curlydx,curlydy
! damping rates for emittance
	alphax_rad = 2*alpha0_rad*(1-curlydx)	!alphax_rad = 2*alphax
	alphay_rad = 2*alpha0_rad*(1-curlydy)
	alphap_rad = 4*alpha0_rad*(1+curlydx/2+curlydy/2)
! nominal dp/p and rms horizontal emit
! reduced compton wavelength of electron
	flamcomp=2.426e-12/(2*pi)
	dponp_rad2 = (55./(32*sqrt(3.)))*flamcomp*gamma**2 *fsynch3/(2*fsynch2+fsynch4x+fsynch4y)
	epsx_rad = (55./(32*sqrt(3.)))*flamcomp*gamma**2 *fsynch5x/(fsynch2-fsynch4x)
	epsy_rad = (55./(32*sqrt(3.)))*flamcomp*gamma**2 *fsynch5y/(fsynch2-fsynch4y)

! d epsx / dt = - alphax_rad*eps + growx_rad
! d epsx / dt = -2 * alphax *eps + growx_rad
	growx_rad = alphax_rad*epsx_rad*dampon
	growy_rad = alphay_rad*epsy_rad*dampon
	write(6,*)' alfarad x,y,p = ', alphax_rad, alphay_rad,alphap_rad

! rf parameter
	betas = sigs/dponp
! rms longitudinal emittance, sigs = sqrt(epsp*betas), dponp = sqrt(epsp/betas)
	epsp = sigs*dponp

! d dponp**2 /dt = -alphap_rad *dponp**2 + growp_rad
	growp_rad = alphap_rad*dponp_rad2*dampon

	open(unit=60,file='emitt.dat',status='unknown')
	open(unit=66,file='ibs.dat',status='unknown')
	open(unit=67,file='ibs_lattice.dat',status='unknown')

!classical radius
	r0 = (zatom**2/aatom)*1.535e-18	
	!r0 = 2.8e-15
!############ Loop for calculate the emittance ###############
	do kloop=1,nsim
	! average the emittances assuming strong coupling
		!epsx = 0.5*(epsx+epsy)
		!epsy=epsx
		call updateemit
		write(60,100)kloop*dtsim*frev,epsx,epsy,sigs,dponp

	enddo
	write(6,100)kloop*dtsim*frev,epsx,epsy,sigs,dponp
	

100   format(10e14.6)
! get the cooling rates

! normalized rms emittance of electrons
	emit=epsx*gamma

	dgamma = dponp*gamma
	re = 2.828e-15
	rp = re/1836
	frac = scool/3834.
	frac2= scool/circ

! longitudinal number density in lab frame
	betacool = sqrt(betax0*betay0)
	dene = pnum/((2*pi)**1.5*sigs*sqrt(epsx*epsy)*betacool)
! rms transverse speed in beam frame
	dbetax = gamma*sqrt(epsx/betacool)
        dbetaY = gamma*sqrt(epsy/betacool)

	write(6,*)'rms beam radius = ',sqrt(epsx*betacool)
! Alexei's rates for protons
	write(6,*)'coulomb = ',coulomb
	coolz = 4*pi*frac*coulomb*dene*clight*re*rp/(gamma**2 * dponp**2 * dbetax)
	coolx = 4*pi*frac*coulomb*dene*clight*re*rp/(gamma**2 * dponp * dbetax**2)
	write(6,*)' Alexei coolz,coolx',coolz,coolx
! Nagaitsev's rates for cold longitudinal beams
	coolzn = coolz*(10./12.56)*dponp/dbetax
	coolxn = coolzn*(0.25*pi)*dponp/dbetax


	write(6,*)'Nagaitsev coolz,coolx',coolzn,coolxn
! cooling time described in pcdr
	coolmmb =  3.342*frac*coulomb*dene*clight*re*rp/(gamma**2 * dponp * dbetax*dbetay)
	write(6,*)' spherical cooling time (s) =',1/coolmmb


! space charge tune shift

!	dqsc = pnum*re*circ/(sigs*sqrt(2*pi)*4*pi*epsx*beta**2 * gamma**3)

	dqscx = pnum*re*circ/(sigs*sqrt(2*pi)*2*pi*sqrt(epsx)*(sqrt(epsx)+sqrt(epsy))*beta**2 * gamma**3)
	dqscy = pnum*re*circ/(sigs*sqrt(2*pi)*2*pi*sqrt(epsy)*(sqrt(epsx)+sqrt(epsy))*beta**2 * gamma**3)

	write(6,*)' small amplitude electron space charge tune shifts = ',dqscx,dqscy

	write(6,*)' rms radii in cooling section = ',sqrt(betacool*epsx),sqrt(betacool*epsy)
	write(6,*)' sigurest2,sig(betax_rest) =  ',dbetax
	write(6,*)' sigurest,dponp = ',dponp
	write(6,*)' qbeamcool = ',pnum*1.602e-19
	write(6,*)' slencool = ', scool
	write(6,*)' ebeams = ', sigs
	peaki = 3.e8*pnum*1.6e-19/(sigs*sqrt(2*pi))
	write(6,*)' Ipeak = ',peaki
! get threshold z/n
	zdivn = slipfactor*2*pi*gamma*erest*dponp**2/peaki
	write(6,*)'threshold z/n = ',zdivn
! get space charge tune shift of ions due to electrons
! electron peak current
	peaki = pnum*clight*1.602e-19/(sqrt(2*pi)*sigs)
! average beta function
	betai=100

	etotdq = 940.e6*gamma

	dqsc = 377*peaki*scool*betai/((4*pi)**2 *betacool*epsx*gamma**2 *etotdq)

	write(6,*)'ion tune shift = ',dqsc
! voltage on 225 MHz
	etote = erest*gamma
	epshat =dponp*etote
	tauhat = sigs/3.e8
	t0 = circ/3.e8

	vdot = (epshat/tauhat)**2 *(t0*slipfactor)/etote
	volt = vdot/(2*pi*225.e6)
	write(6,*)' electron containment voltage = ',volt


	stop
	end

!-------------------------------------------------------------
	subroutine updateemit
	! updates emittances for one dt
	parameter(clight=2.998e8,erest = 0.511e6,pi=3.1415926536,qe=1.602e-19)
	parameter(nb=100000)
	common/ibslatt/bxa(nb),bxap(nb),byap(nb),bya(nb),dxa(nb),dxpa(nb),fraca(nb),dya(nb),dypa(nb)
	common/iparms/nslice,kloop,nsim,ibstype
	common/beam/dtsim,betas,frev
	common/parms/pnum,zatom,aatom,epsx,epsy,gamma,sigs,dponp,coulomb
	common/rad/growp_rad,growx_rad,growy_rad,alphax_rad,alphay_rad,alphap_rad
	common/betafunc/ betax0,betay0,scool,circ,betaxI0,betayI0,dispex,dispix
	common/ionbeam/epsxI,epsyI,dponpI,sigsI,pnumI
	common/switch/ibson,heaton,dampon,intmeth
	frac2= scool/circ

! write(6,*)' in updateemit',nslice
	sumx=0
	sumy=0
	sump=0
	den=0
	beta = sqrt(1-1/gamma**2)
	!circ = 429.727834240
	frev = beta*clight/circ
	do k=1,nslice
        ! read(10,*)bx,bxp,by,dx,dxp,frac
		bx = bxa(k)
		! bxp is d beta/ds
		bxp=bxap(k)
		byp=byap(k)
		by = bya(k)
		dx = dxa(k)
		dxp = dxpa(k)
		frac = fraca(k)
		dy = dya(k)
		dyp = dypa(k)
		!write(6,100)bx,bxp,dx,dxp,by,frac
! get intrabeam scattering rates for emittances
!!! getrates:Nagaitsev without vertical dispertion
!!! getrates2:Borjken-Mtingua Model with vertical dispersion (Betacool)
!!! getrates3:Nagaitsev with vertical dispertion
		call getrates3(alfx,alfy,alfp,bx,bxp,by,byp,dx,dxp,dy,dyp,ibstype)
		sumx = sumx + frac*alfx
		sumy = sumy + frac*alfy
		sump = sump + frac*alfp
		den = den+frac
!save the IBS heating rate
		if(kloop.eq.1.)write(67,100) frac,den,bx,by,dx,dy,alfx,alfy,alfp,sumx,sumy,sump
	enddo
    ! normalize
	sumx = sumx/den
	sumy = sumy/den
	sump = sump/den
         !     write(6,*)' emittance growth rates = ',sumx,sumy,sump
	if(kloop.eq.nsim)write(6,*)' den = ',den
    ! rates are  for emittance growth

! get the electron beam heating rate by ions 
!!! getheating0:Spitzer simple formula
!!! getheating11,12:Based on Bolzmann Transportation equation and coulomb collision operator
!!! getheating111,112: ion cooling from electron in cooling section, for testing the cooling 
!!! heatingIBS: IBS calculation only in cooling section, for benchmark
!!! heatingIBSnew: IBS calculation with dispersion only in cooling section, for benchmark
!!! getheating13:heating with De and Di

	!call getheating0(teqx,teqy,teqs,epsxI,epsyI,dponpI,sigsI,pnumI)
	!call getheating11(heatx,heaty,heatp,epsxI,epsyI,dponpI,sigsI,pnumI,200)
	!write(*,*) "heating11: ",heatx,heaty,heatp
	!call getheating12(heatx,heaty,heatp,epsxI,epsyI,dponpI,sigsI,pnumI,200)
	!write(*,*) "heating12: ",heatx,heaty,heatp
	call getheating13(heatx,heaty,heatp,epsxI,epsyI,dponpI,sigsI,pnumI,200)
	!write(6,*) "heating13: ",dispex,dispix,heatx,heaty,heatp


	!call getheating111(heatx1,heaty1,heatp1,epsxI,epsyI,dponpI,sigsI,pnumI,500)
	!call getheating112(heatx2,heaty2,heatp2,epsxI,epsyI,dponpI,sigsI,pnumI,500)
	!call heatingIBS(heatxIBS,heatyIBS,heatpIBS,295.0,300)
	!write(*,*) "heatingIBS:",heatxIBS,heatyIBS,heatpIBS

	!call heatingIBSnew(heatxIBS,heatyIBS,heatpIBS,295.0,5.0,400)
	!write(*,*) "heatingIBSnew:",heatxIBS,heatyIBS,heatpIBS
	!write(*,*) "part1",heatx,heaty,heatp
	!write(*,*) "part2-1",heatx1,heaty1,heatp1

	!call MotoCarlo(heatx,heaty,heatp,epsxI,epsyI,dponpI,sigsI,pnumI,20000)

	kvalue=kb/me/gamma**2/beta**2/clight**2
	kk = -teq*log(temI/kb-teme/kb)

	!epsx = epsx + epsx*dtsim*heatx
	!epsy = epsy + epsy*dtsim*heaty
	!dponp = dponp + dponp*dtsim*heatp/2.
	

	if(intmeth == 0) then
		epsx = epsx*(1.0 + dtsim*sumx*ibson + dtsim*heatx*frac2*heaton - dtsim*alphax_rad*dampon) + 1.0*growx_rad*dtsim
		epsy = epsy*(1.0 + dtsim*sumy*ibson + dtsim*heaty*frac2*heaton - dtsim*alphay_rad*dampon) + 1.0*growy_rad*dtsim
		dponp = dponp*(1.0 + 0.5*sump*dtsim*ibson + 0.5*dtsim*heatp*frac2*heaton - dtsim*0.5*alphap_rad*dampon) +dtsim*0.5*growp_rad/dponp*1.0
	else
	! here to speed up the calculation, which is based on the estimate of the final emittance at t=infinity,
	!d eps/dt = - rad_damping *eps + rad_excitation + ibs_grwoth(eps) +ion_heating(eps)

	!if the last 2 terms would not depend on eps, there is an analytical solution for t=infinity. But the dependence on eps is small against rad_damping.

	!I iterate by assuming the terms are constant and set eps_new = (eps_old + eps_analytical)/2. This converges in less than 30 iterations, reducing the run time from hours to minutes. I can use the time savings to cut the slices finer.

		epsx =  0.5*( sumx*ibson*epsx + heatx*frac2*heaton*epsx  + growx_rad) / alphax_rad + 0.5*epsx
		epsy =  0.5*( sumy*ibson*epsy + heaty*frac2*heaton*epsy  + growy_rad) / alphay_rad + 0.5*epsy
		dponp = 0.5*( sump*ibson*dponp + heatp*frac2*heaton*dponp  + growp_rad/dponp ) / alphap_rad + 0.5*dponp
		if(epsx*epsy*dponp.lt.0.0) then
			epsx =  abs(epsx)
			epsy =  abs(epsy)
			dponp = abs(dponp)
		endif
			
	endif


	!epsx = epsx*(1.0 + dtsim*sumx*ibson + dtsim*heatx*frac2*heaton - dtsim*alphax_rad*dampon) + 1.0*growx_rad*dtsim
	!epsy = epsy*(1.0 + dtsim*sumy*ibson + dtsim*heaty*frac2*heaton - dtsim*alphay_rad*dampon) + 1.0*growy_rad*dtsim
	!dponp = dponp*(1.0 + 0.5*sump*dtsim*ibson + 0.5*dtsim*heatp*frac2*heaton - dtsim*0.5*alphap_rad*dampon)+dtsim*0.5*growp_rad/dponp*1.0
100 format(20e14.6)
	

	sigs = betas*dponp
	!write(6,*)alphax_rad,alphay_rad,growx_rad
    !write(66,101) kloop,1./sumx,1./sumy,2./sump,1./teqx,1./teqy,2./teqs,1./heatx,1./heaty,2./heatp,"unit=s"
	write(66,101) kloop*dtsim*frev,sumx,sumy,sump,alphax_rad/2.0,alphay_rad/2.0,alphap_rad/2.0,heatx*scool/circ,heaty*scool/circ,heatp*scool/circ,"unit=s^-1"

101 format(10(e13.3),A10)
	return
	end


