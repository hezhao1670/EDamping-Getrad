 

subroutine MotoCarlo(alfxC,alfyC,alfpC,ionemittx,ionemitty,iondponp,ionsigs,ionnum,steps)
!alfx,alfy,alfp are emittance growth rates in per sec
!the vx, vy, vp are the real beam velocity, which are timed beta*clight already.

	parameter(pi=3.141592,clight=2.998e8,qe=1.602e-19,dcool=200.)
	parameter(aatomi=1.,zatomi=1.)
	common/parms/pnum,zatom,aatom,epsx,epsy,gamma,sigs,dponp,coulomb
	common/betafunc/ betax0,betay0,scool,betaxI0,betayI0
	real ionemittx,ionemitty,iondponp,ionsigs,ionnum
	real sumx,sumy,sump
	real iv2x,iv2y,iv2p,ev2x,ev2y,ev2p
	real sumx1,sumx2,sumx3,kvalue,gab,me
	integer steps,nparticles
	
	nparticles = steps

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
	
	ev2x = epsx / betaxI0 * gamma**2	*(beta*clight)**2  *1.
	ev2y = epsy / betayI0 * gamma**2	*(beta*clight)**2 *1.
	ev2p = dponp**2	*(beta*clight)**2*1.

	uspanx = 5.*sqrt(iv2x+ev2x)
	uspany = 5.*sqrt(iv2y+ev2y)
	uspanp = 5.*sqrt(iv2p+ev2p)
	!write(*,*) uspanx,uspany,uspanp
	
	do jj=1,5000	
		vix = gauss_random(0,sqrt(iv2x))
		viy = gauss_random(0,sqrt(iv2y))
		vip = gauss_random(0,sqrt(iv2p))
		sumx=0.0
		sumy=0.0
		sump=0.0
		do ii=1,nparticles
			vex = gauss_random(0,sqrt(ev2x))
			vey = gauss_random(0,sqrt(ev2y))
			vep = gauss_random(0,sqrt(ev2p))
			ux = vix - vex
			uy = viy - vey
			up = vip - vep
			uu = sqrt(ux**2 + uy**2 + up**2)
			
			sumx1 = sumx1 + (uu**2 - ux**2)/uu**3 * vex**2/sqrt(ev2x)
			sumy1 = sumy1 + ux*uy/uu**3 * vex*vey / sqrt(ev2y)
			sump1 = sump1 + ux*up/uu**3 * vex*vep / sqrt(ev2p)
			
		end do
		sumx1 = sumx1/nparticles
		sumy1 = sumy1/nparticles
		sump1 = sump1/nparticles
		sumtotal = sumx1-sumy1-sump1

		write(443,*) jj,sumx1,sumy1,sump1,sumtotal
	end do
	end

	
		

subroutine MotoCarlo1(alfxC,alfyC,alfpC,ionemittx,ionemitty,iondponp,ionsigs,ionnum,steps)
!alfx,alfy,alfp are emittance growth rates in per sec
!the vx, vy, vp are the real beam velocity, which are timed beta*clight already.

	parameter(pi=3.141592,clight=2.998e8,qe=1.602e-19,dcool=200.)
	parameter(aatomi=1.,zatomi=1.)
	common/parms/pnum,zatom,aatom,epsx,epsy,gamma,sigs,dponp,coulomb
	common/betafunc/ betax0,betay0,scool,betaxI0,betayI0
	real ionemittx,ionemitty,iondponp,ionsigs,ionnum
	real sumx,sumy,sump
	real iv2x,iv2y,iv2p,ev2x,ev2y,ev2p
	real sumx1,sumx2,sumx3,kvalue,gab,me
	integer steps,nparticles
	
	nparticles = steps

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
	
	ev2x = epsx / betaxI0 * gamma**2	*(beta*clight)**2  *1.
	ev2y = epsy / betayI0 * gamma**2	*(beta*clight)**2 *1.
	ev2p = dponp**2	*(beta*clight)**2*1.

	uspanx = 5.*sqrt(iv2x+ev2x)
	uspany = 5.*sqrt(iv2y+ev2y)
	uspanp = 5.*sqrt(iv2p+ev2p)
	!write(*,*) uspanx,uspany,uspanp
	
	do jj=1,5000	
		vex = gauss_random(0,sqrt(ev2x))
		vey = gauss_random(0,sqrt(ev2y))
		vep = gauss_random(0,sqrt(ev2p))
		sumx=0.0
		sumy=0.0
		sump=0.0
		do ii=1,nparticles
			
			vix = gauss_random(0,sqrt(iv2x))
			viy = gauss_random(0,sqrt(iv2y))
			vip = gauss_random(0,sqrt(iv2p))
			ux = vix - vex
			uy = viy - vey
			up = vip - vep
			uu = sqrt(ux**2 + uy**2 + up**2)
			
			sumx1 = sumx1 + (uu**2 - ux**2)/uu**3 * vex**2/sqrt(ev2x)
			sumy1 = sumy1 + ux*uy/uu**3 * vex*vey / sqrt(ev2y)
			sump1 = sump1 + ux*up/uu**3 * vex*vep / sqrt(ev2p)
			
		end do
		sumx1 = sumx1/nparticles
		sumy1 = sumy1/nparticles
		sump1 = sump1/nparticles
		sumtotal = sumx1-sumy1-sump1

		write(443,*) jj,sumx1,sumy1,sump1,sumtotal
	end do
	end

	
		

		

	FUNCTION gauss_random(average,sigma)
              REAL gauss_random1,gauss_random2
              REAL average,sigma,mm,nn,w
2             CALL random(xsl)
              mm=2.0*xsl-1.0
              CALL random(xsl)
              nn=2.0*xsl-1.0
              w=mm*mm+nn*nn;
              IF(w.gt.1.0) GOTO 2
              IF(w.eq.0.0) GOTO 2
              gauss_random=mm*sqrt((-2.0*log(w))/w)*sigma+average;
    	END



