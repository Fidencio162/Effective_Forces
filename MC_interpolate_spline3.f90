! The diameter of the particles is used as unit length.
! Fidencio Pérez-Hernández 25.09.2022
        program MC_Wolf
        implicit double precision(a-h,o-z)
        parameter(mp=1000,mr=2**11,nmq=2**6,nvq=200,ns=879)
        dimension x(mp),y(mp),z(mp),zz(mp),diam(mp)
        dimension qi(nmq),Sq(nmq),s(nmq)
        dimension qx(nmq,nvq),qy(nmq,nvq),qz(nmq,nvq)
        dimension r(mr),g(mr)
        dimension ueff(ns),rr(ns),y2(ns)
        common/box/boxl,rc,np
        common/electric/np1,np2,nz1,nz2,alpha,sigma1d,sigma2d,sigma12d
        common/parameters/rho

        pi = 4.d0 * datan(1.d0)

        np = 8**3
        sigmad = 10.d0

        ! packing fraction
        rho = 0.00019099d0
        d = (1.d0 / rho)**(1./3.)

        ! box length
        boxl = (dfloat(np) / rho)**(1.d0 / 3.d0)
        rc = boxl / 2.d0
        dr = rc / float(mr)

        print *, 'The length of the box is: ', boxl
        
        alpha = 4.d0 / boxl
        !xi=1.d0/101.d0

        phi = pi * rho * sigmad**3 / 6.d0

        print *, 'phi=', phi
        print *, 'rhoT=', rho
        !stop
        call iniconf(x, y, z, d)
        !! Numerical potential u_eff
        open(unit=1,file='ueffZ60.dat')
        do j=1, ns
         read(1,*) rr(j), ueff(j)
        enddo
        y1=(ueff(2)-ueff(1))/(rr(2)-rr(1))
        yn=0.d0
        call spline(rr, ueff, ns, y1, yn, y2)
        open(10, file='iniconf.xyz', status='unknown')
        write(10, 49) np
        write(10, 50)
        do i = 1, np
            write(10, 51) x(i), y(i), z(i)
        enddo
        close(10)

        open(30, file='energyMC.dat', status='unknown')
        del = 0.1d0
        nattemp = 0
        nacc = 1
        iseed = 123456789

        do i = 1, 10000000
        call mcmove(x, y, z, rr, ueff, y2, ener, nattemp, nacc, del, iseed)
            call adjust(nattemp, nacc, del)
            if (mod(i, 100) .eq. 0) write(30, *) i*1.d0, ener / np
            if (mod(i, 1000) .eq. 0) print*, i, ener / np, del
        enddo

        print *, 'The system has thermalized :)'

        open(20, file='finalconf.xyz', status='unknown')
        do i = 1, np
            write(20, *) x(i), y(i), z(i)
        enddo
        close(20)

        open(unit=13, file='gr_bb.dat')
        nacco = nacc
        ncp = 0
        do i = 1, mr
            g(i) = 0.d0
        enddo
        print *, 'Calculating g(r)...'
        do j = 1, 5000000
        call mcmove(x, y, z, rr, ueff, y2, ener, nattemp, nacc, del, iseed)
            call adjust(nattemp, nacc, del)
            if (mod(j, np) .eq. 0) then
              call average(x,y,z,g,nattemp,nacc,del,dr,iseed)
              ncp = ncp + 1
            endif
        enddo

        print *, 'Average number:', ncp

        ! Calculating g(r)
        do i = 2, mr
            r(i) = (i-1) * dr
            dv = (4.d0 * pi * r(i)**2 * dr) * rho
            g(i) = g(i) / (np * ncp * dv)
            write(13, *) r(i), g(i)
        enddo
        close(13)

        print *, 'Calculating g(r) finished'
        print *, 'Program finished!'

50      format(x, 'CONFIGURACION INICIAL')
51      format(3X, 'C', 3f15.7)
52      format(3X, 'N', 3f15.7)
49      format(2x, i5)
98      format(x, 'CONFIGURACION', 11x, i9)
99      format(20x, i5)
100     format(3f15.3)
200     format(2f15.7)
201 	 format(5x,'1',2x,'Sol',4x,'C',2x,i5,3f8.3,2x,'0.0000',2x, &
         '0.0000',2x,'0.0000')
202 	 format(5x,'1',2x,'Sol',4x,'N',2x,i5,3f8.3,2x,'0.0000',2x, &
         '0.0000',2x,'0.0000')
        end program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This subroutine calculates the initial configuration in 3D.		 
      	subroutine iniconf(xc,yc,zc,d)
      	implicit double precision(a-h,o-z)
	  parameter(mp=1000,mr=2**11,nmq=2**6,nvq=200)
	  dimension xc(mp),yc(mp),zc(mp)
	  dimension xr(mp),yr(mp),zr(mp)
	  common/box/boxl,rc,np

	  xc(1)=-(boxl-d)/2.d0
	  yc(1)=-(boxl-d)/2.d0
	  zc(1)=-(boxl-d)/2.d0
     	do i=2,np
         xc(i)=xc(i-1)+d
         yc(i)=yc(i-1)
         zc(i)=zc(i-1)
         if (xc(i) .gt. boxl/2.) then
            xc(i)=xc(1)
            yc(i)=yc(i-1)+d
            if (yc(i) .gt. boxl/2.) then
               xc(i)=xc(1)
               yc(i)=yc(1)
               zc(i)=zc(i-1)+d
            endif
         endif
      	enddo

	  return
      	end

! This configuration calculates the energy of a given configuration
	subroutine energy(x,y,z,rr,ueff,y2,xj,yj,zj,ener,j)
	  implicit double precision(a-h,o-z)
	  parameter(mp=1000,ns=879)
	dimension x(mp),y(mp),z(mp)
	dimension ueff(ns),rr(ns),y2(ns)
      	  common/box/boxl,rc,np
	  common/electric/np1,np2,nz1,nz2,alpha,sigma1d,sigma2d,sigma12d
	  common/parameters/rho
	  pi=4.d0*datan(1.d0)
	  ener=0.d0
      	do i=1,j-1
            !uij=0.d0
            xij=x(i)-xj
            yij=y(i)-yj
	    zij=z(i)-zj
            xij=xij-boxl*dnint(xij/boxl) 
            yij=yij-boxl*dnint(yij/boxl) 
	    zij=zij-boxl*dnint(zij/boxl)
            rij2=xij*xij+yij*yij+zij*zij
            rij=dsqrt(rij2) 
            if (rij .lt. rc) then
            if (rij .lt. rr(1)) then
              uij=1.d10
            else if (rij .gt. rr(ns)) then
              uij=0.d0
            else
               call Wolf(rr,ueff,y2,rij,uij)
            endif
            ener=ener+uij
           endif  
      	enddo
      	
      	 do i=j+1,np
            xij=x(i)-xj
            yij=y(i)-yj
	    zij=z(i)-zj
            xij=xij-boxl*dnint(xij/boxl) 
            yij=yij-boxl*dnint(yij/boxl) 
	    zij=zij-boxl*dnint(zij/boxl)
            rij2=xij*xij+yij*yij+zij*zij
            rij=dsqrt(rij2) 
            if (rij .lt. rc) then
            if (rij .lt. rr(1)) then
              uij=1.d10
            else if (rij .gt. rr(ns)) then
              uij=0.d0
            else
               call Wolf(rr,ueff,y2,rij,uij)
            endif
            ener=ener+uij
           endif 
      	enddo
	  return
	  end

! This subroutine calculates the pair potential between particles i & j
        subroutine Wolf(rr,ueff,y2,rij, uij)
          implicit double precision(a-h,o-z)
          parameter(mp=1000, mr=2**11, nmq=2**6, nvq=200,ns=879)
          dimension ueff(ns),rr(ns),y2(ns)
          common/electric/np1,np2,nz1,nz2,alpha,sigma1d,sigma2d,sigma12d
          common/box/boxl, rc, np
          call splint(rr, ueff, y2, ns, rij, y_final)          
          ! Calculate the pair potential uij
          uij = y_final
          !print*,rij,uij
          return
        end

! This subroutine displace the system to a new configuration
	  subroutine mcmove(x,y,z,rr,ueff,y2,ener,nattemp,nacc,del,iseed)
	  implicit double precision(a-h,o-z)
	  parameter(mp=1000,ns=879)
	  dimension x(mp),y(mp),z(mp),zz(mp),diam(mp)
	  dimension ueff(ns),rr(ns),y2(ns)
      	  common/box/boxl,rc,np
      	  
	  nattemp=nattemp+1
	  no=int(ranf(iseed)*np)+1
	  xo=x(no)
	  yo=y(no)
	  zo=z(no)
	  call energy(x,y,z,rr,ueff,y2,xo,yo,zo,enero,no)
	  x(no)=x(no)+(ranf(iseed)-0.5d0)*del
	  y(no)=y(no)+(ranf(iseed)-0.5d0)*del
	  z(no)=z(no)+(ranf(iseed)-0.5d0)*del
! periodic boundary conditions
	  x(no)=x(no)-boxl*dnint(x(no)/boxl)
	  y(no)=y(no)-boxl*dnint(y(no)/boxl)
	  z(no)=z(no)-boxl*dnint(z(no)/boxl)
	  call energy(x,y,z,rr,ueff,y2,x(no),y(no),z(no),enern,no)
	  if (ranf(iseed) .lt. dexp(-(enern-enero))) then
	     ener=ener+0.5d0*(enern-enero)
	     nacc=nacc+1
	     else
	     x(no)=xo
	     y(no)=yo
	     z(no)=zo
	  endif
	  return
	  end

	 subroutine average(x,y,z,g,nattemp,nacc,del,dr,iseed)
	 implicit double precision(a-h,o-z)
	 parameter(mp=1000,mr=2**11)
        dimension x(mp),y(mp),z(mp),g(mr)
	 common/box/boxl,rc,np
	 common/electric/np1,np2,nz1,nz2,alpha,sigma1d,sigma2d,sigma12d
	 nattemp=nattemp+1
	 nbin=0
       ! g11
       do i=1,np-1
          do j=i+1,np
             xij=x(j)-x(i)
             yij=y(j)-y(i)
             zij=z(j)-z(i)
             xij=xij-boxl*dnint(xij/boxl)
             yij=yij-boxl*dnint(yij/boxl)
             zij=zij-boxl*dnint(zij/boxl)
             rij2=xij*xij+yij*yij+zij*zij
             rij=dsqrt(rij2)
             if (rij .lt. rc) then
                nbin=dint(rij/dr)
                if (nbin .lt. mr) then
                   g(nbin)=g(nbin)+2.
                endif
             endif
          enddo
       enddo
	     nacc=nacc+1
	  return
	  end   

! Random generator algorithm
! Numerical Recipes
      FUNCTION ranf(Idum)
      implicit double precision(a-h,o-z)

      PARAMETER (IM1=2147483563, IM2=2147483399, &
      AM=1./IM1, IMM1=IM1-1,IA1=40014, IA2=40692, &
      IQ1=53668, IQ2=52774, IR1=12211,IR2=3791, &
      NTAb=32, NDIv=1+IMM1/NTAb, EPS=1.2E-7, RNMx=1.-EPS)
!!
	dimension iv(ntab)
      SAVE iv, iy, idum2
      DATA idum2/123456789/, iv/NTAb*0/, iy/0/
      IF (Idum.LE.0) THEN
         Idum = MAX(-Idum, 1)
         idum2 = Idum
         DO j = NTAb + 8, 1, -1
            k = Idum/IQ1
            Idum = IA1*(Idum-k*IQ1) - k*IR1
            IF (Idum.LT.0) Idum = Idum + IM1
            IF (j.LE.NTAb) iv(j) = Idum
         END DO
         iy = iv(1)
      END IF
      k = Idum/IQ1
      Idum = IA1*(Idum-k*IQ1) - k*IR1
      IF (Idum.LT.0) Idum = Idum + IM1
      k = idum2/IQ2
      idum2 = IA2*(idum2-k*IQ2) - k*IR2
      IF (idum2.LT.0) idum2 = idum2 + IM2
      j = 1 + iy/NDIv
      iy = iv(j) - idum2
      iv(j) = Idum
      IF (iy.LT.1) iy = iy + IMM1
      ranf = MIN(AM*iy, RNMx)
      RETURN
      END

! This subroutine adjusts the displacement of particles
	  subroutine adjust(nattemp,nacc,dr)
	  implicit double precision(a-h,o-z)
	  if (mod(nattemp,nacc) .eq. 0) then
	     ratio=real(nacc)/real(nattemp)
	     if (ratio .gt. 0.5) then
	        dr=dr*1.05
	     else
	        dr=dr*0.95
	     endif
	  endif
      return
      end
      
      SUBROUTINE spline(x, y, n, y1, yn, y2)
       ! =====================================================
       ! Input x and y=f(x), n (dimension of x,y), (Ordered)
       ! y1 and yn are the first derivatives of f in the 1st point and the n-th
       ! Output: array y2(n) containing second derivatives of f(x_i)
       ! =====================================================
       IMPLICIT NONE
       INTEGER:: n, i, j
       INTEGER, PARAMETER:: n_max = 500
       REAL*8, INTENT(in):: x(n), y(n), y1, yn
       REAL*8, INTENT(out):: y2(n)
       REAL*8:: p, qn, sig, un, u(n)

       IF (y1 > .99e30) THEN    ! natural spline conditions
         y2(1) = 0
         u(1) = 0
       ELSE
         y2(1) = -0.5
         u(1) = (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-y1)
       END IF

       DO i = 2, n-1                            ! tridiag. decomposition
         sig = (x(i)-(i-1))/(x(i+1)-x(i-1))
         p = sig*y2(i-1)+2.
         y2(i) = (sig-1.)/p
         u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
       END DO

       IF (yn > .99e30) THEN   ! natural spline conditions
         qn = 0
         un = 0
       ELSE 
         qn = 0.5
         un=(3./(x(n)-x(n-1)))*(yn-(y(n)-y(n-1))/(x(n)-x(n-1)))
       END IF

         y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

       DO j = n-1, 1, -1          !  backwards substitution tri-diagonale
         y2(j) = y2(j)*y2(j+1)+u(j)
       END DO
       RETURN

       END SUBROUTINE spline

       SUBROUTINE splint(x_in, y_in, spline_res, n, x_0, y_final)
       ! =====================================================
       ! Subroutine that does the actual interpolation
       ! Input arrays of x_in and y_in=f(x), spline_res is the result of 
       ! the 'spline' subroutine, x_0 is the corresponding value we are looking for
       ! i.e. (time_at_max in hubble), y_final is the output result
       ! =====================================================
       IMPLICIT NONE
       INTEGER:: n, k, k_low, k_high
       REAL*8, INTENT(in):: x_in(n), y_in(n), spline_res(n), x_0
       REAL*8, INTENT(out):: y_final
       REAL*8:: a, b, h
 
         k_low = 1
         k_high = n
 
  99   IF (k_high - k_low > 1) THEN
         k = (k_high + k_low) / 2
       IF (x_in(k) > x_0) THEN
         k_high = k
       ELSE
         k_low = k
       END IF
       GOTO 99
       ENDIF

         h = x_in(k_high) - x_in(k_low) 
       IF (h == 0) STOP "Bad x_in input"
         a = (x_in(k_high)-x_0)/h
         b = (x_0 - x_in(k_low))/h
         y_final = a*y_in(k_low)+b*y_in(k_high)+((a**3-a)*spline_res(k_low)+(b**3-b)*spline_res(k_high))*(h**2)/6

       RETURN 
 
       END SUBROUTINE splint
