c****************************************************************************
c Shooting code based on HYBRD for homoclinic orbits of the dimensionless
c equations for a rod in a magnetic field.
c
c The order of the quantities in the integration is:
c
c     y(1)=theta    y(2)=psi    y(3)=p_theta    y(4)=p_psi    
c
c The order of the unknown quantities in the shooting method is:
c
c           x(1)=xl     x(2)=delta
c
c     April 2013 (Gert van der Heijden)
c
c****************************************************************************
      program shoot
      implicit real*8 (a-h,o-z)
c
      parameter (node=4,nrd=4)                           ! DOPRI (integrator)
      parameter (lws=11*node+8*nrd+21,liws=nrd+21)
      dimension ws(lws),iws(liws)
c
      parameter (nx=2,lnx=(nx*(nx+1))/2)              ! HYBRD (Newton solver)
      dimension diag(nx),fjac(nx,nx),r(lnx),qtf(nx)
      dimension wa1(nx),wa2(nx),wa3(nx),wa4(nx)
c
      dimension x(nx),fv(nx),ys(node),par(20),con(8*node),icomp(node)
      common /init/xl,delta,fp(4)
      common /pars/xepirou,xlamda,vv,xalpha,xmu
      external fcn,f,output
      logical newton
c
      newton=.true.                                      ! Newton iterations?
c
      open(unit=1,file='shoot4d.dat',status='unknown')
      read(1,*)xepirou,xlamda,vv,xalpha,xmu                      ! parameters
      read(1,*)(x(j),j=1,nx)                                 ! guess unknowns
c
c constants:
c
*      xepirou=0.85d0                              ! case 1
*      xlamda=0.1d0
*      vv=0.001d0
*      xalpha=0.6d0
*      xmu=0.01d0
c
c      xepirou=0.01d0   ! epsilon                  ! case 2 (J. Phys. A paper)
c      xlamda=0d0       ! gamma
c      vv=0.01d0        ! nu
c      xalpha=0.5d0
c      xmu=0.1d0
c
      PAR(1)=xepirou
      PAR(2)=xlamda
      PAR(3)=vv
      PAR(4)=xalpha
      PAR(5)=xmu
c
c shooting parameters:
c
c      xl=24d0                                     ! for case 2
c      delta=1d0
c
*       xl=23d0                                    ! for case 1
*       delta=5d0
c
       xl=x(1)
       delta=x(2)
c
*      x(1)=xl
*      x(2)=delta
c
      twopi=8d0*datan(1d0)
      if (newton) then                  ! skip for simple forward integration
      do 5 i=1,1                               ! loop for simple continuation
         xtol=1d-12
         maxfev=50
         ml=nx-1
         mu=nx-1
         epsfcn=1d-12
         mode=1
         factor=100d0
         nprint=0
         ldfjac=nx
         call hybrd(fcn,nx,x,fv,xtol,maxfev,ml,mu,epsfcn,diag,
     +              mode,factor,nprint,info,nfev,fjac,ldfjac,r,lnx,
     +              qtf,wa1,wa2,wa3,wa4)
         if (info.eq.1) then                                    ! convergence
            if (x(2).lt.0d0) x(2)=x(2)+twopi
            if (x(2).gt.twopi) x(2)=x(2)-twopi
            write(6,'("shooting converged")')
            write(6,'("xl: ",f12.8,"   delta:",f12.8)') (x(j),j=1,nx)
            write(10,*) xepirou,xlamda,vv,xalpha,xmu       ! write parameters
            write(10,*) fp(1),fp(2),fp(3),fp(4)
            write(10,*) x(1),x(2)
         endif
*         write(16,13) (x(j),j=1,nx)
5     continue 
      endif
c
c final integration:
c
      xa=0d0
      xend=1d0
      call initial(node,ys)
      xtol=1d-12
      iout=2
      itol=0
      rtol=xtol
      atol=xtol
      do 7 i=1,liws
         iws(i)=0
7     continue
      do 8 i=1,lws
         ws(i)=0d0
8     continue
      iws(3)=-1
      iws(4)=1
      iws(5)=nrd
      call dop853(node,f,xa,ys,xend,rtol,atol,itol,output,iout,ws,
     +            lws,iws,liws,rpar,ipar,idid)
      nr=1                                                      ! final point
      call output(nr,xold,xa,ys,node,con,icomp,nd,rpar,ipar,irtrn)
*      write(15,13)ys(1),ys(2),ys(3),ys(4)
*13    format(4f17.8)
      stop
      end
c
c****************************************************************************
c Set up the residual equations (to be zeroed in the shooting method)
c****************************************************************************
      subroutine fcn(n,x,fv,iflag)
      implicit real*8 (a-h,o-z)
c
      parameter (node=4,nrd=4)                           ! DOPRI (integrator)
      parameter (lws=11*node+8*nrd+21,liws=nrd+21)
      dimension ws(lws),iws(liws)
c
      dimension fv(n),x(n),ys(node)
      external f,output
      common /init/xl,delta,fp(4)
c
c enter the unknowns:
c
      xl=x(1)                                                 ! (half) length
      delta=x(2)                                                      ! delta
c
      call initial(node,ys)
      xa=0d0
      xend=1d0
      tol=1d-12
c
      iout=0
      itol=0
      rtol=tol
      atol=tol
      do 7 i=1,liws
         iws(i)=0
7     continue
      do 8 i=1,lws
         ws(i)=0d0
8     continue
      iws(3)=-1
      iws(4)=1
      iws(5)=nrd
      call dop853(node,f,xa,ys,xend,rtol,atol,itol,output,iout,ws,
     +            lws,iws,liws,rpar,ipar,idid)
c
c set up the residuals:
c
c --- R2 reversible solutions (shoot into symmetric section):
c
      fv(1)=dsin(ys(2)/2d0)                                 ! psi=0, 2pi, ...
      fv(2)=ys(3)                                                   ! psita=0
c
      write(*,13)x(1),x(2),fv(1),fv(2)
13    format(4e15.6)
c     write(*,*)' '
      return
      end
c
c****************************************************************************
c Subroutine for setting the initial conditions
c****************************************************************************
      subroutine initial(n,ys)
      implicit real*8 (a-h,o-z)
      parameter (eps=1d-5)
      dimension ys(n),a(n,n),v(2,4),x(4),PAR(20)
      common /init/xl,delta,fp(4)
c
c fixed point (for xepirou=0, xlamda=0, vv=0):
c
*      fp(1)=9.9668652492E-02
*      fp(2)=0d0
*      fp(3)=0d0
*      fp(4)=9.9503719021E-01
c
*      fp(1)=atan(sqrt(xmu)/xalpha)
*      fp(2)=0d0
*      fp(3)=0d0
*      fp(4)=xalpha/sqrt(xalpha*xalpha+xmu)
c
c find fixed point
c
      PAR(1)=xepirou
      PAR(2)=xlamda
      PAR(3)=vv
      PAR(4)=xalpha
      PAR(5)=xmu
c
       call solvefix(x,fv2,PAR)
       fp=x
*       print*,'testp',fp(1),fp(2),fp(3),fp(4)
c
c find eigenvectors:
c
      call xjac(n,fp,a)
      call eigs(a,v)
c
      do 50 i=1,4
         ys(i)=fp(i)+eps*(v(1,i)*dcos(delta)+v(2,i)*dsin(delta))
50    continue
c
      return
      end
c
c****************************************************************************
c Subroutine for evaluating the vector field
c****************************************************************************
      subroutine f(n,r,u,du,rpar,ipar)
      implicit real*8 (a-h,o-z)
      dimension u(n),du(n)
      common /init/xl,delta,fp(4)
      common /pars/xepirou,xlamda,vv,xalpha,xmu
c
c equilibrium equation:
c
      du(1)= u(3)                                                    ! theta
c
      du(2)= (u(4)-dcos(u(1)))/dsin(u(1))/dsin(u(1))-                ! psi
     +        vv*(1+xlamda*dcos(u(1)))*
     +      dsin(u(1))*dcos(u(2))/dsqrt(xmu-2*vv*u(4))
     +      -xlamda*vv/xalpha*dsin(u(1))*dcos(u(2))*
     +       dsin(u(1))*dcos(u(2))
     +      -xepirou
c
      du(3)= -(u(4)-dcos(u(1)))*(1-u(4)*dcos(u(1)))/dsin(u(1))       ! p_theta
     +         /dsin(u(1))/dsin(u(1))
     +      +xalpha*dsin(u(1))*(1+xlamda*dcos(u(1)))
     +        -(dcos(u(1))+xlamda*(dcos(u(1))*dcos(u(1))
     +         -dsin(u(1))*dsin(u(1))))*
     +       dcos(u(2))*dsqrt(xmu-2*vv*u(4))
     +       -xlamda/xalpha*dsin(u(1))*dcos(u(1))*dcos(u(2))*
     +       dcos(u(2))*(xmu-2*vv*u(4))
c
      du(4)= (1+xlamda*dcos(u(1)))*dsin(u(1))*dsin(u(2))*            ! p_psi
     +         dsqrt(xmu-2*vv*u(4))+
     +        xlamda/xalpha*dsin(u(2))*(xmu-2*vv*u(4))*
     +         dsin(u(1))*dsin(u(1))*dcos(u(2))
c
c rescale the equations:
c
      do 100 i=1,n
         du(i)=xl*du(i)
100   continue
c
      return
      end
c
c****************************************************************************
c Subroutine for evaluating the vector field
c****************************************************************************
      subroutine xjac(n,u,a)
      implicit real*8 (a-h,o-z)
      dimension u(n),a(n,n)
      common /init/xl,delta,fp(4)
      common /pars/xepirou,xlamda,vv,xalpha,xmu
c
c the variational equations about a fixed point:
c
*        a(1,1)=0
*        a(1,2)=0
*        a(1,3)=1
*        a(1,4)=0
c
*        a(2,1)=sqrt((xalpha*xalpha+xmu)/xmu)
*        a(2,2)=0
*        a(2,3)=0
*        a(2,4)=(xalpha*xalpha+xmu)/xmu
c
*        a(3,1)=-1+sqrt(xalpha*xalpha+xmu)
*        a(3,2)=0
*        a(3,3)=0
*        a(3,4)=-sqrt((xalpha*xalpha+xmu)/xmu)
c
*        a(4,1)=0
*        a(4,2)=xmu/sqrt(xalpha*xalpha+xmu)
*        a(4,3)=0
*        a(4,4)=0
c
        sita=u(1)
        psi=u(2)
        Ppsi=u(4)
c
*        print*,'test',sita,psi,ppsi
c
        Ss1=sin(sita)
        Sc1=cos(sita)
        Ss2=sin(psi)
        Sc2=cos(psi)
        Sp=Ppsi
c
        a(1,1)=0
        a(1,2)=0
        a(1,3)=1
        a(1,4)=0
c
        a(2,1)=1/Ss1-2*Sc1*(Sp-Sc1)/Ss1/Ss1/Ss1
     +        -(vv*Sc1*Sc2*(1+xlamda*Sc1)
     +          -xlamda*vv*Ss1*Ss1*Sc2)/
     +         sqrt(xmu-2*vv*Sp)
     +        -2*xlamda*vv*Ss1*Sc1*Sc2*Sc2/xalpha
c
        a(2,2)=vv*(1+xlamda*Sc1)*Ss1*Ss2/sqrt(xmu-2*vv*Sp)
     +         +2*xlamda*vv*Ss1*Ss1*Ss2*Sc2/xalpha
c
        a(2,3)=0
c
        a(2,4)=1/Ss1/Ss1-vv*vv*(1+xlamda*Sc1)*Ss1*Sc2/
     +        sqrt((xmu-2*vv*Sp)*(xmu-2*vv*Sp)*(xmu-2*vv*Sp))
c
        a(3,1)=3*Sc1*(Sp-Sc1)*(1-Sp*Sc1)/Ss1/Ss1/Ss1/Ss1
     +         -(1-2*Sp*Sc1+Sp*Sp)/Ss1/Ss1+xalpha*Sc1*(1+xlamda*Sc1)
     +         -xalpha*xlamda*Ss1*Ss1+(Ss1+4*xlamda*Ss1*Sc1)*Sc2*
     +         sqrt(xmu-2*vv*Sp)+xlamda*(Ss1*Ss1-Sc1*Sc1)*Sc2*Sc2*
     +         (xmu-2*vv*Sp)/xalpha
c
        a(3,2)=Ss2*(Sc1+xlamda*(Sc1*Sc1-Ss1*Ss1))*sqrt(xmu-2*vv*Sp)
     +        +2*xlamda*Ss2*Sc2*Ss1*Sc1*(xmu-2*vv*Sp)/xalpha
c
        a(3,3)=0
c
        a(3,4)=-a(2,1)
c
        a(4,1)=a(3,2)
c
        a(4,2)=(1+xlamda*Sc1)*Ss1*Sc2*sqrt(xmu-2*vv*Sp)
     +        +xlamda*Ss1*Ss1*(Sc2*Sc2-Ss2*Ss2)*(xmu-2*vv*Sp)/
     +         xalpha
c
        a(4,3)=0
c
        a(4,4)=-a(2,2)
c
      return
      end
c
c****************************************************************************
c Subroutine for finding the fix point
c****************************************************************************
      subroutine solvefix(x,fv2,PAR)
      implicit real*8 (a-h,o-z)
c
      parameter (nx=4)
      parameter (lnx=(nx*(nx+1))/2)
      dimension diag(nx),fjac(nx,nx),r(lnx),qtf(nx)
      dimension wa1(nx),wa2(nx),wa3(nx),wa4(nx)
      dimension x(nx),fv(nx), PAR(20), PAR1(20)
      external fcnn
c
c  initial guesses:
c
      x(1)=0.5d0
      x(2)=0d0
      x(3)=0d0
      x(4)=1d0
c
      xtol=1d-10                                                   ! tolerance
      maxfev=200                                ! maximum number of iterations
      ml=nx-1
      mu=nx-1
      epsfcn=1d-12
      mode=1
      factor=100d0
      nprint=0
      ldfjac=nx
      call hybrd(fcnn,nx,x,fv,xtol,maxfev,ml,mu,epsfcn,diag,
     +           mode,factor,nprint,info,nfev,fjac,ldfjac,r,lnx,
     +           qtf,wa1,wa2,wa3,wa4)
      if (info.eq.1) then                                        ! convergence
*       write(6,'(/"iterations converged")')
*        write(6,'("u1: ",f16.12,"   u2:",f16.12,"   u3: ",f16.12,
*    +"   u4:",f16.12)') (x(j),j=1,nx)
      else
         write(6,'(/" no convergence -- error:",i3/)') info
      endif

c      stop
      end
c
c*****************************************************************************
c Set up the residual equations (to be zeroed in the Newton iterations)
c*****************************************************************************
      subroutine fcnn(nx,x,fv2,iflag)
      implicit real*8 (a-h,o-z)
c
      dimension fv2(nx),x(nx), PAR(20)
      common /pars/xepirou,xlamda,vv,xalpha,xmu
c
      U1=x(1)
      U2=x(2)
      U3=x(3)
      U4=x(4)
c
      eq1= U3
c
      eq2= (U4-dcos(U1))/dsin(U1)/dsin(U1)-
     +        vv*(1d0+xlamda*dcos(U1))*
     +      dsin(U1)*dcos(U2)/dsqrt(xmu-2*vv*U4)
     +      -xlamda*vv/xalpha*dsin(U1)*dcos(U2)*
     +       dsin(U1)*dcos(U2)
     +      -xepirou           
c
      eq3= -(U4-dcos(U1))*(1d0-U4*dcos(U1))/dsin(U1)
     +         /dsin(U1)/dsin(U1)
     +      +xalpha*dsin(U1)*(1d0+xlamda*dcos(U1))
     +        -(dcos(U1)+xlamda*(dcos(U1)*dcos(U1)
     +         -dsin(U1)*dsin(U1)))*
     +       dcos(U2)*dsqrt(xmu-2*vv*U4)
     +       -xlamda/xalpha*dsin(U1)*dcos(U1)*dcos(U2)*
     +       dcos(U2)*(xmu-2*vv*U4)
c
      eq4= (1d0+xlamda*dcos(U1))*dsin(U1)*dsin(U2)*
     +         dsqrt(xmu-2*vv*U4)+
     +        xlamda/xalpha*dsin(U2)*(xmu-2*vv*U4)*
     +         dsin(U1)*dsin(U1)*dcos(U2)
c
      fv2(1)=eq1
      fv2(2)=eq2
      fv2(3)=eq3
      fv2(4)=eq4
c
*      write(*,13) x(1),x(2),x(3),x(4),fv2(1),fv2(2),fv2(3),fv2(4)
c
*13    format(4e15.6,"   ",4e15.6)
c     write(*,*)' '
      return
      end
c
c       --------------------
        subroutine eigs(a,v)
c       --------------------
c
c       -calculates eigenvectors spanning the unstable manifold 
c        (using the NAG routine F02AGF).
c
	implicit double precision (a-h,o-z)
	parameter(ia=4,n=4,ivr=4,ivi=4,iou=4,zero=1d-7)
        dimension a(4,4),rr(4),ri(4),vi(4,4),vr(4,4),v(2,4)
        dimension rrdum(4),ridum(4),vrdum(4,4),vidum(4,4)
	dimension intger(4)
c
        ifail=0
c
	call f02agf(a,ia,n,rr,ri,vr,ivr,vi,ivi,intger,ifail)
c
        if (ifail.ne.0) then   
           print*,'failed  !!!'
           print*,'ifail=',ifail
        endif
*	print*,'number of iterations for each eigenvalue:'
*	print 10,intger(1),intger(2),intger(3),intger(4)
10	format(2x,i10,2x,i10,2x,i10,2x,i10,2x,i10,2x,i10)
c
c  ORDER THE EIGENVECTORS/VALUES ACCORDING SIZE OF REAL PART OF EIGENVALUE.
c
        do 200 i=1,3
           do 100 j=i+1,4
              if (rr(i).gt.rr(j)) then
                 rrdum(i)=rr(i)
                 ridum(i)=ri(i)
                 rr(i)=rr(j)
                 ri(i)=ri(j)
                 rr(j)=rrdum(i)
                 ri(j)=ridum(i)
                 do 70 k=1,4
                    vrdum(k,i)=vr(k,i)
                    vidum(k,i)=vi(k,i)
                    vr(k,i)=vr(k,j)
                    vi(k,i)=vi(k,j)
                    vr(k,j)=vrdum(k,i)
                    vi(k,j)=vidum(k,i)
70               continue
              endif
100        continue
200     continue
c
c  USE REAL AND IMAGINARY PARTS OF COMPLEX EIGENVECTORS.
c
        do 300 i=1,3
           if (dabs(vr(1,i)-vr(1,i+1)).lt.zero) then
              if (dabs(vr(2,i)-vr(2,i+1)).lt.zero) then
                 if (dabs(vr(3,i)-vr(3,i+1)).lt.zero) then
                    if (dabs(vr(4,i)-vr(4,i+1)).lt.zero) then
                             do 250 j=1,4
                                vr(j,i+1)=vi(j,i)
250                          continue
                          endif
                       endif
                    endif
                 endif
300     continue

        check_re1=a(1,1)*vr(3,1)+a(1,2)*vr(3,2)+a(1,3)*vr(3,3)
     +        +a(1,4)*vr(3,4)
     +        -rr(3)*vr(3,1)+ri(3)*vr(4,1)
        check_im1=a(1,1)*vr(4,1)+a(1,2)*vr(4,2)+a(1,3)*vr(4,3)
     +        +a(1,4)*vr(4,4)
     +        -rr(3)*vr(4,1)-ri(3)*vr(3,1)
        check_re2=a(1,1)*vr(4,1)+a(1,2)*vr(4,2)+a(1,3)*vr(4,3)
     +        +a(1,4)*vr(4,4)
     +        -rr(4)*vr(4,1)+ri(4)*vr(3,1)
        check_im2=a(1,1)*vr(3,1)+a(1,2)*vr(3,2)+a(1,3)*vr(3,3)
     +        +a(1,4)*vr(3,4)
     +        -rr(4)*vr(3,1)-ri(4)*vr(4,1)
*        write(6,*) 'check:',check_re1,check_im1,check_re2,
*     +             check_im2                             ! check eigenvectors
c
	write(iou,*) '(ordered) eigenvalues: '
	write(iou,*) '   (',rr(1),ri(1),')'
	write(iou,*) '   (',rr(2),ri(2),')'
	write(iou,*) '   (',rr(3),ri(3),')'
	write(iou,*) '   (',rr(4),ri(4),')'
	write(iou,*) 'corresponding eigenvectors: '
	write(iou,500) vr(1,1),vr(2,1),vr(3,1),vr(4,1)
	write(iou,500) vr(1,2),vr(2,2),vr(3,2),vr(4,2)
	write(iou,500) vr(1,3),vr(2,3),vr(3,3),vr(4,3)
	write(iou,500) vr(1,4),vr(2,4),vr(3,4),vr(4,4)
c
c  LET v CONTAIN THE TWO APPROPRIATE 4D EIGENVECTORS.
c
	do 350 i=1,2
	   do 330 j=1,4
	      v(i,j)=0.0
330	   continue
350	continue
	do 400 i=1,4
	   do 390 j=1,2
	      v(j,i)=vr(i,j+2)
390	   continue
400	continue

*	print*,v(1,1),v(1,2),v(1,3),v(1,4)
*	print*,v(2,1),v(2,2),v(2,3),v(2,4)
c
        return
c
500	format(2x,'(',f16.12,',',1x,f16.12,',',1x,f16.12,',',1x,
     +   f16.12,')')
c
        end
c
c****************************************************************************
c Subroutine for writing the output
c****************************************************************************
      subroutine output(nr,xold,x,y,n,con,icomp,nd,rpar,ipar,irtrn)
      implicit real*8 (a-h,o-z)
      parameter(np=4,xplot=2d-3)
      dimension y(n),z(np),con(8*nd),icomp(nd)
      common /init/xl,delta,fp(4)
      common/intern/xout
c
      if (nr.eq.1) then
         write(12,100) xl*x,(y(i),i=1,n)
         xout=x+xplot
      else
10       continue
         if (x.ge.xout) then
c
            do 50 i=1,n
               z(i)=contd8(i,xout,con,icomp,nd)
50          continue
c
            write(12,100) xl*xout,(z(i),i=1,n)
c            print*,z(1),z(2),z(3)
            xout=xout+xplot
            goto 10
         endif
      endif
c
      return
100   format(1x,7(e16.8,1x))
      end
