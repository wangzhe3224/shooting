c****************************************************************************
c Shooting code based on HYBRD for homoclinic orbits of the dimensionless
c equations for an initially curved rod. Code includes fragments of
c 'curvenew.f' (e.g., the subroutine 'eigs') and is initialised to give 
c the same solution as 'curvenew.f', namely the solution of the Phil. Trans.
c paper.
c
c The order of the quantities in the integration is:
c
c     y(1)=F1    y(2)=F2    y(3)=F3    y(4)=M1    y(5)=M2    y(6)=M3
c
c The order of the unknown quantities in the shooting method is:
c
c           x(1)=xl     x(2)=delta
c
c     January 2006 (Gert van der Heijden)
c
c****************************************************************************
      program shoot
      implicit real*8 (a-h,o-z)
c
      parameter (node=6,nrd=6)                           ! DOPRI (integrator)
      parameter (lws=11*node+8*nrd+21,liws=nrd+21)
      dimension ws(lws),iws(liws)
c
      parameter (nx=2,lnx=(nx*(nx+1))/2)              ! HYBRD (Newton solver)
      dimension diag(nx),fjac(nx,nx),r(lnx),qtf(nx)
      dimension wa1(nx),wa2(nx),wa3(nx),wa4(nx)
c
      dimension x(nx),fv(nx),ys(node)
      common /init/xl,delta
      common /pars/enu,rho,em,ekappa
      external fcn,f,output
      logical newton
c
      newton=.true.                                      ! Newton iterations?
c
*      open(unit=1,file='shoot.dat',status='unknown')
*      read(1,*)enu,rho,em                                       ! parameters
*      read(1,*)(x(j),j=1,nx)                                ! guess unknowns
c
c constants:
c
      enu=1d0/3d0
      rho=0d0
      ekappa=0.02d0
      em=1.7d0
c
c shooting parameters:
c
c          primary solution:   (xl,delta) = (41.29942737, 5.15691739)
c          bi-modal solution:  (xl,delta) = (52.12980358, 5.52777736)
c
      xl=41.30d0
      delta=5.16d0
c
      x(1)=xl
      x(2)=delta
c
      if (newton) then                  ! skip for simple forward integration
      do 5 i=1,1                               ! loop for simple continuation
         xtol=1d-12
         maxfev=20
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
            write(6,'("shooting converged")')
            write(6,'("xl: ",f12.8,"   delta:",f12.8)') (x(j),j=1,nx)
*            write(10,*) eps,enu,rho,em,ekappa             ! write parameters
*            write(10,*) xl,delta
         endif
         write(16,13) (x(j),j=1,nx)
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
*      write(15,13)ys(1),ys(2),ys(3),ys(4),ys(5),ys(6)
13    format(6f17.8)
      stop
      end
c
c****************************************************************************
c Set up the residual equations (to be zeroed in the shooting method)
c****************************************************************************
      subroutine fcn(n,x,fv,iflag)
      implicit real*8 (a-h,o-z)
c
      parameter (node=6,nrd=6)                           ! DOPRI (integrator)
      parameter (lws=11*node+8*nrd+21,liws=nrd+21)
      dimension ws(lws),iws(liws)
c
      dimension fv(n),x(n),ys(node)
      external f,output
      common /init/xl,delta
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
      fv(1)=ys(2)                                                      ! F2=0
      fv(2)=ys(5)                                                      ! k2=0
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
      dimension ys(n),fp(n),a(n,n),v(3,6)
      common /init/xl,delta
      common /pars/enu,rho,em,ekappa
c
c fixed point (for enu=1/3, rho=0, em=1.7, ekappa=0.02):
c
      fp(1)=0.0337259272127d0
      fp(2)=0d0
      fp(3)=0.999431119104d0
      fp(4)=0.0249794593288d0
      fp(5)=0d0
      fp(6)=0.999687964623d0
c
c find eigenvectors:
c
      call xjac(n,fp,a)
      call eigs(a,v)
c
      do 50 i=1,6
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
      common /init/xl,delta
      common /pars/enu,rho,em,ekappa
c
      eP=1/em**2
      epsilon=0d0
      sigma=0d0
c
c force and moment balance:
c
      du(1)= (1+enu)*(1+rho/2)*u(2)*u(6) - u(3)*u(5)
      du(2)= (1+rho)*u(3)*(u(4)+ekappa) - (1+enu)*(1+rho/2)*u(1)*u(6)
      du(3)= u(1)*u(5) - (1+rho)*u(2)*(u(4)+ekappa)
      du(4)= ((1+enu)*(1+rho/2)-1)*u(5)*u(6) + eP*(u(2)+epsilon*
     +       gamma*u(2)*u(3))
      du(5)= (1+rho)*(u(4)+ekappa)*u(6) - (1+enu)*(1+rho/2)*u(4)*u(6)
     +       + eP*(epsilon*(sigma-gamma)*u(1)*u(3)-u(1))
      du(6)= u(4)*u(5) - (1+rho)*(u(4)+ekappa)*u(5) - eP*epsilon*
     +       sigma*u(1)*u(2)
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
      common /init/xl,delta
      common /pars/enu,rho,em,ekappa
c
        eP=1/em**2
        epsilon=0d0
        sigma=0d0
c
c the variational equations about a fixed point:
c
        a(1,1)=0
        a(1,2)=(1+enu)*(1+rho/2)*u(6)
        a(1,3)=-u(5)
        a(1,4)=0
        a(1,5)=-u(3)
        a(1,6)=(1+enu)*(1+rho/2)*u(2)
c
        a(2,1)=-(1+enu)*(1+rho/2)*u(6)
        a(2,2)=0
        a(2,3)=(1+rho)*(u(4)+ekappa)
        a(2,4)=(1+rho)*u(3)
        a(2,5)=0
        a(2,6)=-(1+enu)*(1+rho/2)*u(1)
c
        a(3,1)=u(5)
        a(3,2)=-(1+rho)*(u(4)+ekappa)
        a(3,3)=0
        a(3,4)=-(1+rho)*u(2)
        a(3,5)=u(1)
        a(3,6)=0
c
        a(4,1)=0
        a(4,2)=eP*(1+epsilon*gamma*u(3))
        a(4,3)=eP*epsilon*gamma*u(2)
        a(4,4)=0
        a(4,5)=((1+enu)*(1+rho/2)-1)*u(6)
        a(4,6)=((1+enu)*(1+rho/2)-1)*u(5)
c
        a(5,1)=eP*(epsilon*(sigma-gamma)*u(3)-1)
        a(5,2)=0
        a(5,3)=eP*epsilon*(sigma-gamma)*u(1)
        a(5,4)=(rho-((1+enu)*(1+rho/2)-1))*u(6)
        a(5,5)=0
        a(5,6)=(1+rho)*(u(4)+ekappa)-(1+enu)*(1+rho/2)*u(4)
c
        a(6,1)=-eP*epsilon*sigma*u(2)
        a(6,2)=-eP*epsilon*sigma*u(1)
        a(6,3)=0
        a(6,4)=-rho*u(5)
        a(6,5)=u(4)-(1+rho)*(u(4)+ekappa)
        a(6,6)=0
c
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
	parameter(ia=6,n=6,ivr=6,ivi=6,iou=6,zero=1d-7)
        dimension a(6,6),rr(6),ri(6),vi(6,6),vr(6,6),v(3,6)
        dimension rrdum(6),ridum(6),vrdum(6,6),vidum(6,6)
	dimension intger(6)
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
*	print 10,intger(1),intger(2),intger(3),intger(4),intger(5),
*     +           intger(6)
10	format(2x,i10,2x,i10,2x,i10,2x,i10,2x,i10,2x,i10)
c
c  ORDER THE EIGENVECTORS/VALUES ACCORDING SIZE OF REAL PART OF EIGENVALUE.
c
        do 200 i=1,5
           do 100 j=i+1,6
              if (rr(i).gt.rr(j)) then
                 rrdum(i)=rr(i)
                 ridum(i)=ri(i)
                 rr(i)=rr(j)
                 ri(i)=ri(j)
                 rr(j)=rrdum(i)
                 ri(j)=ridum(i)
                 do 70 k=1,6
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
        do 300 i=1,5
           if (dabs(vr(1,i)-vr(1,i+1)).lt.zero) then
              if (dabs(vr(2,i)-vr(2,i+1)).lt.zero) then
                 if (dabs(vr(3,i)-vr(3,i+1)).lt.zero) then
                    if (dabs(vr(4,i)-vr(4,i+1)).lt.zero) then
                       if (dabs(vr(5,i)-vr(5,i+1)).lt.zero) then
                          if (dabs(vr(6,i)-vr(6,i+1)).lt.zero) then
                             do 250 j=1,6
                                vr(j,i+1)=vi(j,i)
250                          continue
                          endif
                       endif
                    endif
                 endif
              endif
           endif
300     continue
c
*	write(iou,*) '(ordered) eigenvalues: '
*	write(iou,*) '   (',rr(1),ri(1),')'
*	write(iou,*) '   (',rr(2),ri(2),')'
*	write(iou,*) '   (',rr(3),ri(3),')'
*	write(iou,*) '   (',rr(4),ri(4),')'
*	write(iou,*) '   (',rr(5),ri(5),')'
*	write(iou,*) '   (',rr(6),ri(6),')'
*	write(iou,*) 'corresponding eigenvectors: '
*	write(iou,500) vr(1,1),vr(2,1),vr(3,1),vr(4,1),vr(5,1),vr(6,1)
*	write(iou,500) vr(1,2),vr(2,2),vr(3,2),vr(4,2),vr(5,2),vr(6,2)
*	write(iou,500) vr(1,3),vr(2,3),vr(3,3),vr(4,3),vr(5,3),vr(6,3)
*	write(iou,500) vr(1,4),vr(2,4),vr(3,4),vr(4,4),vr(5,4),vr(6,4)
*	write(iou,500) vr(1,5),vr(2,5),vr(3,5),vr(4,5),vr(5,5),vr(6,5)
*	write(iou,500) vr(1,6),vr(2,6),vr(3,6),vr(4,6),vr(5,6),vr(6,6)
c
c  LET v CONTAIN THE THREE APPROPRIATE 6D EIGENVECTORS.
c
	do 350 i=1,3
	   do 330 j=1,6
	      v(i,j)=0.0
330	   continue
350	continue
	do 400 i=1,6
	   do 390 j=1,2
	      v(j,i)=vr(i,j+4)
390	   continue
400	continue
*	v(3,3)=1.0
*	print*,v(1,1),v(1,2),v(1,3),v(1,4),v(1,5),v(1,6)
*	print*,v(2,1),v(2,2),v(2,3),v(2,4),v(2,5),v(2,6)
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
      parameter(np=6,xplot=2d-3)
      dimension y(n),z(np),con(8*nd),icomp(nd)
      common /init/xl,delta
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
            xout=xout+xplot
            goto 10
         endif
      endif
c
      return
100   format(1x,7(e16.8,1x))
      end
