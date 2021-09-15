      PROGRAMfree
      implicit double precision(a-h,o-z)
c      LOGICAL CHECK
c      logical sign,sign1
      REAL*8 xxx(1),FVEC(1),rmup,rmun,z(1)
      real*8 X0(1),X1(1),XMIN(1),XMAX(1),fout(1)
      REAL*8 rhoinf,rhosup, a1, a2,pi,pi2,hc,drho,rho
      real*8 yp,dens,rhop,rhon,rmmu,rm
      real*8 rlambm,dmupdrhop,dmupdrhon,drndrp,rnrp,qui
      integer npoint,i,j,ipar
      common/temp/pi,pi2,gamma,rho0
c      common/coup/rmv,rms,rmrho,gs,gv,grho,rka,rlambda,xsi,gwr
      COMMON/PDN/X0,X1,XMIN,XMAX,stepmax,acc
      COMMON/PDNa/ X0a,X1a,XMINa,XMAXa
      common/fields/rv0,b0
      COMMON/MA/rme
      COMMON/ME/rmelect,xkfe
      common/mu/rmmu,xkfmu
      COMMON/ROS/rhop,rhon,xkfp,xkfn,rmun,rmup
      common/vec/fvec
      common/ro/rhob
      external thermo,f3p,f3n,thermo1,y,yd
      pi=dacos(-1.d0)
      pi2=pi*pi
      hc=197.326d0
      RM=939.D0
      RMELECT=0.511d0/rm
      rmmu=108.d0/rm
      gamma=2
      pi=dacos(-1.d0)
      pi2=pi*pi
      hc=197.326d0
      RM=939.D0
      rho0=0.148d0/((RM/hc)**3)		
        acc=.1d-11  
        x0(1)=-0.89D0
        x1(1)=-0.87D0
        xmin(1)=-1.D0
	XMAX(1)=0.0d0
        acc=.1d-11  
        stepmax=0.02d0
        ierr=0
       do i=1,150
c          xlog=-1.d0+0.01d0*(i-1)
c          rho=10**xlog
          rho=0.000001+0.01*i
          rhob=rho*(hc/rm)**3          
c----------------------------------------------------------
c     calculo de derivadas
c----------------------------------------------------------
      call thermo(delta,rme,esy,press,dener,ener)
      rhoe=xkfe**3/(3.d0*pi2)
      rhomu=xkfmu**3/(3.d0*pi2)
      write(6,*)'rho=',rho,delta
      write(70,99)rhob*(rm/hc)**3,ener*rm*(rm/hc)**3,press*rm*(rm/hc)**3
     1    ,esy*rm,rhon/rhob,rhop/rhob,
     2     rhoe/rhob,rhomu/rhob
      write(71,99)rhob*(rm/hc)**3,ener*(rm/hc)**4,press*(rm/hc)**4
      enddo
C-----------------------------------------------------------------
 99   format(2x,30(d12.4,3x))
 20   format(2x,10d13.4)
 32   format(2x,10d13.4)
      stop
      end      
c#######################################################################
c     SUBROTINAS
C#######################################################################
      subroutine thermo(delta,rma,esy,press,dener,ener)
      implicit none
      real*8 mup,mun,mue,nup,nun,nue,rho,v0,rhot,pi,rhoa
      real*8 delta,dens,press,rmup,rmun,rmue,rnup,rnun,rnue,T
      real*8 PI2,GS,GV,GRHO,RKA,RLAMBDA,XSI,RMA,RINFP,RINFN,RME
      real*8 RHOB,RHOP,RHON,RHOE,RHO3,rhoin,del1
      real*8 acc,stepmax,RPHI,B0,RV0,RM,RMS,RMV,RMRHO,y,yd
      real*8 X0(1),X1(1),XMIN(1),XMAX(1),x(1)
      real*8 xkfp,xkfn,xkfe,gwr,dener,ener,esy,rmrho2,rho0
      real*8 yp,ef,xkf,res4n,res4p,res4e,res3n,res3p,res3e
      real*8 f4p,f4n,f4e,f3p,f3n,f3e,gamma,rmelect,rmmu
      real*8 res3mu,res4mu,f3mu,f4mu,xkfmu
      integer ierr,iact,ii,i      

      COMMON/ROS/RHOP,RHON,xkfp,xkfn,rmun,rmup
      common/temp/pi,pi2,gamma,rho0 
      COMMON/PDN/X0,X1,XMIN,XMAX,stepmax,acc
      common/fields/rv0,b0
      COMMON/MA/RMe
      common/mu/rmmu,xkfmu
      common/ro/dens
      COMMON/ME/rmelect,xkfe
      external dnewto,y,yd,f4p,f4n,f4e,f3p,f3n,f3e,f3mu,f4mu
c-----------------------------------------------------------------------
C     DADOS PARA DNEWTON
C-----------------------------------------------------------------------
      rhob=dens
c      write(6,*)'2 rho=', rhob
      ierr=0
      acc=1.d-10
      call dnewto(x0,x1,xmin,xmax,1,y,yD,acc,
     1     stepmax,90,x,iact,ierr)
      write(6,*)'x=',x
      rma=1.d0
      delta=x(1)
      rmue=rmun-rmup
      rho3=delta*rhob
      RHO3=RHOP-RHON
c----------------------------------------------------------------------
C     DEFININDO ENERGIA, PRESSAO, ENTROPIA
C----------------------------------------------------------------------
      CALL GAUSS(F3P,0.D0,xkfP,10,RES3P,II)
      CALL GAUSS(F3N,0.D0,xkfN,10,RES3N,II)
      CALL GAUSS(F3E,0.D0,xkfE,10,RES3E,II)
      CALL GAUSS(F3mu,0.D0,xkfmu,10,RES3mu,II)
      ENER=GAMMA/(2.D0*PI2)*(RES3P+RES3N+RES3E+res3mu) 
      DENER=ENER/RHOB-1.D0
c-------------------------------------------------------------------
C     DEFININDO PRESSAO
C------------------------------------------------------------------
      CALL GAUSS(F4P,0.D0,xKfP,10,RES4P,II)
      CALL GAUSS(F4N,0.D0,xKfN,10,RES4N,II)
      CALL GAUSS(F4E,0.D0,xKfE,10,RES4E,II)
      CALL GAUSS(F4mu,0.D0,xKfmu,10,RES4mu,II)
      PRESS=GAMMA/(6.D0*PI2)*(RES4P+RES4N+RES4E+res4mu)
      xkf=(3.d0*pi2*rhob/2.d0)**(1.d0/3.d0)
      ef=dsqrt(xkf**2+rme**2)
      esy=xkf**2/(6.d0*ef)
      return
      end
c-------------------------------------------------------
      function y(i,z,idim,ierr)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Z(1)
      external F1P,F1N,ya,yda
      COMMON/MA/RMA
      COMMON/ROS/rhop,rhon,xkfp,xkfn,rmun,rmup
      common/temp/pi,pi2,gamma,rho0 
      common/ro/rhob
      COMMON/ME/rmelect,xkfe
      common/fields/rv0,b0
      common/mu/rmmu,xkfmu

      delta=z(1)
      rho3=delta*rhob
      rhop=(rhob+rho3)/2.d0
      rhon=(rhob-rho3)/2.d0
      xkfp=(3.D0*PI2*rhop)**(1.D0/3.D0)
      xkfn=(3.D0*PI2*rhon)**(1.D0/3.D0)
      rmun=(xkfn**2.D0+rma**2.D0)**(1.D0/2.D0)
      rmup=(xkfp**2.D0+rma**2.D0)**(1.D0/2.D0)
      rmue=rmun-rmup
      xkfe2=(rmue**2.D0-rmelect**2.D0)
      xkfmu2=(rmue**2.D0-rmmu**2.D0)
      if(xkfmu2.le.0.d0)then
         xkfmu=0.d0
      else
         xkfmu=dsqrt(xkfmu2)
      endif
      if(xkfe2.le.0.d0)then
         xkfe=0.d0
      else
         xkfe=dsqrt(xkfe2)
      endif
      IF(I.EQ.1)Y=xkfe**3+xkfmu**3-xkfp**3
      return
      end
c-------------------------------------------------------------------
      function yD(i,j,x,idim,ierr)
      implicit real *8 (A-H,O-Z) 
      yd=0.D0
      return
      end
c-------------------------------------------------------------------
      FUNCTION F1P(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/MA/RMA
      E=DSQRT(X*X+RMA*RMA)
      F1P=X*X/E
      RETURN
      END
C------------------------------------------------------------
      FUNCTION F1N(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/MA/RMA
      E=DSQRT(X*X+RMA*RMA)
      F1N=X*X/E
      RETURN
      END
C------------------------------------------------------------
      FUNCTION F2(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/MA/RMA
      E=DSQRT(X*X+RMA*RMA)
      F2=X*X
      RETURN
      END
C----------------------------------------------------------------
      FUNCTION F3P(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/MA/RMA
      RME=RMA
      E=DSQRT(X*X+RME*RME)
      F3P=X*X*E
      RETURN
      END
C-----------------------------------------------------------------
      FUNCTION F3N(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/MA/RMA
      RME=RMA
      E=DSQRT(X*X+RME*RME)
      F3N=X*X*E
      RETURN
      END
C-----------------------------------------------------------------
      FUNCTION F3E(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ME/RMELECT,xkfe
      RME=RMELECT
      E=DSQRT(X*X+RME*RME)
      F3E=X*X*E
      RETURN
      END
C-----------------------------------------------------------------
      FUNCTION F4mu(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      common/mu/rmmu,xkfmu
      RME=RMmu
      E=DSQRT(X*X+RME*RME)
      F4mu=X**4/E
      RETURN
      END
C-----------------------------------------------------------------
      FUNCTION F3mu(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      common/mu/rmmu,xkfmu
      RME=RMmu
      E=DSQRT(X*X+RME*RME)
      F3mu=X*X*E
      RETURN
      END

C------------------------------------------------------------------
      FUNCTION F4P(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/MA/RMA
      RME=RMA
      E=DSQRT(X*X+RME*RME)
      F4P=X**4/E
      RETURN
      END
c------------------------------------------------------------------
      FUNCTION F4N(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/MA/RMA
      RME=RMA
      E=DSQRT(X*X+RME*RME)
      F4N=X**4/E
      RETURN
      END
c------------------------------------------------------------------
      FUNCTION F4E(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ME/RMELECT,xkfe
      RME=RMELECT
      E=DSQRT(X*X+RME*RME)
      F4E=X**4/E
      RETURN
      END
C---------------------------------------------------------------
        SUBROUTINE GAUSS(F,UG,OG,NN,FXINT,II)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION WG(10),ZG(10)
        DATA WG/0.6667134430869 D-01,
     C  0.1494513491506D00,
     C  0.2190863625160D00,
     C  0.2692667193100D00,
     C  0.2955242247148D00,
     C  0.2955242247148D00,
     C  0.2692667193100D00,
     C  0.2190863625160D00,
     C  0.1494513491506D00,
     C           0.6667134430869D-01/
        DATA ZG/-0.9739065285172D00,
     C  -0.8650633666890D00,
     C  -0.6794095682990D00,
     C  -0.4333953941292D00,
     C  -0.1488743389816D00,
     C  +0.1488743389816D00,
     C  +0.4333953941292D00,
     C  +0.6794095682990D00,
     C  +0.8650633666890D00,
     C          +0.9739065285172D00/

        FXINT=0.D0
        HH=(OG-UG)/DBLE(FLOAT(NN))
        U=UG
        O=U+HH
        KK=1
24      OU=O+U
        RI=0.D0
        DO 26 I=1,10
        X=0.5D0*(ZG(I)*HH+OU)
        FUNCAO=F(X)
26      RI=RI+WG(I)*FUNCAO
        FXINT=RI*HH/2.D0+FXINT
        KK=KK+1
        IF(KK-NN)28,28,9999
28      U=O
        O=O+HH
        GO TO 24
9999    RETURN
        END
c---------------------------------------------------------------
        include 'dnewton.f'
c        include 'broydn.f'



