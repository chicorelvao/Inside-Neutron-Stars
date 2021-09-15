c-----------------------------------------------------------
c     em Dropbox/lambdaSR/t0eoswrb1corrAGR1.f
c------------------------------------------------------------
C	PROGRAM PARA DETERMINAR A EQUACAO DE ESTADO NO MODELO
C	DO WALECKA NAO LINEAR COM MESONS RHO 
C       PARA T=0 e beta-eq
c-------------------------------------------------------
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	DIMENSION X(5),rml(2),RNL(2),fvec1(5)
        real*8 X0a(1),X1a(1),XMINa(1),XMAXa(1),Za(1)
	EXTERNAL F3P,F3N,F4P,F4N,f4,f3
        logical check
        CHARACTER*3 ichoice,ipne,ipnc,ipn
        COMMON/MA/RMA
        COMMON/MUN/Rhob,RMUP,rkp,rkn
	COMMON/ROS/RHOP,RHON
        common/field/RNL,RV0,B0,rho0,phi0
	COMMON/CT2/GAMMA,RMS,RMV,rmrho,RML,rm,hc
	COMMON/CT/PI2,GS,GV,GRHO,RKA,RLAMBDA,XSI,gwr
	COMMON/MAIN/RME
        COMMON/PDN/stepmax,acc,X0a,X1a,XMINa,XMAXa
        DATA RML/0.511d0,105.66d0/
C-----------------------------------------------------------
C	INPUT PARAMETERS AND CONSTANTS
C-----------------------------------------------------------
	PI=DACOS(-1.D0)
        PI2=PI*PI
        hc=197.326d0
        write(6,*)'ip: 1=NL3,2-NL3wr,3-FSU2R'
        read(5,*)ip
c        ip=1
c---------------------------------------------------------
c$$$        open(unit=73,file='AGR7f2rT0.73',status='unknown')
c$$$        write(73,*) '#1-rho(fm**-3), 2-rhop, 3-rhon, 4-Pb(MeV.fm**-3), 
c$$$     &5-Pe(MeV.fm**-3), 6-Eb/A(MeV), 7-Ee(MeV), 8-mue, 9-mun, 10-mup, 11
c$$$     &-gwr'
C-------------------------------------------------------
c       NL3 SET OF PARAMETERS (PRC55,1997,540), L=118 MeV
C--------------------------------------------------------
        if(ip.eq.1)then
           RM=939.D0
           RMS=508.194D0/RM
           RMV=782.501D0/RM
           RMRHO=763.D0/RM
           GS=10.217D0
           GV=12.868D0
           GRHO=8.948D0
           RKA=2.D0*10.431D0*hc/RM
           RLAMBDA=-6.D0*28.885D0
           XsI=0.0
           gwr=0.d0
        endif
C---------------------------------------------------------
c     ctes Horowitz+ Piekarewicz --NL3wr
c------------------------------------------------------
        if(ip.eq.2)then
         RM=939.D0
         Rms=508.194D0/RM
         RMV=782.501D0/RM
         RMRHO=763.D0/RM
         GS=10.217D0
         GV=12.868D0
         RKA=2.D0*10.431D0*hc/RM
         RLAMBDA=-6.D0*28.885D0 
         XSI=0.0d0
         write(6,*)'ipar=?'
         read(5,*)ipar
        if (ipar.eq.1)then !L=88 MeV
            GRHO=dsqrt(90.941066817193501d0)
            gwr=0.01d0
         endif
         if (ipar.eq.2)then !L=68 MeV
            GRHO=dsqrt(106.03566641594703d0)
            gwr=0.02d0
         endif
         if (ipar.eq.3)then !L=61 MeV
            GRHO=dsqrt(115.63210260639855d0)
            gwr=0.025d0
         endif
         if (ipar.eq.4)then !L=55 MeV
            GRHO=dsqrt(127.13837543370576d0)
            gwr=0.03d0
         endif
       endif
c------------------------------------------
c     FSU2R (Tolos et al Astrophys J 834 3 2017)
c------------------------------------------
      if(ip.eq.3)then
	 RM=939.d0
         Rms=497.479D0/rm
         RMV=782.50D0/rm
         RMRHO=763.D0/rm
         GS=dsqrt(107.5751D0)
         GV=dsqrt(182.3949d0)
	 grho=dsqrt(247.3409d0)
         rka=3.0911*gs**3/rm
         rLAMBDA=-0.00168*gs**4
         XSI=0.024d0         
         gwr=0.05d0
         rho0=0.1505d0/((rM/HC)**3)
      endif
c---------------------------------------------------------
c     ENTRADA DE DADOS
c---------------------------------------------------------
c        write(6,*)'rho=?'
c        read(5,*)rho1
c---------------------------------------------------------
C-----------------------------------------------------------------
        do i=1,2
          rml(i)=rml(i)/rm
        enddo
	IGAMMA=2
	GAMMA=DFLOAT(IGAMMA)
        a1=0.1D0
        a2=0.6d0
        a3=0.98d0
        a4=-0.1d0*(hc/rm)**3
c-----------------------------------------------
         X(1)=dASIN(dSQRT(a1))
        x(2)=dASIN(dSQRT(a2))
        X(4)=DSQRT(a3)
        X(5)=dsqrt(-a4)
         RMUE=dsin(X(1))**2
        rma=dsin(X(2))**2       
        rmun=X(4)**2             
        rho3=-X(5)**2          
        r0=0.08122826*(hc/rm)**3
        DO J=1,300!655 !1
c$$$           xlog=-3.d0+0.005d0*(j-1)
c$$$           rhob=10**xlog
c$$$           rhob=rhob*(hc/rm)**3
c--------------------------------------------------------------------
          rhob=0.001d0*(hc/rm)**3+DFLOAT(J-1)*0.005*(hc/rm)**3
c-----------------------------------------------------------------
c$$$           rhob=r0+DFLOAT(J-1)*0.001*(hc/rm)**3 
c----------------------------------------------------------
c$$$           rhob=rho1*(hc/rm)**3
c---------------------------------------------------------------
           write(6,*)'rho=',rhob/(hc/rm)**3
           x(3)=rhob*gv**2/rmv**2
C--------------------------------------------------------
        CALL newt(X,5,CHECK)
        IF(CHECK)PAUSE 'There is no root in broydn...'
         RMUE=dsin(X(1))**2
         rma=dsin(X(2))**2      
         rmun=X(4)**2           
         rho3=-X(5)**2  
         rv0=x(3)/gv
c-----------------------------------------------------------
C     DEFININDO AS OUTRAS QUANTIDADES (AINDA ADIMENSIONAIS)
C--------------------------------------------------
         RME=RMA
         PHI0=(-RMA+1.D0)/GS
         RHO=rhob
         rmrho2=rmrho**2+2*gwr*grho**2*gv**2*rv0**2 
         B0=GRHO/(2.D0*RMRHO2)*(rhop-rhon)
         rkp=(3.d0*pi2*rhop)**(1.d0/3.d0)
         rkn=(3.d0*pi2*rhon)**(1.d0/3.d0)
c-------------------------------------------------------------------
C     DEFININDO ENERGIA, PRESSAO, ENTROPIA
C------------------------------------------------------------------
           CALL GAUSS(F3P,0.D0,RKP,10,RES3P,II)
           CALL GAUSS(F3N,0.D0,RKN,10,RES3N,II)
           ENER=1.d0/PI2*(RES3P+RES3N) 
     &          +RMS**2*phi0**2/2.D0+RKA*(phi0)**3/6.D0 
     &          + RLAMBDA*(phi0)**4/24.D0
     &          -RMV**2*rV0**2/2.D0-XSI*GV**4*rV0**4/24.D0  
     &          -RMRHO**2*B0**2/2.D0
     &          +GV*rV0*rhob+GRHO/2.D0*B0*RHO3
     &          -gwr*grho**2*gv**2*b0**2*rv0**2
           DENER=ENER/RHOB-1.D0
           DENER2=ENER/RHOB
           CALL GAUSS(F4P,0.D0,RKP,10,RES4P,II)
           CALL GAUSS(F4N,0.D0,RKN,10,RES4N,II)
           PRESS=1.d0/(3.D0*PI2)*(RES4P+RES4N)
     &          -RMS**2*phi0**2/2.D0
     &          - RKA*phi0**3/6.D0-RLAMBDA*phi0**4/24.D0
     &          +RMV**2*rV0**2/2.D0+XSI*GV**4*rV0**4/24.D0
     &          +RMRHO**2*B0**2/2.D0
     &          +gwr*grho**2*gv**2*b0**2*rv0**2
           pressl=0.d0
           enerlep=0.d0
           DO I=1,2 
              RMe=RML(i)
              RNU=RMUE
              rkfl2=rmue**2-rme**2              
              IF(rkfl2.GT.0.d0)THEN                                  
                 rkfl=dSQRT(rkfl2)     
              else
                 rkfl=0.d0
              endif     
              CALL GAUSS(F4,0.D0,rkfl,10,RE4l,II)
              PRESSl=PRESSl+RE4l/(3.D0*PI2)             
              CALL GAUSS(F3,0.D0,Rkfl,10,RE3L,II)
              ENERLEP=ENERLEP+RE3L/(PI2)
c             write(6,*)'rme',i,rkfl,rmue,rme,rkfl,pressl
           ENDDO
C-----------------------------------------------------------------
C     ACERTANDO AS DIMENSOES
C---------------------------------------------------------------
           bind=dener+enerlep/rhob
           ENER1=ENER*RM*(RM/197)**3
           ENERe1=ENERlep*RM*(RM/197)**3
           RHOP=RHOP*(RM/hc)**3
           RHON=RHON*(RM/hc)**3
           DENS=RHOP+RHON
           RMUP=RMUP*RM
           RMUN=RMUN*RM	
           rmue=rmue*rm
           PRESS1=PRESS*RM*(RM/hc)**3
           PRESSl1=PRESSl*RM*(RM/hc)**3
           DENER=DENER*RM
           denere=enerlep/rhob*rm
           DENER2=DENER2*RM
           bind=bind*RM
            ENERtov=(ENER+enerlep)*(RM/hc)**4
           PRESStov=(PRESS+pressl)*(RM/hc)**4
           pressureh=presstov*hc
           write(70,*)dens,bind,ener1,enere1,press1,pressl1         
c---------------------------------------------------------------------------
        write(73,321)dens,rhop,rhon,press1,pressl1,dener,denere,rmue,
     1       rmun,rmup,gwr,bind,press1+pressl1
c---------------------------------------------------------------------
c----------------------------------------------------------------------        
c----This is the file that will be the input of your TOV code to get M(R)---------------------------------------------------------------------
        write(6,*)dens,enertov,presstov
        write(71,*)dens,enertov,presstov
c-------------------------------------------------------------------------
      write(75,*)dens,rhon,rhop,rhop/dens,pressureh
c----------------------------------------------------------------------
 32   format(2x,10(d11.4,2x))
 321  format(2x,3(f11.8,2x),4(f11.7,2x),10(e11.4E2,2x))
	ENDDO
c$$$        close(73)
	STOP
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
C----------------------------------------------------------------
        FUNCTION F1P(X)
        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON/MA/RMA
        E=DSQRT(X*X+RMA*RMA)
        F1P=X*X/E
        RETURN
        END
c---------------------------------------------------------------------------
        FUNCTION F1N(X)
        IMPLICIT REAL*8 (A-H,O-Z)
        COMMON/MA/RMA
        E=DSQRT(X*X+RMA*RMA)
        F1N=X*X/E
        RETURN
        END
C------------------------------------------------------------------
	FUNCTION F3(X)
	IMPLICIT REAL*8 (A-H,O-Z)
	COMMON/MAIN/RME
        E=DSQRT(X*X+RME*RME)
        F3=X*X*E
        RETURN
	END
C------------------------------------------------------------------
	FUNCTION F3P(X)
	IMPLICIT REAL*8 (A-H,O-Z)
	COMMON/MAIN/RME
	E=DSQRT(X*X+RME*RME)
        F3P=X*X*E
        RETURN
	END
C-----------------------------------------------------------------
	FUNCTION F3N(X)
	IMPLICIT REAL*8 (A-H,O-Z)
	COMMON/MAIN/RME
	E=DSQRT(X*X+RME*RME)
        F3N=X*X*E
        RETURN
	END
C------------------------------------------------------------------
	FUNCTION F4P(X)
	IMPLICIT REAL*8 (A-H,O-Z)
	COMMON/MAIN/RME
	E=DSQRT(X*X+RME*RME)
        F4P=X**4/E
        RETURN
        END
C------------------------------------------------------------------
	FUNCTION F4(X)
	IMPLICIT REAL*8 (A-H,O-Z)
	COMMON/MAIN/RME
        E=DSQRT(X*X+RME*RME)
        F4=X**4/E
        RETURN
        END
c------------------------------------------------------------------
 	FUNCTION F4N(X)
	IMPLICIT REAL*8 (A-H,O-Z)
	COMMON/MAIN/RME
	E=DSQRT(X*X+RME*RME)
        F4N=X**4/E
        RETURN
        END
C-------------------------------------------------------------
        SUBROUTINE funcv(N,X,FVEC)
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION x(5),RML(2),RNL(2),fvec(5),fvec1(5)
         real*8 X0a(1),X1a(1),XMINa(1),XMAXa(1),Za(2)
        external F1P,F1N,ya,yda
        COMMON/MA/RMA
        COMMON/MUN/Rhob,RMUP,rkp,rkn
	COMMON/ROS/RHOP,RHON
        common/field/RNL,RV0,B0,rho0,phi0
	COMMON/CT2/GAMMA,RMS,RMV,rmrho,RML,rm,hc
        COMMON/CT/PI2,GS,GV,GRHO,RKA,RLAMBDA,XSI,gwr
        COMMON/PDN/stepmax,acc,X0a,X1a,XMINa,XMAXa
        rmue=dsin(X(1))**2      
        rma=dsin(X(2))**2       
        rmun=X(4)**2            
        rho3=-X(5)**2           
        rv0=x(3)/gv
        rho=rhob
        rhop=0.5d0*(rhob+rho3)
        rhon=0.5d0*(rhob-rho3)
        rmup=rmun-rmue
c----------------------------------------------------------
c     R0 and B0 fields
C-----------------------------------------------------
        PHI0=(-RMA+1.D0)/GS
        rmrho2=rmrho**2+2*gwr*grho**2*gv**2*rv0**2
	B0=GRHO/(2.D0*RMRHO2)*rho3
C----------------------------------------------------------------        
C     DEFININDO POTENCIAIS QUIMICOS efectivos
C------------------------------------------------------------------
        RnUP=RmUP-GV*RV0-GRHO*B0/2.D0
        RnUN=RmUN-GV*RV0+GRHO*B0/2.D0
        RKp2=RnUp**2-RMA**2
        IF(rkp2.GT.0.d0)THEN                                  
           rkp=dSQRT(rkp2)     
        else
           rkp=0.d0
        endif     
        rhop=rkp**3/3.d0/pi2
        RKn2=RnUn**2-RMA**2
        IF(rkn2.GT.0.d0)THEN                                  
           rkn=dSQRT(rkn2)     
        else
           rkn=0.d0
        endif     
        rhon=rkn**3/3.d0/pi2
c---------------------------------------
c     leptons
c     
        DO ii=1,2
           rkfl2=rmue**2-rml(ii)**2
           IF(rkfl2.GT.0.d0)THEN                                  
              rkfl=dSQRT(rkfl2)     
           else
              rkfl=0.d0
           endif     
           RNL(ii)=rkfl**3/3.d0/pi2
        ENDDO
        chargel=rnl(1)+rnl(2)
C-----------------------------------------------------------
        CALL GAUSS(F1P,0.D0,RKP,10,RES1P,II)
        CALL GAUSS(F1N,0.D0,RKN,10,RES1N,II)
        RHOS=GAMMA/(2.D0*PI2)*RMA*(RES1P+RES1N) 
        Fvec(1)=GS/(RMS**2)*RHOS-RKA/(2.D0*RMS**2)*
     &       PHI0**2-RLAMBDA/(6.D0*RMS**2)*PHI0**3-phi0
        fvec(2)=gv**2/rmv**2
     &       *(rho-xsi/6.d0*(x(3))**3
     &       -2*gwr*x(3)*grho**2*b0**2)-x(3)
        fvec(3)=rhop-(rho+rho3)/2.d0 
        fvec(4)=rhon-(rho-rho3)/2.d0
        fvec(5)=chargel-rhop
         return
        end
c-----------------------------------------------------------------------
c--------------------------------------------------------------------------
c        include 'newt.f'
c--------------------------------------------------------------------------
c--------------------------------------------------------------------------
C=======================================================================
C     Other subroutines from Numerical Recipes
C=======================================================================

C ... Given an initial guess x(1:n) for a root in n dimensions, find the root by a
C ... globally convergent Newton's method. The vector of functions to be zeroed,
C ... called fvec(1:n) in the routine below, is returned by a user-supplied
C ... subroutine that must be called funcv and have the declaration subroutine
C ... funcv(n,x,fvec). The output quantity check is false on a normal return and
C ... true if the routine has converged to a local minimum of the function fmin
C ... defined below. In this case try restarting from a different initial guess.
C ... Parameters: NP is the maximum expected value of n; MAXITS is the maximum
C ... number of iterations; TOLF sets the convergence criterion on function
C ... values; TOLMIN sets the criterion for deciding whether spurious convergence
C ... to a minimum of fmin has occurred; TOLX is the convergence criterion on ffix;
C ... STPMX is the scaled maximum step length allowed in line searches. 
      SUBROUTINE newt(x,n,check) 
      INTEGER n,nn,NP,MAXITS 
      LOGICAL check 
      REAL*8 x(n),fvec,TOLF,TOLMIN,TOLX,STPMX 

c      PARAMETER(NP=40,MAXITS=2000,TOLF=1.e-10,TOLMIN=1.e-12,TOLX=1.e-12, 
c     &          STPMX=100.)
c-------------------------------------------------------------------------------
c      PARAMETER(NP=40,MAXITS=200000,TOLF=1.e-8,TOLMIN=1.e-10,TOLX=1.e-10 
c     &          ,STPMX=100.)
c-----------------------------------------------------------------------------------
      PARAMETER(NP=40,MAXITS=200000,TOLF=1.e-8,TOLMIN=1.e-10,TOLX=1.e-10 
     &          ,STPMX=100.)
c      PARAMETER(NP=40,MAXITS=2000,TOLF=1.e-10,TOLMIN=1.e-12,TOLX=1.e-12, 
c     &          STPMX=100.)
c      PARAMETER (NP=40,MAXITS=200,TOLF=1.e-4,TOLMIN=1.e-6,TOLX=1.e-7, 
c     &          STPMX=100.)

      COMMON /newtv/ fvec(NP),nn 
      SAVE /newtv/ 

c ... test

      common/iter/its 

c ... 



C USES fdjac,fmin,lnsrch,lubksb,ludcmp


      INTEGER i,its,j,indx(NP) 
      REAL*8 d,den,f,fold,stpmax,sum,temp,test,fjac(NP,NP), 
     & g(NP),p(NP),xold(NP),fmin

      EXTERNAL fmin 

      nn=n 

      f=fmin(x) 
      test=0. 

      do 11 i=1,n 
      if(abs(fvec(i)).gt.test) test=abs(fvec(i)) 
 11   continue

      if(test.lt..01*TOLF) then
      check=.false. 
      return 
      endif 
      sum=0. 

      do 12 i=1,n
      sum=sum+x(i)**2 
 12   continue
 
      stpmax=STPMX*max(sqrt(sum),float(n)) 
      do 21 its=1,MAXITS 

      call fdjac(n,x,fvec,NP,fjac) ! Numerical
c      call fdjac(n,x,fjac) ! Analytical

      do 14 i=1,n 
      sum=0. 
      do 13 j=1,n
      sum=sum+fjac(j,i)*fvec(j) 
 13   continue
      g(i)=sum 
 14   continue
      do 15 i=1,n 
      xold(i)=x(i) 
 15   continue
      fold=f 
      do 16 i=1,n 
      p(i)=-fvec(i) 
 16   continue

      call dludcmp(fjac,n,NP,indx,d)
      call dlubksb(fjac,n,NP,indx,p) 
      call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)

      test=0. 
      do 17 i=1,n

      if(abs(fvec(i)).gt.test) test=abs(fvec(i)) 

 17   continue 

      if(test.lt.TOLF) then
       check=.false. 
      return 
      endif 
      if(check)then 
      test=0. 
      den=max(f,.5*n) 
      do 18 i=1,n
      temp=abs(g(i))*max(abs(x(i)),1.)/den 
      if(temp.gt.test)test=temp 
18    continue
      if(test.lt.TOLMIN)then
      check=.true. 
      else
      check=.false. 
      endif 
      return

      endif 
      test=0. 

      do 19 i=1,n
      temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.) 
      if(temp.gt.test)test=temp 
 19   continue
      if(test.lt.TOLX) return 
 21   continue 
      pause 'MAXITS exceeded in newt' 
      return
      END
c----------------------------------------------------------------------
	SUBROUTINE DLUBKSB(A,N,NP,INDX,B)
C..............................................................
C	Routine to solve the linear system
C	A. X = B
C	where A is the matrix already in LU decomposition form
C	(see DLUDCMP routine)
C	and B is the rhs vector.
C
C	On output, B contains the solution
C
C	The way of solving a system of linear equations is
C	Call DLUDCMP to change A into its LU decomposition
C	call afterwards DLUBCKS (backsubstitution)
C
C	A	dimensioned at (NP,NP), logical NxN
C	INDX	order of rows (see DLUDCMP), dimensioned to N
C	B	Column vector, dimensioned (at least) to N
C
C	A is not destroyed, so this routine may be called
C	many times, per one call of DLUDCMP
C
C..............................................................
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION A(NP,NP),INDX(N),B(N)
	II=0
	DO 12 I=1,N
	LL=INDX(I)
	SUM=B(LL)
	B(LL)=B(I)
	IF (II.NE.0)THEN
	DO 11 J=II,I-1
	SUM=SUM-A(I,J)*B(J)
11	CONTINUE
	ELSE IF (SUM.NE.0.D0) THEN
	II=I
	ENDIF
	B(I)=SUM
12	CONTINUE
	DO 14 I=N,1,-1
	SUM=B(I)
	IF(I.LT.N)THEN
	DO 13 J=I+1,N
	SUM=SUM-A(I,J)*B(J)
13	CONTINUE
	ENDIF
	B(I)=SUM/A(I,I)

14	CONTINUE
	RETURN
	END
!-----------------------------------------------------------------------               
	SUBROUTINE DLUDCMP(A,N,NP,INDX,D)
C................................................................
C	Given a NxN matrix A, with PHYSICAL dimension NP
C	the routine replaces it by its Lower-Upper LU
C	decomposition, with a row-wise ordering given by
C	INDX (dimensioned at N) and sign exchange D
C	Used with DLUBKSB to solve linear systems or
C	inverting a matrix
C
C	Parameters
C	A	matrix dimensioned at (NP,NP), logical NxN
c	INDX	integer array dimensioned to N at least
C	D	Real*8 number indicating the sign of permutations
C
C	Limitations: N max is 250
C
C	A is destroyed
C................................................................
C
      IMPLICIT REAL*8 (A-H,O-Z)
	     PARAMETER (NMAX=250,TINY=1.0D-20)
	     DIMENSION A(NP,NP),INDX(N),VV(NMAX)
	     D=1.D0
	     DO 12 I=1,N
	     AAMAX=0.D0
	     DO 11 J=1,N
	     IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
   11	CONTINUE
	     IF (AAMAX.EQ.0.D0) PAUSE 'Singular matrix.'
	     VV(I)=1.D0/AAMAX
   12	CONTINUE
	     DO 19 J=1,N
	     IF (J.GT.1) THEN
	     DO 14 I=1,J-1
	     SUM=A(I,J)
	     IF (I.GT.1)THEN
	     DO 13 K=1,I-1
	     SUM=SUM-A(I,K)*A(K,J)
   13	CONTINUE
	     A(I,J)=SUM
	     ENDIF
   14	CONTINUE
	     ENDIF
	     AAMAX=0.D0
	     DO 16 I=J,N
	     SUM=A(I,J)
	     IF (J.GT.1)THEN
	     DO 15 K=1,J-1
	     SUM=SUM-A(I,K)*A(K,J)
   15	CONTINUE
	     A(I,J)=SUM
	     ENDIF
	     DUM=VV(I)*ABS(SUM)
	     IF (DUM.GE.AAMAX) THEN
	     IMAX=I
	     AAMAX=DUM
	     ENDIF
   16	CONTINUE
	     IF (J.NE.IMAX)THEN
	     DO 17 K=1,N
	     DUM=A(IMAX,K)
	     A(IMAX,K)=A(J,K)
	     A(J,K)=DUM
   17	CONTINUE
	     D=-D
	     VV(IMAX)=VV(J)
	     ENDIF
	     INDX(J)=IMAX
	     IF(J.NE.N)THEN
	     IF(A(J,J).EQ.0.D0) A(J,J)=TINY
	     DUM=1.D0/A(J,J)
	     DO 18 I=J+1,N
	     A(I,J)=A(I,J)*DUM
   18	CONTINUE
	     ENDIF
   19	CONTINUE
	     IF(A(N,N).EQ.0.D0)A(N,N)=TINY
	     RETURN
	     END
!-----------------------------------------------------------------------               
!-----------------------------------------------------------------------               
       SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func) 
!-----------------------------------------------------------------------               
C ... Given an n-dimensional point xold(1:n), the value of the function and
C ... gradient there, fold and g(1:n), and a direction p(1:n), finds a new point
C ... x(1:n) along the direction p from xold where the function func has decreased
C ... "sufficiently." The new function value is returned in f. stpmax is an input
C ... quantity that limits the length of the steps so that you do not try to
C ... evaluate the function in regions where it is undefined or subject to
C ... overflow. p is usually the Newton direction. The output quantity check is
C ... false on a normal exit. It is true when x is too close to xold. In a
C ... minimization algorithm, this usually signals convergence and can be ignored.
C ... However, in a zero-finding algorithm the calling program should check
C ... whether the convergence is spurious. Parameters: ALF ensures sufficient
C ... decrease in function value; TOLX is the convergence criterion on \Delta x.

C USES func
!-----------------------------------------------------------------------               
       INTEGER n 
       LOGICAL check 
       REAL*8 f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF,TOLX 
       PARAMETER (ALF=1.e-4,TOLX=1.e-7) 
       EXTERNAL func 

       INTEGER i 
       REAL*8 a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope, 
     & sum,temp,test,tmplam

       check=.false. 
       sum=0. 

       do 11 i=1,n
        sum=sum+p(i)*p(i) 
 11    continue 
       sum=sum**(1.d0/2.d0) 
       if(sum.gt.stpmax) then        
         do 12 i=1,n
           p(i)=p(i)*stpmax/sum 
 12    continue 
       endif 
       slope=0. 
       do 13 i=1,n
          slope=slope+g(i)*p(i) 
 13    continue
      if(slope.ge.0.) pause 'roundoff problem in lnsrch' 
      test=0.                       

      do 14 i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.)
        if(temp.gt.test)test=temp 
 14   continue
      alamin=TOLX/test
      alam=1.               
  1   continue             

      do 15 i=1,n
        x(i)=xold(i)+alam*p(i) 
 15   continue

      f=func(x) 

      if(alam.lt.alamin) then 
                              
      do 16 i=1,n
      x(i)=xold(i) 
 16   continue 
      check=.true. 
      return 
      elseif(f.le.fold+ALF*alam*slope) then  
      return 
      else                                    
      if(alam.eq.1.) then                     
        tmplam=-slope/(2.*(f-fold-slope)) 
      else                                    

      rhs1=f-fold-alam*slope 
      rhs2=f2-fold-alam2*slope
      a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
      b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/ 
     &   (alam-alam2)


      if(a.eq.0.) then
       tmplam=-slope/(2.*b) 
      else
       disc=b*b-3.*a*slope 
      if(disc.lt.0.)then
      tmplam=.5*alam 
      elseif(b.le.0.)then
      tmplam=(-b+disc**(1.d0/2.d0))/(3.*a) 
      else
      tmplam=-slope/(b+disc**(1.d0/2.d0))  
      endif 
      endif 
      if(tmplam.gt..5*alam)tmplam=.5*alam

      endif 
      endif 

      alam2=alam 
      f2=f 
      alam=max(tmplam,.1*alam)

      goto 1 
      return
      END
!-----------------------------------------------------------------------               
C-----------------------------------------------------------------------
      SUBROUTINE fdjac(n,x,fvec,np,df)
      INTEGER n,np,NMAX
      DOUBLE PRECISION df(np,np),fvec(n),x(n),EPS
      PARAMETER (NMAX=40,EPS=1.d-8)
CU    USES funcv
      INTEGER i,j
      DOUBLE PRECISION h,temp,f(NMAX)
      do 12 j=1,n
        temp=x(j)
        h=EPS*abs(temp)
        if(h.eq.0.d0)h=EPS
        x(j)=temp+h
        h=x(j)-temp
        call funcv(n,x,f)
        x(j)=temp
        do 11 i=1,n
          df(i,j)=(f(i)-fvec(i))/h
11      continue
12    continue
      return
      END
C-----------------------------------------------------------------------
      FUNCTION fmin(x)
      INTEGER n,NP
      DOUBLE PRECISION fmin,x(*),fvec
      PARAMETER (NP=40)
      COMMON /newtv/ fvec(NP),n
      SAVE /newtv/
CU    USES funcv
      INTEGER i
      DOUBLE PRECISION sum
      call funcv(n,x,fvec)
      sum=0.d0
      do 11 i=1,n
        sum=sum+fvec(i)**2
11    continue
      fmin=0.5d0*sum
      return
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------




