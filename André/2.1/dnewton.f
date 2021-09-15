c---------------------------------------------------------------------------
c-----dnewto.for---911111 Coimbra 860117 Regensburg---Alex H.Blin-------------
      subroutine dnewto(x00,x11,xmin,xmax,idim,y,yp,acc,stemax,ilim,
     1                  x1,iact,ierr)
c-----Newton iteration procedure to find array x as zero of y(x)
c-----to an accuracy acc; real*8 version.
c-----if y' is a known function, set x0(*)=x1(*) and supply fct. yp;
c-----if y' is not known, supply some fct. yp (it will not be called).
c-----input:  x00(*),x11(*)= initial guesses for solution array, if
cx00=x11
c-----        then y' is calculated using function yp, otherwise
cnumerically;
c-----        xmin(*),xmax(*)= range in which search is to be confined,
cwith
c-----        restriction xmin<xmax;
c-----        idim = dimension of arrays x00,x11,xmin,xmax,x (max.
cidimx, see
c-----               data statement);
c-----        y,yp = external function names of y(x) and y'(x), to be
c-----               called as y(i,x,idim,ierr) (i= ith fct.y_i) and 
c-----               yp(i,j,x,idim,ierr) (j= partial derivative with
crespect 
c-----               to x(j)), sample fcts. at the end of this file;
c-----        acc = accuracy to which y(x)=0 to be found;
c-----        stemax = max. stepsize in searching x to avoid insta-
c-----                 bilities;
c-----        ilim = max. # of iterations.
c-----in/out: ierr= debug/error variable: in: < 0 displays debug data,
c-----        out: > 0 means error occurred.
c-----output: x1 = solution;
c-----        iact = actual number of iterations.
c-----note: this routine written for dimension of up to idimx (see
c-----data statement). To increase it, change idimx and arrays in
c-----dimension statement.
c-----note: This routine must be linked with function ddet for the 
c-----      calculation of the determinant of an idim*idim matrix.
c
      implicit real*8 (a-h)
      implicit real*8 (o-z)
c     external y,yp,ddet
      logical anal(6),ldeb
      dimension x00(idim),x11(idim),x1(idim),xmin(idim),xmax(idim),
     1     x0(6),xd(6),y1(6),d(6),xste(6),dydx(6,6),aa(6,6)
      data dz/1.d-100/,dinf/1.d30/,ideb/-0/,idimx/6/
c
      ldeb=ierr.lt.ideb
      if(idim.gt.idimx) goto 206
c     -300-
      do 300 j=1,idim
         if(xmin(j).ge.xmax(j)) goto 207
         anal(j)=x00(j).eq.x11(j)
         x0(j)=x00(j)
         x1(j)=x11(j)
         xste(j)=x1(j)-x0(j)
  300 continue
c     
c     -100-iteration loop
      do 100 i=1,ilim
         iact=i
         s=0.
c     -101-
         do 101 j=1,idim
            y1(j)=y(j,x1,idim,ierr)
            if(ierr.gt.0) goto 201
            s=s+y1(j)**2
 101     continue
         if(ldeb) write(*,*) ' dnewto #1: x1(*),y1(*),s=',x1,y1,s
         if(dsqrt(s).lt.acc) return
c     
c-----derivatives
c     -104-
         do 104 j=1,idim
c     -114-
            do 114 k=1,idim
c     
c-----analytic y'
               if(anal(j)) then
                  dydx(k,j)=yp(k,j,x1,idim,ierr)
                  if(ierr.gt.0) goto 202
                  if(ldeb) write(*,*) 
     1                 ' dnewto #2: k,j,dydx(k,j)=',k,j,dydx(k,j)
c     
c-----numerical y'
               else
c     -124-first store x1 in xd
                  do 124 l=1,idim
 124                 xd(l)=x1(l)
c-----now replace the j'th component by x0(j)
                     xd(j)=x0(j)
                     y0=y(k,xd,idim,ierr)
                     if(ierr.gt.0) goto 201
c-----new dydx only if xste>0, otherwise keep old value
              if(dabs(xste(j)).gt.dz) dydx(k,j)=(y1(k)-y0)/xste(j)
              if(ldeb) write(*,*) ' dnewto #3: k,j,x0(j),y0,dydx(k,j)=',
     1                    k,j,x0(j),y0,dydx(k,j)
                  endif
c     
                  a=dabs(dydx(k,j))
                  if(a.gt.dinf) goto 204
                  if(a.lt.dz) dydx(k,j)=dsign(dz,dydx(k,j))
                  if(a.eq.0.) dydx(k,j)=dz
 114           continue
 104        continue
c     
c-----idim*idim_determinants
            dn=ddet(dydx,idimx,idim,ierr)
            if(ierr.gt.0) goto 208
            if(dabs(dn).lt.dz) goto 205
c     -400-determinant in the numerator
            do 400 ii=1,idim
c     -410-first prepare matrix of dydx
               do 410 jj=1,idim
                  do 410 kk=1,idim
 410                 aa(jj,kk)=dydx(jj,kk)
c     -420-now overwrite the appropriate column with vector y1
                     do 420 jj=1,idim
 420                    aa(jj,ii)=-y1(jj)
                        d(ii)=ddet(aa,idimx,idim,ierr)
                        if(ierr.gt.0) goto 208
 400                 continue
c-----end of idim*idim_determinants
c     
                     s=0.
c     -103-
                     do 103 j=1,idim
                        x0(j)=x1(j)
                        xste(j)=d(j)/dn
              if(dabs(xste(j)).gt.stemax) xste(j)=dsign(stemax,xste(j))
                        x1(j)=x1(j)+xste(j)
                        if(x1(j).gt.xmax(j)) x1(j)=xmax(j)
                        if(x1(j).lt.xmin(j)) x1(j)=xmin(j)
                        s=s+xste(j)**2
 103                 continue
                     if(dsqrt(s).lt.dz) goto 203
 100              continue
c     
c-----errors
c     
 209  write(*,*) ' error 9 in dnewto: #_of_iterations too large, iact='
     1                 ,iact
         ierr=9
         return
c     
 208    write(*,*) ' error 8 in dnewto: error in ddet, ierr=',ierr
                  ierr=8
                  return
c     
 207    write(*,*) ' error 7 in dnewto: xmin not < xmax, xmin=',
     1                 xmin,' xmax=',xmax
                  ierr=7
                  return
c     
 206    write(*,*) ' error 6 in dnewto: dimension idim=',idim,
     1                 ' too large'
                  ierr=6
                  return
c     
 205    write(*,*) ' error 5 in dnewto: system of eqs. has no solution',
     1                 ' dn=',dn
                  ierr=5
      return
c     
 204  write(*,*) ' error 4 in dnewto: y'' too large, dydx(k,j)=',
     1     dydx(k,j)
      ierr=4
      return
c     
 203  write(*,*) ' error 3 in dnewto: accuracy cannot be reached, s=',s
      ierr=3
      return
c     
 202  write(*,*) ' error 2 in dnewto: error in yp, ierr=',ierr
      return
c     
 201  write(*,*) ' error 1 in dnewto: error in y, ierr=',ierr
      return
      end
c------------------------------------------------------------------------------
c-----program example----------------------------------------------------------
c     implicit real*8 (a-h)
c     implicit real*8 (o-z)
c     external y,yp
c     dimension x0(2),x1(2),x(2),xmin(2),xmax(2)
c     write(*,*) ' initial values x0(1),x0(2),x1(1),x1(2)?'
c     read (*,*) x0(1),x0(2),x1(1),x1(2)
c     write(*,*) ' accuracy?'
c     read(*,*) acc
c     ierr=0
c     xmin(1)=-100.
c     xmin(2)=-100.
c     xmax(1)= 100.
c     xmax(2)= 100.
c     call dnewto(x0,x1,xmin,xmax,2,y,yp,acc,10.,1000,x,iact,ierr)
c     write(*,*) 'x(1),x(2),iact',x(1),x(2),iact
c     stop
c     end
c------------------------------------------------------------------------------
c-----example for idim=2; solutions: x(1)=2., x(2)=-3.
c     function y(i,x,idim,ierr)
c     implicit real*8 (a-h)
c     implicit real*8 (o-z)
c     dimension x(idim)
c     if(i.eq.1) y=x(1)**4+x(2)**3+11.
c     if(i.eq.2) y=x(1)**2+x(2)-1.
c     return
c     end
c-----------------------------------------------------------------------------
c     function yp(i,j,x,idim,ierr)
c     implicit real*8 (a-h)
c     implicit real*8 (o-z)
c     dimension x(idim)
c     if(i.eq.1.and.j.eq.1) yp=4.*x(1)**3
c     if(i.eq.1.and.j.eq.2) yp=3.*x(2)**2
c     if(i.eq.2.and.j.eq.1) yp=2.*x(1)
c     if(i.eq.2.and.j.eq.2) yp=1.
c     return
c     end
c----------------------------------------------------------------------------
c-----ddet.for---911111 Coimbra 860504 Regensburg---Alex H.Blin---------------
      real*8 function ddet(a,idum,idim,ierr)
c-----calculates determinant of idim*idim_submatrix of idum*idum_array a,
c-----starting always with element a(1,1); real*8 version.
c-----input:  a(*,*)= idum*idum array representing the supermatrix
c-----        idum= dim. of array, max. see idumx in data statement;
c-----        idim= dim. of submatrix, max. idum.
c-----in/out: ierr= debug/error variable, in:<0 displays debug data,
c-----              out:>0 means error occurred.
c-----output: ddet= determinant of a(idim,idim) submatrix.
c-----note:   exchanging rows and columns doesn't matter.
      implicit real*8 (a-h)
      implicit real*8 (o-z)
      external ddet3,ddet4,ddet5,ddet6
      logical ldeb
      dimension a(idum,idum)
      data idumx/6/
c
      if(idum.gt.idumx) goto 101
      if(idim.gt.idum ) goto 102
      ldeb=ierr.lt.0
c     
      if(idim.eq.1) then
         ddet=a(1,1)
c     
      elseif(idim.eq.2) then
         ddet=a(1,1)*a(2,2)-a(1,2)*a(2,1)
c     
      elseif(idim.eq.3) then
         ddet=ddet3(a(1,1),a(1,2),a(1,3),
     1                    a(2,1),a(2,2),a(2,3),
     2                    a(3,1),a(3,2),a(3,3))
c     
      elseif(idim.eq.4) then
         ddet=ddet4(a(1,1),a(1,2),a(1,3),a(1,4),
     1                    a(2,1),a(2,2),a(2,3),a(2,4),
     2                    a(3,1),a(3,2),a(3,3),a(3,4),
     3                    a(4,1),a(4,2),a(4,3),a(4,4))
c     
      elseif(idim.eq.5) then
         ddet=ddet5(a(1,1),a(1,2),a(1,3),a(1,4),a(1,5),
     1                    a(2,1),a(2,2),a(2,3),a(2,4),a(2,5),
     2                    a(3,1),a(3,2),a(3,3),a(3,4),a(3,5),
     3                    a(4,1),a(4,2),a(4,3),a(4,4),a(4,5),
     4                    a(5,1),a(5,2),a(5,3),a(5,4),a(5,5))
      else
         ddet=ddet6(a(1,1),a(1,2),a(1,3),a(1,4),a(1,5),a(1,6),
     1                    a(2,1),a(2,2),a(2,3),a(2,4),a(2,5),a(2,6),
     2                    a(3,1),a(3,2),a(3,3),a(3,4),a(3,5),a(3,6),
     3                    a(4,1),a(4,2),a(4,3),a(4,4),a(4,5),a(4,6),
     4                    a(5,1),a(5,2),a(5,3),a(5,4),a(5,5),a(5,6),
     5                    a(6,1),a(6,2),a(6,3),a(6,4),a(6,5),a(6,6))
      endif
      if(ldeb) write(*,*) ' ddet: a(*,*)=',(i,(j,a(i,j),j=1,idim),
     1                    i=1,idim)
      return
c
c-----errors
  101 write(*,*) ' error in ddet: idum too large, idum=',idum
      ierr=1
      ddet=0.
      return
  102 write(*,*) ' error in ddet: idim too large, idim=',idim
      ierr=2
      ddet=0.
      return
      end
c------------------------------------------------------------------------------
c-----function ddet3 called by ddet--------------------------------------------
      real*8 function ddet3(a11,a12,a13,
     1                                a21,a22,a23,
     2                                a31,a32,a33)
c-----determinant of 3x3 matrix
      implicit real*8 (a-h)
      implicit real*8 (o-z)
      ddet3=a11*a22*a33+a12*a23*a31+a13*a21*a32
     1    -a31*a22*a13-a32*a23*a11-a33*a21*a12
      return
      end
c------------------------------------------------------------------------------
c-----function ddet4 called by ddet--------------------------------------------
      real*8 function ddet4(a11,a12,a13,a14,
     1                                a21,a22,a23,a24,
     2                                a31,a32,a33,a34,
     3                                a41,a42,a43,a44)
c-----determinant of 4x4 matrix
      implicit real*8 (a-h)
      implicit real*8 (o-z)
      external ddet3
      ddet4= a11*ddet3(a22,a23,a24,
     1                           a32,a33,a34,
     2                           a42,a43,a44)
      ddet4=-a12*ddet3(a21,a23,a24,
     1                            a31,a33,a34,
     2                            a41,a43,a44)+ddet4
      ddet4= a13*ddet3(a21,a22,a24,
     1                           a31,a32,a34,
     2                           a41,a42,a44)+ddet4
      ddet4=-a14*ddet3(a21,a22,a23,
     1                            a31,a32,a33,
     2                            a41,a42,a43)+ddet4
      return
      end
c-----function ddet5 called by ddet--------------------------------------------
      real*8 function ddet5(a11,a12,a13,a14,a15,
     1                                a21,a22,a23,a24,a25,
     2                                a31,a32,a33,a34,a35,
     3                                a41,a42,a43,a44,a45,
     4                                a51,a52,a53,a54,a55)
c-----determinant of 5x5 matrix
      implicit real*8 (a-h)
      implicit real*8 (o-z)
      external ddet4
      ddet5= a11*ddet4(a22,a23,a24,a25,
     1                           a32,a33,a34,a35,
     2                           a42,a43,a44,a45,
     3                           a52,a53,a54,a55)
      ddet5=-a12*ddet4(a21,a23,a24,a25,
     1                            a31,a33,a34,a35,
     2                            a41,a43,a44,a45,
     3                            a51,a53,a54,a55)+ddet5
      ddet5= a13*ddet4(a21,a22,a24,a25,
     1                           a31,a32,a34,a35,
     2                           a41,a42,a44,a45,
     3                           a51,a52,a54,a55)+ddet5
      ddet5=-a14*ddet4(a21,a22,a23,a25,
     1                            a31,a32,a33,a35,
     2                            a41,a42,a43,a45,
     3                            a51,a52,a53,a55)+ddet5
      ddet5= a15*ddet4(a21,a22,a23,a24,
     1                           a31,a32,a33,a34,
     2                           a41,a42,a43,a44,
     3                           a51,a52,a53,a54)+ddet5
      return
      end
c-----function ddet6 called by ddet--------------------------------------------
      real*8 function ddet6(a11,a12,a13,a14,a15,a16,
     1                                a21,a22,a23,a24,a25,a26,
     2                                a31,a32,a33,a34,a35,a36,
     3                                a41,a42,a43,a44,a45,a46,
     4                                a51,a52,a53,a54,a55,a56,
     5                                a61,a62,a63,a64,a65,a66)
c-----determinant of 6x6 matrix
      implicit real*8 (a-h)
      implicit real*8 (o-z)
      external ddet5
      ddet6= a11*ddet5(a22,a23,a24,a25,a26,
     1                           a32,a33,a34,a35,a36,
     2                           a42,a43,a44,a45,a46,
     3                           a52,a53,a54,a55,a56,
     4                           a62,a63,a64,a65,a66)
      ddet6=-a12*ddet5(a21,a23,a24,a25,a26,
     1                            a31,a33,a34,a35,a36,
     2                            a41,a43,a44,a45,a46,
     3                            a51,a53,a54,a55,a56,
     4                            a61,a63,a64,a65,a66)+ddet6
      ddet6= a13*ddet5(a21,a22,a24,a25,a26,
     1                           a31,a32,a34,a35,a36,
     2                           a41,a42,a44,a45,a46,
     3                           a51,a52,a54,a55,a56,
     4                           a61,a62,a64,a65,a66)+ddet6
      ddet6=-a14*ddet5(a21,a22,a23,a25,a26,
     1                            a31,a32,a33,a35,a36,
     2                            a41,a42,a43,a45,a46,
     3                            a51,a52,a53,a55,a56,
     4                            a61,a62,a63,a65,a66)+ddet6
      ddet6= a15*ddet5(a21,a22,a23,a24,a26,
     1                           a31,a32,a33,a34,a36,
     2                           a41,a42,a43,a44,a46,
     3                           a51,a52,a53,a54,a56,
     4                           a61,a62,a63,a64,a66)+ddet6
      ddet6=-a16*ddet5(a21,a22,a23,a24,a25,
     1                            a31,a32,a33,a34,a35,
     2                            a41,a42,a43,a44,a45,
     3                            a51,a52,a53,a54,a55,
     4                            a61,a62,a63,a64,a65)+ddet6
      return
      end
c------------------------------------------------------------------------------
c-----program example
c     dimension a(4,4)
c     ierr=0
c     k=0
c     do 100 i=1,4
c     do 100 j=1,4
c     k=k+1
c 100 a(i,j)=k
c     a(1,1)=-1.
c     a(4,4)=-16.
c     do 110 i=1,4
c     d=ddet(a,4,i,ierr)
c 110 write(*,*) d
c     write(*,*) ' -1., -16., 8., -256.'
c     stop
c     end

