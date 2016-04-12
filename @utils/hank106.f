cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        This is the end of the debugging code, and the beginning of
c        the Hankel function code proper
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
ccc        subroutine hank106a(r,h0,h1,w)
        subroutine hank106a(r,h0,h1,w,ninterv)
        implicit real *8 (a-h,o-z)
        integer *4 ninterv
c
c       ninterv is number of intervals used in equispaced
c       subdivision of current interval
c
        complex *16 rk,h0,h1
        dimension w(1)
        save
c
ccc        call hank106b(r,h0,h1,w,w(11) )
        call hank106b(r,h0,h1,w,w(11),ninterv)
c
        return
        end
c
c
c
c
c
ccc        subroutine hank106init(rk7,rmin,rmax,w,keep)
        subroutine hank106init(rk7,rmin,rmax,w,keep,ninterv)
        implicit real *8 (a-h,o-z)
        complex *16 rk,h0,h1,ima,z,rk7,u0,u07,com
        dimension w(1),rea(2),rws(10)
        integer ijw(1),iws(20)
        
c
        equivalence (rea(1),com),(rws(1),iws(1))
        save
c
        data ima/(0.0d0,1.0d0)/
c
c        This subroutine evaluates the Hankel functions 
c        H_0^1, H^1_1 of a complex argument, the argument 
c        living on a ray. The subroutine you are looking at is
c        the initialization subroutine; the evaluation 
c        subroutine is hank106 (see).
c
c
c         8/28/02 - added ninterv to calling sequence  (LG)
c
c
c   PLEASE NOTE THAT THE USE OF THIS SUBROUTIONE IS NOT COMPLETELY
C   STRAIGHTFORWARD: A STRAIGHTFORWARD SUBROUTINE TO USE IS HANK103
C   (SEE). HOWEVER, THIS SUBROUTINE IS ABOUT 4 TIMES FASTER THAN 
C   HANK103.
c
c        Recommended pairs rmin,rmax (assuming that rk \sim 1):
c
c
c
c        rmin       rmax
c
c        0.062      0.125
c        0.031      0.062
c        0.062       0.125
c        0.125       0.25
c        0.25        0.5
c        0.5         1
c        1           2
c        2           5
c        5           10
c        10          100
c        10          200
c
c  GENERALLY, WITH ABS(RMIN*RK) > 10, THERE IS NO SHARP LIMIT ON
C  RMAX. HOWEVER, FOR SUFFICIENTLY LARGE (RMAX-RMIN)*RK, THE 
C  code loses speed due to caching problems. THE DETERIORATION
C  BECOMES NOTICEABLE AT SOME POINT AFTER |(RMAX-RMIN)*RK| > 100
C  (ON THE PENTIUM-IV DESKTOP).
C
c
c                   Input parameters:
c
c  rk7 - the Helmholtz coefficient
c  rmin - the minimum r for which this subroutine will evaluate 
c       the Hankel functions
c  rmax - the minimum r for which this subroutine will evaluate 
c       the Hankel functions
c
c                   Output parameters:
c
c  w - contains various data to be used by the entry hank106 (see)
c  keep - the first keep elements of the array w should not be 
c       changed between the call to this entry, and the subsequent
c       calls to the entry hank106.
c  ninterv - number of intervals used in equispaced
c       subdivision of current interval
c
c
        rk=rk7
c
        dfool=rk/abs(rk)
        if(abs(dfool) .lt. 1.0d-30) rk=rk+2*1.0d-30*abs(rk)
c
        d=rk
c
        x1=d*rmin
        x2=d*rmax
c
c        initialize the evaluation of H_0, H_1
c
        n=11
        ddd=(abs(rmax*rk)-abs(rmin*rk))*2
        i=ddd
c

ccc        call prin2('ddd=*',ddd,1)
        if(i .lt. 10) i=10
        ninterv=i
c
ccc        call prinf('ninterv as calculated*',ninterv,1)
c
        d=rk
        d2=-ima*rk
        coef=d2/d
c
c       allocate memory for the initialization
c
        icenters=21
        lcenters=ninterv*2+4
c
        ih0s=icenters+lcenters
        lh0s=ninterv+4
        lh0s=lh0s*2
c
        ih1s=ih0s+lh0s
        lh1s=ninterv+4
        lh1s=lh1s*2
c
        ih0derss=ih1s+lh1s
        lh0derss=(ninterv*n+4)*2
c
        ih1derss=ih0derss+lh0derss
        lh1derss=(ninterv*n+4)*2
c
        keep=ih1derss+lh1derss
c
        call hank106ini(coef,x1,x2,ninterv,n,
     1      w(icenters),w(ih0derss),w(ih0s),w(ih1s),h,
     2      w(ih1derss),u07)
c
        u0=1/u07
c
c       store in the beginning of the array w various types of data 
c
        ix1=1
        ih=2
        iu0=3
        in=5
c
        irk=6
c
        w(ix1)=x1
        w(ih)=h
        w(ih)=1/h
c
        com=u0
        w(iu0)=rea(1)
        w(iu0+1)=rea(2)
c
        w(in)=n+0.1
c
        com=rk
c
        w(irk)=rea(1)
        w(irk+1)=rea(2)
c
c        store integer data in the array w
c
        iws(1)=ih0derss
        iws(2)=ih0s
        iws(3)=ih1derss
        iws(4)=ih1s
c
        do 3200 j=1,8
c
        w(10+j)=rws(j)
 3200 continue
c
        return
c
c
c
c
ccc        entry hank106b(r,h0,h1,w,ijw)
        entry hank106b(r,h0,h1,w,ijw,ninterv)
c
        call hank106eva(r,w(ix1),n,w(icenters),w(ijw(1)),
     1      w(ijw(2)),w(ijw(4)),h0,h1,w(ih),w(ijw(3)),
ccc     2      w(iu0),w(irk) )
     2      w(iu0),w(irk),ninterv)
c
        return
        end
c
c
c
c
c
        subroutine hank106eva(r,x1,n,centers,h0derss,
ccc     1      h0s,h1s,h0,h1,h,h1derss,u0,rk)
     1      h0s,h1s,h0,h1,h,h1derss,u0,rk,ninterv)
        implicit real *8 (a-h,o-z)
        complex *16 h0s(1),h0derss(n,1),h0,h1,
     1      h1derss(n,1),h1s(1),u0,rk,zcom,centers(1),zh,z
c
        save
c 
c   input:
c   ninterv - number of intervals used in equispaced
c       subdivision of current interval
c
c-----------------------------------
c
c       find the subinterval where the point z lives
c
        z=rk*r
c
ccccc        d = (z-x1)*h
ccc        if (d.lt.0) d = 1.0d-12
ccccc       if (d.lt.0) d = 0
ccccc        i = d+1
c
        i=(z-x1)*h +1
ccc        call prinf(' i = *',i,1)
        if (i.lt.0) then 
            i = 1
        else if (i.gt.ninterv) then 
            i = ninterv
        endif
c
c        evaluate the functions h0 and h1 at the point z
c
        zh=z-centers(i)
        t=zh*u0
ccc        call prin2(' zh = *',zh,2)
ccc        call prin2(' t = *',t,1)
c
        h0=(((((((((h0derss(10,i)*t+h0derss(9,i))*t+h0derss(8,i) ) 
     1    * t+h0derss(7,i))*t+h0derss(6,i))*t+h0derss(5,i))*t 
     2    +h0derss(4,i))*t+h0derss(3,i))*t+h0derss(2,i))
     3    *t+h0derss(1,i) ) * t + h0s(i)
c
        h1=(((((((((h1derss(10,i)*t+h1derss(9,i))*t+h1derss(8,i) ) 
     a    * t+h1derss(7,i))*t+h1derss(6,i))*t+h1derss(5,i))*t
     2    +h1derss(4,i) ) *t + h1derss(3,i) )* t +h1derss(2,i))
     3    *t+h1derss(1,i) ) * t +h1s(i)
c
        return
        end
c
c
c
c
c
        subroutine hank106ini(coef,x1,x2,ninterv,n,
     1      centers,h0derss,h0s,h1s,h,h1derss,u0)
        implicit real *8 (a-h,o-z)
        complex *16 h0s(1),h0derss(n,1),h1derss(n,1),ima,u0,
     1      us(22),h1s(1)
        dimension centers(2,1)
c
        data ima/(0.0d0,1.0d0)/
c 
c        construct the subintervals
c
        h=(x2-x1)/ninterv
c
        do 1200 i=1,ninterv
c
        ab1i=(i-1)*h+x1
        ab2i=(i-1)*h+x1 +h
        centers(1,i)=(ab2i+ab1i)/2
        centers(2,i)=coef*centers(1,i)
 1200 continue
c
c        construct the values of Hankel functions and their
c        derivatibes at the centers
c
        do 1400 i=1,ninterv
c
        call hank0ders(centers(1,i),n,h0s(i),h1s(i),
     1      h0derss(1,i),h1derss(1,i) )
 1400 continue
c
c       scale them things by factorials and by complex powers
c
        u0=1+ima*coef
        u0=u0/abs(u0)
c
        us(1)=u0
        do 1500 i=1,20
c
        us(i+1)=us(i)*u0
 1500 continue
c
        do 1800 i=1,ninterv
        fact=1
        do 1600 j=1,n-1
        h0derss(j,i)=h0derss(j,i)*fact * us(j)
        h1derss(j,i)=h1derss(j,i)*fact * us(j)
c
        fact=fact/(j+1)
 1600 continue
c
 1800  continue
c
        return
        end
c
c
c
c
c
        subroutine hank0ders(z,n,h0,h1,h0ders,h1ders)
        implicit real *8 (a-h,o-z)
        complex *16 z,h0,h0ders(1),h1,h1ders(1)
c
        data ifexpon/1/
c 
c        evaluate h0 and h1
c
        call hank103(z,h0,h1,ifexpon)
c
        h0ders(1)=-h1
        h0ders(2)=-(h0ders(1)/z+h0)
        h0ders(3)=-(2*h0ders(2)+h0ders(1)*z+h0)/z
        h0ders(4)=-(3*h0ders(3)+h0ders(2)*z+2*h0ders(1))/z
c
        if(n .le. 4) return
c
        do 1400 m=2,n-2
c
        h0ders(m+2)=-( (m+1)*h0ders(m+1)+z*h0ders(m)+
     1      m*h0ders(m-1) )/z
 1400 continue
c
        do 1600 i=1,n-1
c
        h1ders(i)=-h0ders(i+1)
 1600 continue
c
        return
        end


ccc        subroutine hank106datagen(rk,rmin,rmax,ab,nab,ninterval,
ccc     1             w,lw,istart,ier)
        subroutine hank106datagen(rk,ier)
        implicit real *8 (a-h,o-z)
ccc        integer *4 istart(1)
        integer *4 istart(29),nintervec(28)
ccc        dimension w(1),ab(2,nab)
        dimension w(50000),ab(2,28)
        complex *16 rk,h0,h1,z,rksav
        data nab/28/
        data lw/50000/
        data rmin/1.0d-6/
        data rmax/200/
        save nab,ninterval,lw,istart,nintervec
        save w,ab,rmin,rmax,rksav
c
c     INPUT PARAMETERS:    ----> Now hidden as local vars.....
c
c       create top level (dyadic) intervals for hank106init, which 
c       then uses equisized subintervals to precompute interpolation
c       polynomials
c
c     rk (complex *16) frequency parameter
c
c     rmin, rmax (real *8)  desired range of argument to hank106
c                           [rmin*rk,...,rmax*rk]
c 
c     ab(2,nab) (real *8)   blank array of length 2*nab
c     w(lw)     (real *8)   work array of length lw
c
c     OUTPUT PARAMETERS:
c
c     ninterval (integer *4)   number of subintervals created
c     ab(2,ninterval) (real *8) boundary of ith interval is
c                               (ab(1,i),ab(2,i))
c     nintervec(ninterval)    nuomber of equispaced subintervals
c                             used for ith interval
c     istart  (integer *4)     istart(i) is pointer into workspace for
c                              data pertaining to ith interval
c     ier (integer *4) error flag
c             ier = 0 upon normal execution.
c             ier = 1 if length (nab) of array ab is of 
c                     insufficient length 
c             ier = 2 if length (lw) of workspace w is of 
c                     insufficient length 
c-----------------------------------------------------------------
c
c
        ier = 0
        ninterval = 1
        istart(1) = 1
        rmaxloc = rmin
ccc        rksav = rk
        rksav = rk/cdabs(rk)
        rminsav = rmin
        do i = 1,1000
           rminloc = rmaxloc 
           rmaxloc = 2*rminloc 
           if (rminloc.gt.100) rmaxloc = rminloc+100
           ab(1,i) = rminloc
           ab(2,i) = rmaxloc
ccc           call prinf(' i = *',i,1)
ccc           call prin2(' rminloc = *',rminloc,1)
ccc           call prin2(' rmaxloc = *',rmaxloc,1)
ccc           call prinf(' istart(i) = *',istart(i),1)
ccc           call hank106init(rk,rminloc,rmaxloc,w(istart(i)),keep)
ccc           call hank106init(rksav,rminloc,rmaxloc,w(istart(i)),keep)
           call hank106init(rksav,rminloc,rmaxloc,w(istart(i)),
     1        keep,ninterv)
           nintervec(i) = ninterv
ccc           call prinf(' keep = *',keep,1)
           istart(i+1) = istart(i) + keep + 1
           if (rmaxloc.ge.rmax) goto 1111
           if (i.gt.nab) then
              ier = 1
              return
           endif
           if (istart(i+1).gt.lw) then
              ier = 2
              return
           endif
           ninterval = ninterval + 1
        enddo
1111    continue
ccc        call prinf(' istart(i+1) = *',istart(i+1),1)
ccc        call prinf(' i+1 = *',i+1,1)
ccc        call prinf(' ninterval = *',ninterval,1)
ccc        call prinf(' istart(ninterval) = *',istart(ninterval),1)
        return
c
ccc        entry hank106(z,h0,h1,ab,ninterval,w,istart)
        entry hank106(z,h0,h1,ifexpon)
c
c
c       z   (complex *16) argument for Hankel function evaluation.
c       h0,h1 (complex *16) H_0(z) and H_1(z) where frequency 
c                           parameter is ASSUMED TO BE that from 
c                           previous call to hank106datagen.
c
c       ab,ninterval,w,istart   defined above.
c
c-------------------------------------------------------------
c
c       determine subinterval and call hank106.
c
ccc        call prinf(' ninterval = *',ninterval,1)
ccc        call prin2(' rminsav = *',rminsav,1)
ccc        call prin2(' rksav = *',rksav,2)
ccc        x = dreal(z/(rksav*rminsav))
ccc        call prin2(' x = *',x,1)
ccc        call prinf(' ifexpon = *',ifexpon,1)
        if (ifexpon.eq.0) then
           call hank103(z,h0,h1,ifexpon)
           return
        endif
        r = dreal(z/rksav)
        call findinte(r,ab,ninterval,i)
ccc        call prin2(' z = *',z,2)
cc        call prin2(' rksav = *',rksav,2)
ccc        call prin2(' r = *',r,1)
ccc        call prinf(' i = *',i,1)
        if (i.le.ninterval) then
cc        call prin2(' w is = *',w(istart(i)),10)
cc        call prinf(' int part is = *',w(istart(i+10)),10)
ccc           call hank106a(r,h0,h1,w(istart(i)))
           call hank106a(r,h0,h1,w(istart(i)),nintervec(i))
ccc        call prin2(' 106 gives h0 is = *',h0,2)
ccc        call hank103(z,h0,h1,ifexpon)
ccc        call prin2(' 103 gives h0 is = *',h0,2)
        else
           call hank103(z,h0,h1,ifexpon)
        endif
        return
        end
c
c
c
        subroutine findinte(x,ab,nn,intnum)
        implicit real *8 (a-h,o-z)
        integer *4 intold,ithresh
        dimension ab(2,nn)
c
        data intold/-10/
        data ithresh/10/
c
c       check if the point is on the subinterval as the preceding one
c
        if(intold .le. 0) goto 2000
        if(intold .gt. nn) goto 2000
c
        intnum=intold
ccc        call prinf(' intnum is *',intnum,1)
ccc        call prin2(' ab is *',ab(1,intnum),2)
        if( (x .ge. ab(1,intnum) ) .and. (x .le. ab(2,intnum) ) ) return
c
 2000 continue
       if(x .lt. ab(1,1)) then
           intnum = 777
           return
       else if(x .gt. ab(2,nn)) then
           intnum = 777
           return
       endif
c
c      the point is not on the same subinterval as the preceding one.
c      if nn is less than ithresh, use direct scan to find the proper 
c      interval
c
       if(nn .gt. ithresh) goto 3000
c
c
        do 2200 j=1,nn
c
           intnum=j
c
        if(ab(2,j) .ge. x) goto 2400
 2200 continue
c
 2400 continue
c
        intold=intnum
        return
c
 3000 continue
c
c      The point is not on the same subinterval as the preceding one,
c      and nn is greater than ithresh; use bisection to find the proper 
c      interval
c
       i1=1
       i2=nn
       i3=(i1+i2)/2
c
cccc       nsteps=0
       do 3400 i=1,100
c
       if(x .ge. ab(1,i3)) i1=i3
       if(x .le. ab(2,i3)) i2=i3
c
       if(i2 .eq. i1) goto 3600
c
       i3=(i1+i2)/2
 3400 continue
c
 3600 continue

       if(x .lt. ab(1,i3)) i3=i3-1
       if(x .gt. ab(2,i3)) i3=i3+1

       intnum=i3
       intold=intnum
c       
        return
        end
c
