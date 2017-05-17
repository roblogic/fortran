      program main
! created by   R. Papesch 
! created on   09 Sep 1998 or earlier
! license      GNU GPLV2.0, please see 
!              https://github.com/papesch/fortran/blob/master/LICENSE
! description  Solves LP problem with Revised Simplex Method. 
!              (This code is probably superseded by simplex.f)

      implicit none
      integer o,p
      parameter(o=20,p=100)
      integer i,j,m,n,r,s
      integer apos(o),bpos(o)
      double precision A(o,p),Bi(o,o),b(o),cost(p),ca(o+p),xb(o)

      open(unit=8,file='lp.dat',status='old')
      read(8,*) m,n
      if (m.gt.o.or.n.gt.p) then
        write(*,*) 'Exceeded variable/constraint limits'
      endif
      read(8,*) (cost(j),j=1,n)
      do i=1,m
        read(8,*) (A(i,j),j=1,n),b(i)
      enddo

      write(*,*) 'Commencing Phase 1..'
      do i=1,n                     ! Prepare artificial cost,basis,variables
        ca(i)=0                    ! for phase 1..
        apos(i)=0                  !  all variables are nonbasic
      enddo
      do i=1,m 
        do j=1,m
          Bi(i,j)=0                !  initial Bi=identity matrix
        enddo
        Bi(i,i)=1
        ca(n+i)=1                  !  artificial cost vector now constructed
        bpos(i)=n+i                !  artificials fill the basis
      enddo
      call rsm(m,n,A,Bi,b,ca,bpos,apos,1)

                                   ! Now do RSM with the new basis:Phase 2
      write(*,*) 'Commencing Phase 2..'
      call rsm(m,n,A,Bi,b,cost,bpos,apos,2)
      end

C--------------------------------------------------

      subroutine rsm(m,n,A,Bi,b,c,bpos,apos,phase)
      integer i,icount,imax,r,s,note
      double precision z
      do i=1,imax                  ! iterate until optimality or failure
        icount=i
        call currentsoln(bi,b,bpos,m,xb,z)       ! initially xb=b
        call computedual(Bi,ca,bpos,m,n,pi)      ! initially pi=e
        call reducedcost(A,ca,pi,apos,m,n,s,rc)
        if(s.eq.0) then
          if(z.gt.tiny) then
            write(*,*) 'Infeasibility detected.'
            note=-2
          else
            write(*,*) 'Optimality reached.'
            note=1
          endif
          i=imax
          goto 1
        endif
        call mivitwwbias(A,Bi,m,n,s,bias)
        call lvratiotest(xb,bias,m,r,phase,bpos)
        if (r.eq.0) then
          write(*,*) 'Unboundedness detected.'
          note=-1
          i=imax
          goto 1
        endif
        call basisupdate(Bi,bias,ca,pi,xb,bpos,apos,m,n,r,s)
1    enddo
     write(*,*) 'RSM iteration count: ',icount
     end

C-----------------------------------------------------
 
      subroutine currentsoln(Bi,b,bpos,c,m,xb,z)
      z=0
      do i=1,m
        xb(i)=0
        do j=1,m
          xb(i)=xb(i)+Bi(i,j)*b(j)
        enddo
        write(*,*) 'x',bpos(i),'=',xb(i)
        z=z+c(bpos(i))*xb(i)
      enddo     
      write(*,*) 'Objective value z=',z
      end
C-----------------------------------------------------

      subroutine computedual(Bi,c,bpos,m,n,pi)
      integer i,j,k
      do i=1,m
        k=bpos(i)                    ! k is just an abbreviation
        pi(k)=0
        do j=1,m
          pi(k)=pi(k)+c(k)*Bi(j,k)
        enddo
      enddo
      end
C----------------------------------------------------

      subroutine reducedcost(A,c,pi,apos,m,n,s,bestrc)
      double precision rc,bestrc,pian,tiny
      rc=0
      bestrc=0
      tiny=-10E-10
      do i=1,n                       ! determine entering var x(s)
        if(apos(i).eq.0) then        ! look at nonbasic columns
          pian=0                     
          do j=1,m
            pian=pian+pi(j)*A(j,i)
          enddo
          rc=c(i)-pian               ! red. cost for nonbasic x(i)
          if(rc.lt.bestrc) then
            bestrc=rc
            s=i
          endif
        endif
      enddo
      if(bestrc.ge.tiny) then        ! This means no ev was found,
        s=0                          ! optimality has been reached
      endif                          ! so flag s to 0
      end
C----------------------------------------------------

      subroutine mivitwwbias(A,Bi,m,n,s,bias)
      do i=1,m
        bias(i)=0
        do j=1,m
          bias(i)=bias(i)+Bi(i,j)*A(j,s)
          enddo
        enddo
      end
C----------------------------------------------------

      subroutine lvratiotest(xb,bias,m,r,phase,bpos)
      double precision ratio,rmin,tiny
      r=0                            ! r=0 returned if no lv can be found
      tiny=10E-10
      teeny=-10E-10
      rmin=10E+10
      do i=1,m                    
        if(bpos(i).gt.n.and.phase.eq.2) then  ! artificial in phase2 basis
          if(bias(i).lt.teeny.or.bias(i).gt.tiny) then
            rmin=0                            ! Extended lv routine to force
            r=i                               ! out artificial variables
          endif
        else
          if(bias(i).gt.tiny) then   ! Usual lv routine..
            ratio=xb(i)/bias(i)      !   ratio = rate of change of xb(i), as
            if(ratio.lt.rmin) then   !           xs increases from 0 (enters)
              rmin=ratio
              r=i                    !   the leaving variable is xb(r), the
            endif                    !   basic variable which decreases most
          endif
        endif
      enddo
      end
C----------------------------------------------------

      subroutine basisupdate(Bi,bias,c,pi,x,bpos,apos,m,n,r,s)
      integer i,j
      do i=1,m                       ! do gauss-jordan pivot on rth element
        if(i.eq.r) then              ! of as, and thus also on Bi and xb
          do j=1,m
            Bi(i,j)=Bi(i,j)/bias(r)
            enddo
          xb(i)=xb(i)/bias(r)
        else
          do j=1,m
            Bi(i,j)=Bi(i,j)-bias(i)/bias(r)*Bi(r,j)
            enddo
          xb(i)=xb(i)-bias(i)/bias(r)*xb(i)
        endif
        bias(i)=0
      enddo
      bias(r)=1
      if(bpos(r).le.n) then
        apos(bpos(r))=0    ! xb(r) now nonbasic; 0 in apos
      endif
      bpos(r)=s            ! xb(r) replaced by xs in bpos
      apos(s)=r            ! xs located in rth basic position (in apos)
      end
