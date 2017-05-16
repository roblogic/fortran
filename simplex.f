      program main

      implicit none                ! Set up variables..
      integer o,p
      parameter(o=20,p=100)
      integer i,j,m,n,r,s
      integer apos(o),bpos(o)
      real A(o,p),Bi(o,o),b(o),xb(o),cost(p),ca(o+p)
                                   ! Read data..
      open(unit=8,file='lp.dat',status='old')
      open(unit=9,file='lp.out',status='new')
      read(8,*) m,n
      if (m.gt.o.or.n.gt.p) then
        write(9,*)'Exceeded variable/constraint limits'
      endif
      read(8,*) (cost(j),j=1,n)
      write(9,*) 'm    ',m
      write(9,*) 'n    ',n
      write(9,*) 'c''  ',(cost(j),j=1,n)
      write(9,*)
      write(9,*) 'A.x = b'
      do i=1,m
        read(8,*) (A(i,j),j=1,n),b(i)
        write(9,*) (A(i,j),j=1,n),' =',b(i)
      enddo
c ------------------extra bit..-----------------
      do i=1,m
        read(8,*) (Bi(i,j),j=1,m)
        write(9,*) (Bi(i,j),j=1,m)
      enddo
      read(8,*) (bpos(i),i=1,m)
      write(9,*) (bpos(i),i=1,m)
      read(8,*) (apos(i),i=1,n)
      write(9,*) (apos(i),i=1,n)
      do i=1,m
        do j=1,m
          xb(i)=xb(i)+Bi(i,j)*b(j)
        enddo
      enddo
      goto 5
c -----------------------miss out for now------------------------
      write(9,*) 
      write(9,*) 'Phase 1'
      do i=1,n                     ! Phase 1
        ca(i)=0                    !  Prepare artificial cost,basis,variables 
        apos(i)=0                  !  all original variables are nonbasic
      enddo
      do i=1,m 
        do j=1,m
          Bi(i,j)=0                !  initially Bi = identity matrix, I
        enddo
        Bi(i,i)=1
        ca(n+i)=1                  !  artificial cost vector now constructed
        bpos(i)=n+i                !  artificial variables fill the basis
        xb(i)=b(i)                 !  initially xb = b (because xb=Bi.b, Bi=I)
      enddo
      call rsm(m,n,A,Bi,b,xb,ca,bpos,apos,1)
c ---------------------------------------------------------------
5     write(9,*)
      write(9,*) 'Phase 2'
                                   ! Phase 2
                                   !  Now do RSM with the new basis
      call rsm(m,n,A,Bi,b,xb,cost,bpos,apos,2)
      end
C--------------------------------------------------

      subroutine rsm(m,n,A,Bi,b,xb,c,bpos,apos,phase)
      parameter(o=20,p=100)
      integer apos(n),bpos(m),phase      ! variables passed to rsm
      real A(o,p),Bi(o,o),b(m),c(n+m)    !     "
      integer i,j,imax,r,s,note          ! introducing new variables
      real bestrc,bias(m),pi(m),xb(m),z  !     "
      character retval

      imax=2*(m+n)
      tiny=1E-10
      i=0
      call currentsoln(Bi,b,bpos,c,m,xb,z)
1     if (i.lt.imax) then                       ! Iterate Revised Simplex Method..
        i=i+1
        call computedual(Bi,c,bpos,m,n,pi)      ! find pi (initially pi=e)
        call reducedcost(A,c,pi,apos,m,n,s,bestrc)
        if(s.eq.0) then
          if(z.gt.tiny.and.phase.eq.1) then
            write(9,*) '  Infeasibility detected.'
            note=-2
          else
            write(9,*) '  Optimality reached.'
            note=1
          endif
          goto 3
        endif
        call mivitwwbias(A,Bi,m,s,bias)
        call lvratiotest(xb,bias,m,n,r,phase,bpos)
        if (r.eq.0) then
          write(9,*) '  Unboundedness detected.'
          note=-1
          goto 3
        endif
      write(9,*) 'basisupdate:'
      write(9,*) '  bpos',(bpos(j),j=1,m),'   apos',(apos(j),j=1,n)
        call basisupdate(Bi,bias,c,pi,xb,bpos,apos,m,n,r,s)
      write(9,*) '    ->',(bpos(j),j=1,m),'     ->',(apos(j),j=1,n)
        call currentsoln(Bi,b,bpos,c,m,xb,z)
       else
        goto 3
      endif
2     goto 1
3     write(9,*) '  RSM iterations executed: ',i
      end
C-----------------------------------------------------

      subroutine currentsoln(Bi,b,bpos,c,m,xb,z)
      parameter(o=20)
      integer bpos(*)                     ! old
      real Bi(o,o),b(*),c(*),xb(*),z      ! old
      integer i,j                         ! new
      z=0
      do i=1,m
        z=z+c(bpos(i))*xb(i)
      enddo     

      write(9,*) 'currentsoln:  Bi, xb'
      do i=1,m
      write(9,*) ' ',(Bi(i,j),j=1,m),'     x',bpos(i),' =',xb(i)
      enddo
      write(9,*)
      write(9,*) '  --------------------------- Objective value z =',z
      write(9,*)
      end
C-----------------------------------------------------

      subroutine computedual(Bi,c,bpos,m,n,pi)
      parameter(o=20)
      integer bpos(*)                     ! old
      real Bi(o,o),c(*),pi(*)             ! old
      integer i,j,k                       ! new
      do i=1,m
        pi(i)=0
        do j=1,m
          k=bpos(j)
          pi(i)=pi(i)+c(k)*Bi(j,i)
        enddo
      enddo
      write(9,*) 'computedual:    pi ',(pi(i),i=1,m)
      end
C----------------------------------------------------

      subroutine reducedcost(A,c,pi,apos,m,n,s,bestrc)
      parameter(o=20,p=100)
      integer apos(*),s                   ! old
      real A(o,p),c(*),pi(*),bestrc       ! old
      integer i,j                         ! new
      real rc,pian,tiny                   ! new
      rc=0
      bestrc=0
      tiny=-1E-10
      do i=1,n                     !  determine entering var x(s)
        if(apos(i).eq.0) then      !< look at nonbasic columns
          pian=0                   !  (never price artificials)
          do j=1,m
            pian=pian+pi(j)*A(j,i)
          enddo
          rc=c(i)-pian             !< red. cost for nonbasic x(i)
          if(rc.lt.bestrc) then
            bestrc=rc
            s=i
          endif
        endif
      enddo
      if(bestrc.ge.tiny) then      !< This means no ev was found,
        s=0                        !  optimality has been reached,
      endif                        !  so flag s to 0
      write(9,*) 'reducedcost:    x',s,' enters with rc',bestrc
      end
C----------------------------------------------------

      subroutine mivitwwbias(A,Bi,m,s,bias)
      parameter(o=20,p=100)
      real A(o,p),Bi(o,o),bias(*)
      integer i,j,s
      do i=1,m                     !  a quick calculation of bias,
        bias(i)=0                  !  aka. the most important vector
        do j=1,m                   !  in the whole world.
          bias(i)=bias(i)+Bi(i,j)*A(j,s)
          enddo
        enddo
      write(9,*) 'mivitwwbias:    bias',(bias(j),j=1,m)
      end
C----------------------------------------------------

      subroutine lvratiotest(xb,bias,m,n,r,phase,bpos)
      integer bpos(*),r,phase      ! old
      real xb(*),bias(*)           ! old
      integer i                    ! new
      real ratio,rmin,tiny         ! new
      r=0                          ! r=0 returned if no lv can be found
      tiny=1E-10
      teeny=-1E-10
      rmin=10E+10
      do i=1,m                    
        if(bpos(i).gt.n.and.phase.eq.2) then  
                                   ! Artificial present in phase2 basis: 
                                   ! Extended lv routine to force out
                                   !   artificial variables..
          if(bias(i).lt.teeny.or.bias(i).gt.tiny) then
            rmin=0
            r=i
          endif
        else                       ! Usual lv routine..
          if(bias(i).gt.tiny) then !   ratio = rate of change of xb(i), as
            ratio=xb(i)/bias(i)    !           xs increases from 0 (enters)
            if(ratio.lt.rmin) then 
              rmin=ratio           !   the leaving variable is xb(r), the
              r=i                  !   basic variable which decreases most
            endif                  !   rapidly as xs enters.
          endif
        endif
      enddo
      write(9,*) 'lvratiotest:    xb',r,' leaves.'
      end
C----------------------------------------------------

      subroutine basisupdate(Bi,bias,c,pi,xb,bpos,apos,m,n,r,s)
      parameter(o=20)
      integer apos(*),bpos(*),m,n,r,s        ! old
      real Bi(o,o),bias(*),c(*),pi(*),xb(*)  ! old
      integer i,j                            ! new
      do i=1,m                     ! do gauss-jordan pivot on rth element
        if(i.eq.r) then            ! of bias, and thus also on Bi and xb
          do j=1,m
            Bi(r,j)=Bi(r,j)/bias(r)
            enddo
          xb(r)=xb(r)/bias(r)
        else
          do j=1,m
            Bi(i,j)=Bi(i,j)-bias(i)/bias(r)*Bi(r,j)
            enddo
          xb(i)=xb(i)-bias(i)/bias(r)*xb(r)
        endif
      enddo
      if(bpos(r).le.n) then
        apos(bpos(r))=0              ! xb(r) now nonbasic; 0 in apos
      endif
      bpos(r)=s                      ! xb(r) replaced by xs in bpos
      apos(s)=r                      ! xs located in rth basic position 
      end
