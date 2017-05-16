      program minimax
! created by   R. Papesch 
! created on   09 Sep 2009
! license      GNU GPLV2.0, please see 
!              https://github.com/papesch/fortran/blob/master/LICENSE
! description  Sets up an LP file to be read by sim.g and minimised..
      implicit none
      integer j,m,n,p
      double precision x(0:40),b(40),i,k,l
      character*10 flnm

      print *,'Minimax Polynomial interpolation'
      print *
      print *,'Enter output filename:'
      read *,flnm
      print *,'Enter p:'
      read *,p
      print *,'Enter l:'
      read *,l
      m=p+2
      n=2*l
      k=l-1.0
      x(0)=0
      do i=1,k
        x(i)=2*i/k
        b(i)=exp(-2*x(i))
      enddo
      
      open(unit=8,file=flnm,status='new')
      write(8,*) m
      write(8,*) n
      write(8,*) (-b(i),i=0,k),(b(i),i=0,k)  ! cost vector
      write(8,*) (1,i=1,n),1                 ! 1st constraint, e vectors, rhs=1
      write(8,*) (1,i=0,k),(-1,i=0,k),0      ! 2nd constraint, +/- x(i)^0, rhs=0
      do j=1,p                               ! iterate thru other constraints
        write(8,*) (x(i)**j,i=0,k),(-x(i)**j,i=0,k),0,0
      enddo
      end
