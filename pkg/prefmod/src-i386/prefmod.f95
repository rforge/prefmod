!gfortran  -O3 -Wall  -c gens3.f95 -o gens3.o
!gfortran -shared -s  -o gens3.dll gens3.o

! gfortran --shared -o gens.dll gens.f95
! gfortran maxi.f95 -o maxi.exe


       subroutine calcs(ncomp,notnaidx,iout)

       integer(kind=4) izahl, ncomp, notnaidx
       integer(kind=4), dimension(2**ncomp)::iout

       do izahl=0,2**ncomp-1
         newbyte=0
         j=0
         do i=0,ncomp
           if(IBITS(notnaidx,i,1).eq.1) then
              call MVBITS(izahl,i,1,newbyte,j)
              j=j+1
           end if
         end do
         iout(izahl+1)=newbyte+1
       end do

       end subroutine

!------------------------------------------------------------

       subroutine calcs3(ncomp,notnaidx,iout)

       integer(kind=4) izahl, ncomp, notnaidx
       integer(kind=4), dimension(ncomp)::tern,terns
       integer(kind=4), dimension(3**ncomp)::iout
       integer, dimension(3)::nbd
       data nbd/0,1,2/

       do izahl=0,3**ncomp-1

         ! init
         number=izahl
         do j=1,ncomp
            tern(j)=0
            terns(j)=0
         end do

         ! convert to ternary
         ipos=1
         do while (number>0)
             idx=floor(((number/3.)-floor(number/3.))*3+1.5)
             tern(ncomp+1-ipos)=nbd(idx)
             number=floor(number/3.)
             ipos=ipos+1
         end do


         ! remove NAs
         ishift=0
         do j=0,ncomp-1
           if (IBITS(notnaidx,j,1).eq.1) then
              terns(ncomp-j+ishift)=tern(ncomp-j)
           else
              ishift=ishift+1
           end if
         end do


         ! convert to integer
         idec=0
         do j=0,ncomp-1
           idec=idec+terns(ncomp-j)*3**j
         end do

         iout(izahl+1)=idec+1


       end do

       end subroutine

!------------------------------------------------------------
       subroutine pattp(ncomp, nobj, lambda, ppatt)

       integer(kind=4) ncomp, nobj, row
       integer(kind=4), dimension(ncomp,nobj)::b
       double precision, dimension(nobj-1)   :: lambda
       double precision, dimension(2**ncomp) :: ppatt
       double precision                      :: sumppatt
       integer(kind=4), dimension(ncomp)     :: a,aa,x


       do i=1,ncomp
         do j=1,nobj
           b(i,j)=0
         end do
       end do

       row=1
       do j=2,nobj
         do i=1,j-1
            b(row, i)=1
            b(row, j)=-1
            row=row + 1
         end do
       end do

       sumppatt=0D0
       do izahl=0,2**ncomp-1

         ppatt(izahl+1)=0D0
         do j=1,ncomp
           a(ncomp-j+1) = IBITS(izahl,j-1,1)
         enddo


         do j=1,nobj
           x(j)=0
           do i=1,ncomp
              aa(i)= 1-2*a(i)
              x(j)=x(j)+aa(i)*b(i,j)
           end do
         end do

         do j=1,nobj-1
           ppatt(izahl+1)=ppatt(izahl+1)+x(j)*lambda(j)
         enddo
         sumppatt=sumppatt+dexp(ppatt(izahl+1))
       end do

       do i=1,2**ncomp
          ppatt(i)=dexp(ppatt(i))/sumppatt
       end do
       end subroutine
