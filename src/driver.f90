   program driver
   use mod_input, only : readinput, fbase, nrow, irodmap
   use mod_nurbsmesh, only : init_nurbsmesh
   implicit none
!=======================================================================
!
!  Main driver routine
!
!  Create MFEM NURBS mesh for an array of pincells (PWR lattice)
!
!  Scott Palmtag
!  August 2026
!
! @version CVS $Id: driver.f90,v 1.6 2026/08/22 03:03:50 palmtag Exp $
!
!=======================================================================

      integer :: ia
      integer :: i, j

      integer :: isqbnd(8)
      integer :: isqbnd0(8)

      character(len=200) :: fname    ! input file name

!--- input (temp)

      integer :: itype

!--- read input file name from command line

      ia=iargc()
      if (ia.ne.1) then
        stop 'usage: makenurbs [input]'
      endif
      call getarg(1,fname)

      fbase=fname
      i=len_trim(fbase)    ! remove suffix
      if (i.gt.4) then
        if (fbase(i-3:i).eq.'.inp') fbase(i-3:i)=' '
      endif

!--- read input

      call readinput(fname)

!--- initialize nurbs mesh, allocate memory

      call init_nurbsmesh(nrow)

!--- fill vertices and edges for the square mesh first

   ! todo: flag the boundary edges

      call filledge()

!--- loop over pincells and fill mesh

   ! isqbnd is a map of the vertices surrounding each pincell
   ! it is not very elegant
   ! see top of pincell_octant for numbering convention
   ! we don't actually need all 8 values, only the 3 on the left

      isqbnd0(1)=0    ! lower left
      isqbnd0(2)=isqbnd0(1)+1
      isqbnd0(3)=isqbnd0(1)+2
      isqbnd0(8)=isqbnd0(1)+2*nrow+1
      isqbnd0(4)=isqbnd0(8)+2
      isqbnd0(7)=isqbnd0(8)+2*nrow+1
      isqbnd0(6)=isqbnd0(7)+1
      isqbnd0(5)=isqbnd0(7)+2

      do j=nrow, 1, -1        ! start in lower left corner
        isqbnd=isqbnd0
        do i=1, nrow
          itype=irodmap(i,j)
          call pincell_octant(isqbnd, itype)
          isqbnd=isqbnd+2
        enddo
        isqbnd0=isqbnd0+4*nrow+2   ! increment for new row
      enddo

!--- check that every edge is used (debug)

      call checkedge()

!--- fill 3D information

      call make_xielem()    ! generate 3D element data

!--- write mesh file

      call writemesh(fbase)

      write (*,*) 'done'
      end program driver
