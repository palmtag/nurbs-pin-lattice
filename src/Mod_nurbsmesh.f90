   module mod_nurbsmesh
   implicit none
!-----------------------------------------------------------------------
!
!  Module to store NURBS mesh information
!
!  Scott Palmtag
!  August 2026
!
! @version CVS $Id: Mod_nurbsmesh.f90,v 1.2 2026/08/22 03:01:37 palmtag Exp $
!
!-----------------------------------------------------------------------

!--- mesh data

      integer, parameter :: iorder=4       ! "square" elements

      integer, protected :: maxvert        ! maximum sizes
      integer, protected :: maxedge
      integer, protected :: maxelem

      integer :: nvert    ! number of vertices
      integer :: nedge    ! number of edges
      integer :: nelem    ! number of elements
      integer :: nbound   ! number of boundary elements
      integer :: nknot    ! number of knot vectors

      real(8), allocatable :: xi(:,:)      ! (3,nvert)       x,y,z for each node
      real(8), allocatable :: xielem(:,:)  ! (2,nelem)   *** allocated in code  x,y,z of element 3D
      real(8), allocatable :: ximid(:,:)   ! (2,nedge)       x,y,z of midpoint
      real(8), allocatable :: weight(:)    ! (nvert)
      real(8), allocatable :: wmid(:)      ! (nedge)         edge midpoint weights
      integer, allocatable :: gg(:,:)      ! (iorder,nelem)  vertice numbers per element
      integer, allocatable :: matl(:)      ! (nelem)     material numbers / element attributes
      integer, allocatable :: gedge(:,:)   ! (3,nedge)   knot, vert1, vert2

   contains

!=======================================================================
      subroutine init_nurbsmesh(nrow)
      implicit none
      integer, intent(in) :: nrow

!--- max values  (may not be enough if using large number of rings)

      if (nrow.le.3) then
         maxvert=1000
         maxedge=1000
         maxelem=1000
      else
         maxvert=40*nrow*nrow
         maxedge=70*nrow*nrow
         maxelem=40*nrow*nrow
      endif

      write (*,*) 'allocating NURBS mesh'
      write (*,*) '  maxvert  ', maxvert
      write (*,*) '  maxedge  ', maxedge
      write (*,*) '  maxelem  ', maxelem

!--- initialize

      nvert=0
      nedge=0
      nelem=0    ! 2 rings of 4 + center
      nbound=0

      nknot=4    ! number of knot vectors

!--- allocate mesh data

      allocate(gg(iorder,maxelem))   ! 4 vertices
      allocate(matl(maxelem))        ! materials / attributes
      gg=999999      ! use big number to catch errors easily
      matl=999999

      allocate(gedge(3,maxedge))   ! knot, vert1, vert2
      gedge=999999

      allocate (weight(maxvert))
      allocate (xi(3,maxvert))
      weight=1.0d0
      xi=1.0d20

      allocate (wmid(maxedge))
      allocate (ximid(2,maxedge))
      wmid=1.0d0
      ximid=1.0d20

      return
      end subroutine init_nurbsmesh
!=======================================================================

   end module
