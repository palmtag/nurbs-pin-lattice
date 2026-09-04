!=======================================================================
!
!  Mesh utilities for NURBS mesh generator
!
!  Scott Palmtag
!  August 2026
!
! @version CVS $Id: meshutil.f90,v 1.3 2026/08/22 03:01:37 palmtag Exp $
!
!=======================================================================
!
!  Subroutine to fill edge mesh for entire problem
!   (not the interior of the pincells themselves)
!
!=======================================================================
      subroutine filledge()
      use mod_nurbsmesh
      use mod_input, only : nrow, apitch, ppitch
      implicit none

!--- local

      integer :: i, j
      integer :: iknot
      integer :: ioff
      integer :: isize
      integer :: j1, j2

      real(8) :: gap
      real(8) :: yy
      real(8), allocatable :: xval(:)  ! (nrow)
      real(8), allocatable :: xmid(:)  ! (nrow)

!--- fill x-value array to include assembly gap

      write (*,*)
      write (*,*) 'subroutine filledge'

      gap=(apitch-nrow*ppitch)*0.5d0
      write (*,'(1x,a,f12.6)') 'half gap = ', gap
      if (gap.lt.0.0d0 .or. gap.gt.0.2) then
         write (*,*) 'apitch = ', apitch
         write (*,*) 'ppitch = ', ppitch
         write (*,*) 'nrow   = ', nrow
         write (*,*) 'product= ', ppitch*dble(nrow)
         write (*,*) '***** invalid gap ******'
         stop 'invalid gap'
      endif

!  xval is the edge of the squares
!  xmid is the midpoint of the squares
!    need to correctly account for gaps on first and last cell

      allocate(xval(nrow))   ! xval(0)=0.0
      do i=1, nrow
        xval(i)=dble(i)*ppitch+gap
      enddo
      xval(nrow)=apitch

      allocate(xmid(nrow))
      do i=1, nrow
        xmid(i)=xval(i)-0.5d0*ppitch
      enddo
      if (nrow.gt.1) then
        i=nrow
        xmid(i)=xval(i-1)+0.5d0*ppitch  ! account for gap on last cell
      endif

!--- fill vertices and edges for the square problem first

!  start on bottom row and number up and to the right

      j=0   ! bottom row
        yy=0.0d0   ! bottom row
        i=0
        nvert=nvert+1
        xi(1,nvert)=0.0d0
        xi(2,nvert)=yy
        do i=1, nrow
          nvert=nvert+1
          xi(1,nvert)=xmid(i)   ! half-way in x-direction
          xi(2,nvert)=yy
          nvert=nvert+1
          xi(1,nvert)=xval(i)
          xi(2,nvert)=yy
        enddo

      do j=1, nrow

        yy=xmid(j)  ! midpoint row in y-direction
        i=0
        nvert=nvert+1
        xi(1,nvert)=0.0d0
        xi(2,nvert)=yy
        do i=1, nrow
          nvert=nvert+1
          xi(1,nvert)=xmid(i)   ! half-way in x-direction
          xi(2,nvert)=yy
          nvert=nvert+1
          xi(1,nvert)=xval(i)
          xi(2,nvert)=yy
        enddo

        yy=xval(j)    ! full row
        i=0
        nvert=nvert+1
        xi(1,nvert)=0.0d0
        xi(2,nvert)=yy
        do i=1, nrow
          nvert=nvert+1
          xi(1,nvert)=xmid(i)   ! half-way in x-direction
          xi(2,nvert)=yy
          nvert=nvert+1
          xi(1,nvert)=xval(i)
          xi(2,nvert)=yy
        enddo

      enddo
      write (*,*) 'final nvert ', nvert

      deallocate (xmid)
      deallocate (xval)

!--- fill outer edges
!      start along bottom, work up on right side, then across top, then up on left side

      isize=2*nrow+1

      iknot=1   ! horizontal across bottom
      ioff=0
      do j=0, nrow
        iknot=1   ! horizontal across top of pincells
        do i=1, nrow     ! 2 edges per pincell with octant splitting
          nedge=nedge+1
          gedge(1,nedge)= iknot
          gedge(2,nedge)= 2*i-2+ioff
          gedge(3,nedge)= 2*i-1+ioff
          nedge=nedge+1
          gedge(1,nedge)= iknot
          gedge(2,nedge)= 2*i-1+ioff
          gedge(3,nedge)= 2*i  +ioff
        enddo
        ioff=ioff+2*isize
      enddo

      iknot=0   ! vertical up left side
      ioff=0
      do i=0, nrow
        do j=1, nrow     ! 2 edges per pincell with octant splitting
          nedge=nedge+1
          gedge(1,nedge)= iknot
          gedge(2,nedge)= (2*j-2)*isize +2*i
          gedge(3,nedge)= (2*j-1)*isize +2*i
          nedge=nedge+1
          gedge(1,nedge)= iknot
          gedge(2,nedge)= (2*j-1)*isize +2*i
          gedge(3,nedge)=     2*j*isize +2*i
        enddo
      enddo

      write (*,*) 'final nedge ', nedge

      write (*,*)
      write (*,*) 'debug: vertices for square mesh'
      do i=1, nvert
        write (*,'(i4,2f10.6)') i-1, xi(1:2,i)
      enddo
      write (*,*)
      write (*,*) 'debug: edges for square mesh'
      do i=1, nedge
        write (*,'(4i4)') i, gedge(:,i)
      enddo

!--- define edge midpoints and weights
!    for flat edges, midpoint is simple average
!    weights of 1 have already been assigned

      do i=1, nedge
        j1=gedge(2,i)+1   ! edges are zero based numbering
        j2=gedge(3,i)+1   ! edges are zero based numbering
        ximid(1,i)=0.5*(xi(1,j1)+xi(1,j2))  ! midpoint
        ximid(2,i)=0.5*(xi(2,j1)+xi(2,j2))
      enddo

!--- done

      return
      end subroutine filledge

!=======================================================================
      subroutine make_xielem
      use mod_nurbsmesh
      implicit none

! add 3D element data - take average for now ??

!  generate vertice at center of element
!  I do not think this is used in 2D

!  for now, just take the midpoint of the four vertices of the node

!  todo:  not really sure what the 3D data is for?

      integer :: j, k
      integer :: v1, v2, v3, v4
      real(8) :: wvert1, wvert2, wvert3, wvert4
      real(8) :: wmid1, wmid2, wmid3, wmid4
      real(8) :: xmid1, xmid2, xmid3, xmid4
      real(8) :: ymid1, ymid2, ymid3, ymid4

      allocate (xielem(2,nelem))
      allocate (wcenter(nelem))

      do j=1, nelem

	! Extract 0-based vertex IDs for the element
        v1=gg(1,j)
        v2=gg(2,j)
        v3=gg(3,j)
        v4=gg(4,j)

        ! Initialize mid-weights to 1.0
        wmid1 = 1.0d0; wmid2 = 1.0d0; wmid3 = 1.0d0; wmid4 = 1.0d0

        ! Initialize mid-pints to 0.0
        xmid1 = 0.0d0; xmid2 = 0.0d0; xmid3 = 0.0d0; xmid4 = 0.0d0;
        ymid1 = 0.0d0; ymid2 = 0.0d0; ymid3 = 0.0d0; ymid4 = 0.0d0;

        ! Edge 1: v1 to v2 (bottom)
        k = local2global(j, 1)
        wmid1 = wmid(k)
        xmid1 = ximid(1,k)
        ymid1 = ximid(2,k)

        ! Edge 2: v2 to v3 (right)
        k = local2global(j, 2)
        wmid2 = wmid(k)
        xmid2 = ximid(1,k)
        ymid2 = ximid(2,k)

        ! Edge 3: v3 to v4 (top)
        k = local2global(j, 3)
        wmid3 = wmid(k)
        xmid3 = ximid(1,k)
        ymid3 = ximid(2,k)

        ! Edge 4: v4 to v1 (left)
        k = local2global(j, 4)
        wmid4 = wmid(k)
        xmid4 = ximid(1,k)
        ymid4 = ximid(2,k)

        ! Extract vertex weights (arrays are 1-based, so add 1)
        wvert1 = weight(v1+1)
        wvert2 = weight(v2+1)
        wvert3 = weight(v3+1)
        wvert4 = weight(v4+1)

        ! Apply Transfinite Interpolation formula
        wcenter(j) = 0.5d0*(wmid1 + wmid2 + wmid3 + wmid4) - &
              0.25d0*(wvert1 + wvert2 + wvert3 + wvert4)

        xielem(1,j) = (0.5d0*(xmid1*wmid1 + xmid2*wmid2 + xmid3*wmid3 + xmid4*wmid4) - &
              0.25d0*(xi(1,v1+1)*wvert1 + xi(1,v2+1)*wvert2 + xi(1,v3+1)*wvert3 + xi(1,v4+1)*wvert4))/wcenter(j)
        xielem(2,j) = (0.5d0*(ymid1*wmid1 + ymid2*wmid2 + ymid3*wmid3 + ymid4*wmid4) - &
              0.25d0*(xi(2,v1+1)*wvert1 + xi(2,v2+1)*wvert2 + xi(2,v3+1)*wvert3 + xi(2,v4+1)*wvert4))/wcenter(j)

      enddo

! weight is printed out explicitly in writemesh array

      return
      end subroutine make_xielem

!=======================================================================
      subroutine writemesh(fbase)
      use mod_nurbsmesh
      use mod_input, only : title, maxmatl, matname, matcool
      implicit none

      character(len=*), intent(in) :: fbase

!--- data

      integer :: i, j
      integer :: iatt
      integer, allocatable :: iused(:)
      !real(8) :: wcon1

!--- open file

      write (*,'(3a)') 'creating mfem file: ', trim(fbase), '.mesh'
      open (12,file=trim(fbase)//'.mesh')
      write (12,'(a)') 'MFEM NURBS mesh v1.0'
      write (12,*)
      write (12,'(2a)') '# ', trim(title)
      write (12,'(a)') '#'
      do i=1, maxmatl
        if (matname(i).ne.' ') then
          write (12,24) i, trim(matname(i))
        endif
      enddo
      write (12,24) matcool,'coolant'
   24 format ('#  material ',i0,1x,a)
      write (12,'(a)') '#'

!--- dimensions

      write (12,'(/,a,/,i0)') 'dimension', 2

!--- elements

   ! first number is attribute.
   ! allows you to assign attributes to mesh.
   ! assign material number as attribute

   ! second number is type.  3=square

      write (12,'(/,a,/,i0)') 'elements', nelem
      do i=1, nelem
        write (12,240) matl(i), 3, gg(:,i)
      enddo
 240 format (6(1x,i0))
!--- print boundary
!     boundary attribute, geom, vert1, vert2
!     geometry=1 for segment

      ! Find exterior edges (edges used by exactly 1 element)
      allocate(iused(nedge))
      iused = 0
      do i=1, nelem
        do j=1, 4
          iused(local2global(i,j)) = iused(local2global(i,j)) + 1
        enddo
      enddo

      ! Count actual boundary edges
      nbound = 0
      do i=1, nedge
        if (iused(i) == 1) nbound = nbound + 1
      enddo

      iatt = 3 ! boundary attribute, assume it's reflecting

      write (12,'(/,a,/,i0)') 'boundary', nbound
      do i=1, nedge
        if (iused(i) == 1) then
          write (12,240) iatt, 1, gedge(2,i), gedge(3,i)
        endif
      enddo

      deallocate(iused)
!--- print edges
!     knot vector, vert1, vert2

      write (12,'(/,a,/,i0)') 'edges', nedge
      do i=1, nedge
        write (12,240) gedge(1,i), gedge(2,i), gedge(3,i)
      enddo

!--- vertices

      write (12,'(/,a,/,i0)') 'vertices', nvert

!--- knotvectors

   ! we only use one knot vector, but listing 4 here

      nknot=4
      write (12,'(/,a,/,i0)') 'knotvectors', nknot
      do i=1, nknot
         write (12,'(a)') '2 3 0 0 0 1 1 1'
      enddo

!--- print weights

      !wcon1=0.853553d0   ! **** 3D element weight   todo:

      write (12,'(/,a)') 'weights'
      do i=1, nvert
        if (weight(i).eq.1.0d0) then   ! vertice should always be 1?
          write (12,'(a)') '1'
        else
          write (12,'(f9.6)') weight(i)
        endif
      enddo
      do i=1, nedge
        if (wmid(i).eq.1.0d0) then
          write (12,'(a)') '1'
        else
          write (12,'(f9.6)') wmid(i)
        endif
      enddo
      do i=1, nelem   ! for xielem array
        write (12,'(f9.6)') wcenter(i)
      enddo

!--- data

      write (12,*)
      write (12,'(a)') 'FiniteElementSpace'
      write (12,'(a)') 'FiniteElementCollection: NURBS2'
      write (12,'(a)') 'VDim: 2'
      write (12,'(a)') 'Ordering: 1'
      write (12,*)
!x    write (12,*) '# vertices:'   ! don't include comments in data
      do i=1, nvert
         write (12,250) xi(1,i), xi(2,i)
      enddo
!x    write (12,*) '# midpoints:'   ! don't include comments in data
      do i=1, nedge
         write (12,250) ximid(1,i), ximid(2,i)
      enddo
!x    write (12,*) '# element faces'   ! don't include comments in data
      do i=1, nelem   ! saved as separate array
         write (12,250) xielem(1,i), xielem(2,i)
      enddo
 250  format (2f10.6)

!--- end of file

      close (12)

      write (*,*) 'nbound=', nbound
      write (*,*) 'nedge= ', nedge
      write (*,*) 'nelem= ', nelem
      write (*,*) 'nvert= ', nvert
      write (*,*) 'nknot= ', nknot

      return
      end subroutine writemesh

!=======================================================================
!--- check if every edge is defined (debug)
!=======================================================================
      subroutine checkedge()
      use mod_nurbsmesh
      implicit none

      integer :: i, j
      integer :: j1, j2, j3, j4
      integer, allocatable :: iused(:)

      logical :: iflost

!--- start

      allocate (iused(nedge))
      iused=0   ! number of times an edge is used

      do i=1, nelem

          iflost=.true.
          j1=gg(1,i)
          j2=gg(2,i)
          do j=1, nedge
            j3=gedge(2,j)
            j4=gedge(3,j)
            if (j1.eq.j3 .and. j2.eq.j4) then
              iflost=.false.
              iused(j)=iused(j)+1
              local2global(i,1) = j
            endif
            if (j1.eq.j4 .and. j2.eq.j3) then
              iflost=.false.
              iused(j)=iused(j)+1
              local2global(i,1) = j
            endif
          enddo
          if (iflost) then
             write (0,*) '**** edge not defined ', j1, j2
          endif

          iflost=.true.
          j1=gg(2,i)
          j2=gg(3,i)
          do j=1, nedge
            j3=gedge(2,j)
            j4=gedge(3,j)
            if (j1.eq.j3 .and. j2.eq.j4) then
              iflost=.false.
              iused(j)=iused(j)+1
              local2global(i,2) = j
            endif
            if (j1.eq.j4 .and. j2.eq.j3) then
              iflost=.false.
              iused(j)=iused(j)+1
              local2global(i,2) = j
            endif
          enddo
          if (iflost) then
             write (0,*) '**** edge not defined ', j1, j2
          endif

          iflost=.true.
          j1=gg(3,i)
          j2=gg(4,i)
          do j=1, nedge
            j3=gedge(2,j)
            j4=gedge(3,j)
            if (j1.eq.j3 .and. j2.eq.j4) then
              iflost=.false.
              iused(j)=iused(j)+1
              local2global(i,3) = j
            endif
            if (j1.eq.j4 .and. j2.eq.j3) then
              iflost=.false.
              iused(j)=iused(j)+1
              local2global(i,3) = j
            endif
          enddo
          if (iflost) then
             write (0,*) '**** edge not defined ', j1, j2
          endif

          iflost=.true.
          j1=gg(4,i)
          j2=gg(1,i)
          do j=1, nedge
            j3=gedge(2,j)
            j4=gedge(3,j)
            if (j1.eq.j3 .and. j2.eq.j4) then
              iflost=.false.
              iused(j)=iused(j)+1
              local2global(i,4) = j
            endif
            if (j1.eq.j4 .and. j2.eq.j3) then
              iflost=.false.
              iused(j)=iused(j)+1
              local2global(i,4) = j
            endif
          enddo
          if (iflost) then
             write (0,*) '**** edge not defined ', j1, j2
          endif

        enddo

!d      write (*,*)
!d      write (*,*) 'number of times each edge is used'
!d      do j=1, nedge
!d        write (*,*) j, iused(j)
!d      enddo
        deallocate (iused)

      return
      end subroutine checkedge
