      subroutine pincell_octant(iouter, itype)
      use mod_nurbsmesh
      use mod_input, only : pintype, matcool
      implicit none
!=======================================================================
!
!  Create mesh for the interior of a single pincell 
!
!  Scott Palmtag
!  August 2026
!
! @version CVS $Id: pincell8.f90,v 1.7 2026/08/22 19:39:55 palmtag Exp $
!
!=======================================================================

!  create mesh for a single pincell with octant divisions and 2 rings
!    input: vertices for outer square that have already been defined

      integer, intent(in) :: iouter(8)   ! vertices on outer edge of pincell
      integer, intent(in) :: itype    ! pin type number

!  pass in the edge vertices: in array iouter
!   7 - 6 - 5
!   8       4
!   1 - 2 - 3

! outer vertices and edges have already been filled

! @version CVS $Id: pincell8.f90,v 1.7 2026/08/22 19:39:55 palmtag Exp $

      integer :: i
      integer :: i0, i2
      integer :: nr        ! ring counter
      integer :: ip
      integer :: istart
      integer :: ivcent    ! center vertice in pincell
      integer :: nesave
      integer :: nring

      real(8) :: radx      ! temp radius
      real(8) :: xc, yc    ! center of cell
      real(8) :: w0
      real(8) :: x0, x1, x2
      real(8) :: y0, y1, y2

      real(8), allocatable :: edgerad(:)   ! flag for rounded edges - edge radius

      real(8), parameter :: sq2=1.0d0/sqrt(2.0d0)

!d    write (*,*) '1/sqrt(2)=', sq2

!--- use new input

      write (*,*)
      write (*,*) 'subroutine pincell_octant'
      write (*,*) 'creating pincell map for type ', itype
      write (*,'(2a)') '   material ', trim(pintype(itype)%pname)

      nring=pintype(itype)%nring

! checks to make sure outer ring is in expected order
!d    write (*,*) 'debug: outer ring'
!d    write (*,*) iouter(7), iouter(6), iouter(5)
!d    write (*,*) iouter(8), 0,         iouter(4)
!d    write (*,*) iouter(1:3)

      if (iouter(1)+1.ne.iouter(2)) stop 'outer ring error pincell8 1'
      if (iouter(2)+1.ne.iouter(3)) stop 'outer ring error pincell8 2'
      if (iouter(8)+2.ne.iouter(4)) stop 'outer ring error pincell8 8'
      if (iouter(7)+1.ne.iouter(6)) stop 'outer ring error pincell8 7'
      if (iouter(6)+1.ne.iouter(5)) stop 'outer ring error pincell8 6'

      ivcent=iouter(8)+1   ! center vertice
      xc=xi(1,ivcent+1)    ! account for 0-based numbering
      yc=xi(2,ivcent+1)    ! account for 0-based numbering
!d    write (*,*) 'debug: center vertice ', ivcent
!d    write (*,*) 'debug: xc= ', xc
!d    write (*,*) 'debug: yc= ', yc

!--- add vertices around rings

  ! (8 vertices per ring)

      istart=nvert-1  ! save index and subtract 1 for 0 based

      do nr=nring, 1, -1
        radx=pintype(itype)%pinrad(nr)

        nvert=nvert+1
        xi(1,nvert)=xc-radx*sq2
        xi(2,nvert)=yc-radx*sq2
        nvert=nvert+1
        xi(1,nvert)=xc
        xi(2,nvert)=yc-radx
        nvert=nvert+1
        xi(1,nvert)=xc+radx*sq2
        xi(2,nvert)=yc-radx*sq2
        nvert=nvert+1
        xi(1,nvert)=xc+radx
        xi(2,nvert)=yc
        nvert=nvert+1
        xi(1,nvert)=xc+radx*sq2
        xi(2,nvert)=yc+radx*sq2
        nvert=nvert+1
        xi(1,nvert)=xc
        xi(2,nvert)=yc+radx
        nvert=nvert+1
        xi(1,nvert)=xc-radx*sq2
        xi(2,nvert)=yc+radx*sq2
        nvert=nvert+1
        xi(1,nvert)=xc-radx
        xi(2,nvert)=yc
      enddo

!--- define elements
!    elements are defined counter-clockwise (CCW)

      ip=istart      ! pointer for starting vertice
!d    write (*,*) 'debug: ip = ', ip
      do i=1, 8    ! outer ring of coolant
        nelem=nelem+1
        gg(1,nelem)=iouter(i)
        if (i.lt.8) then
          gg(2,nelem)=iouter(i+1)
          gg(3,nelem)=ip+i+1
        else
          gg(2,nelem)=iouter(1)  ! close the circle
          gg(3,nelem)=ip+1       ! close the circle
        endif
        gg(4,nelem)=ip+i
        matl(nelem)=matcool      ! coolant
      enddo

   ! 8 elements per fuel ring (except inside ring)

      do nr=nring, 2, -1   ! loop over rings
        do i=1, 8    ! ring of fuel
          nelem=nelem+1
          gg(1,nelem)=ip+i
          if (i.lt.8) then
            gg(2,nelem)=ip+i+1
            gg(3,nelem)=ip+i+9
          else
            gg(2,nelem)=ip+1   ! close the circle
            gg(3,nelem)=ip+9   ! close the circle
          endif
          gg(4,nelem)=ip+i+8
          matl(nelem)=pintype(itype)%pinmat(nr)       ! fuel
        enddo
        ip=ip+8
      enddo ! nr

  ! inside ring is only 4 elements

      nr=1   ! ring
      nelem=nelem+1    ! inside square of fuel
        gg(1,nelem)=ip+1
        gg(2,nelem)=ip+2
        gg(3,nelem)=ivcent
        gg(4,nelem)=ip+8
        matl(nelem)=pintype(itype)%pinmat(nr)       ! fuel
      nelem=nelem+1
        gg(1,nelem)=ip+2
        gg(2,nelem)=ip+3
        gg(3,nelem)=ip+4
        gg(4,nelem)=ivcent
        matl(nelem)=pintype(itype)%pinmat(nr)       ! fuel
      nelem=nelem+1
        gg(1,nelem)=ip+4
        gg(2,nelem)=ip+5
        gg(3,nelem)=ip+6
        gg(4,nelem)=ivcent
        matl(nelem)=pintype(itype)%pinmat(nr)       ! fuel
      nelem=nelem+1
        gg(1,nelem)=ip+6
        gg(2,nelem)=ip+7
        gg(3,nelem)=ip+8
        gg(4,nelem)=ivcent
        matl(nelem)=pintype(itype)%pinmat(nr)       ! fuel

!--- define edges

      allocate (edgerad(maxedge))
      edgerad(:)=-1000.0d0  ! negative is flat

      nesave=nedge   ! save for calculating round midpoints below

      ip=istart   ! vertice pointer
      do nr=nring, 1, -1 
        if (nr.lt.nring) ip=ip+8   ! update for each ring
        radx=pintype(itype)%pinrad(nr)

        nedge=nedge+1
        gedge(1,nedge)= 1   ! 1=horizontal
        gedge(2,nedge)= ip+1
        gedge(3,nedge)= ip+2
        edgerad(nedge)=radx
        nedge=nedge+1
        gedge(1,nedge)= 1
        gedge(2,nedge)= ip+2
        gedge(3,nedge)= ip+3
        edgerad(nedge)=radx
        nedge=nedge+1
        gedge(1,nedge)= 0   !  0=vertical
        gedge(2,nedge)= ip+3
        gedge(3,nedge)= ip+4
        edgerad(nedge)=radx
        nedge=nedge+1
        gedge(1,nedge)= 0
        gedge(2,nedge)= ip+4
        gedge(3,nedge)= ip+5
        edgerad(nedge)=radx
        nedge=nedge+1
        gedge(1,nedge)= 1
        gedge(2,nedge)= ip+6
        gedge(3,nedge)= ip+5
        edgerad(nedge)=radx
        nedge=nedge+1
        gedge(1,nedge)= 1
        gedge(2,nedge)= ip+7
        gedge(3,nedge)= ip+6
        edgerad(nedge)=radx
        nedge=nedge+1
        gedge(1,nedge)= 0
        gedge(2,nedge)= ip+8
        gedge(3,nedge)= ip+7
        edgerad(nedge)=radx
        nedge=nedge+1
        gedge(1,nedge)= 0
        gedge(2,nedge)= ip+1
        gedge(3,nedge)= ip+8
        edgerad(nedge)=radx

      enddo

      ! center cross  (4)
        nedge=nedge+1
        gedge(1,nedge)= 1   ! 1=horizontal
        gedge(2,nedge)= ip+8
        gedge(3,nedge)= ivcent
        nedge=nedge+1
        gedge(1,nedge)= 1
        gedge(2,nedge)= ivcent
        gedge(3,nedge)= ip+4
        nedge=nedge+1
        gedge(1,nedge)= 0   ! 0=vertical
        gedge(2,nedge)= ip+2
        gedge(3,nedge)= ivcent
        nedge=nedge+1
        gedge(1,nedge)= 0
        gedge(2,nedge)= ivcent
        gedge(3,nedge)= ip+6

      ! diagonals going out
        do i=1, 8            ! outer ring of coolang
          nedge=nedge+1
          gedge(1,nedge)= 3   ! diagonals going out
          gedge(2,nedge)= istart+i
          gedge(3,nedge)= iouter(i)
        enddo
        do nr=nring, 2, -1     ! don't include inner ring with 4 elements
          do i=1, 8
            nedge=nedge+1
            gedge(1,nedge)= 3   ! diagonals going out
            gedge(2,nedge)= istart+i+8
            gedge(3,nedge)= istart+i
          enddo
          istart=istart+8
       enddo

!--- calculate vertices and weights
!  weights have already been assigned as 1

!   calculate midpoint and midpoint weight
!   for straight lines, the midpoint is the average and the weight is 1.0
!   for circles, it is a little more complicated

! Equations for mid control point and weight from:
!  https://www.geometrictools.com/Documentation/NURBSCircleSphere.pdf
! Note that you I am using the inverse of the weight that he used.
! Also, make sure to convert the coordinates to a unity circle at the origin,
!  then apply the equation, then convert the circle back to the real circle.

      do i=nesave+1, nedge
!d      write (*,*) 'debug: edge ', i
        i0=gedge(2,i)+1   ! edges are zero based numbering
        i2=gedge(3,i)+1   ! edges are zero based numbering

        if (edgerad(i).lt.0.0d0) then
           ximid(1,i)=0.5*(xi(1,i0)+xi(1,i2))  ! midpoint
           ximid(2,i)=0.5*(xi(2,i0)+xi(2,i2))
        else

!d         write (*,*) 'qtr circle edge ', i
           x0=(xi(1,i0)-xc)/edgerad(i)    ! move so origin is at center
           x2=(xi(1,i2)-xc)/edgerad(i)
           y0=(xi(2,i0)-yc)/edgerad(i)
           y2=(xi(2,i2)-yc)/edgerad(i)
           x1=(y2-y0)/(x0*y2-x2*y0)
           y1=(x0-x2)/(x0*y2-x2*y0)
        ! this is actually the inverse weight that I use
           w0=sqrt(2.0d0*(x1*x1+y1*y1-1)/(1.0d0-x0*x2-y0*y2))
           ximid(1,i)=x1*edgerad(i)+xc   ! unnormalize
           ximid(2,i)=y1*edgerad(i)+yc   ! unnormalize
           wmid(i)=1.0d0/w0

        endif
      enddo

      deallocate (edgerad)

      return
      end subroutine pincell_octant

!=======================================================================


