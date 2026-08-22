   module mod_input
      implicit none
!-----------------------------------------------------------------------
!
!  Module to store input values for NURBS mesh generator
!
!  Scott Palmtag
!  August 2026
!
! @version CVS $Id: Mod_input.f90,v 1.4 2026/08/22 03:01:37 palmtag Exp $
!
!-----------------------------------------------------------------------

      integer, protected :: nrow=0      ! number of rows of pins on one edge
      integer, protected :: matcool         ! coolant material number

      real(8), protected :: ppitch          ! ppitch
      real(8), protected :: apitch          ! total assembly pitch (flat to flat)

      character(len=200) :: fbase    ! base input file name without suffix
      character(len=200), protected :: title    ! title
      character(len=200), protected :: title2   ! title

!  internal

      real(8), parameter :: smallgap=0.1d0  ! gap for inserts

      integer, protected, allocatable :: irodmap(:,:)

      integer, parameter :: maxpintype=10
      integer, parameter :: maxring=10
      integer, parameter :: maxmatl=10

      character(len=20) :: matname(maxmatl)

      integer :: npintype

      type :: pintype_type
         integer :: nring           ! number of rings
         integer :: pinmat(maxring) ! ring material
         real(8) :: pinrad(maxring) ! ring radii
         character(len=4) :: pname  ! pin type name
      end type pintype_type
      type(pintype_type), protected :: pintype(maxpintype)

   contains
!=======================================================================
!
!  Subroutine to read input file for mesh generator
!
!=======================================================================
      subroutine readinput(fname)
      implicit none

      character(len=*), intent(in) :: fname

!--- locals

      integer :: linp=22    ! input unit number
      integer :: ibad       ! error flag counter
      integer :: i, j, k
      logical :: ifrodmap
      logical :: ifxst

      integer :: itmp(maxring)
      real(8) :: xtmp(maxring)

      real(8) :: xlen

      integer, parameter :: maxrod=600
      character(len=4), allocatable :: crodmap(:,:)  ! character rod map

      character(len=4)   :: ctmp
      character(len=12)  :: card
      character(len=100) :: line

!--- initialize

!d    write (0,*) 'trace: start input'

      ppitch=-100.0d0
      apitch=-100.0d0

      matcool=0

      title=' '
      title2=' '

      matname=' '   ! list of materials

      npintype=0
      do i=1, maxpintype
        pintype(i)%nring=0         ! number of rings
        pintype(i)%pinmat(:)=0     ! ring material
        pintype(i)%pinrad(:)=0.0d0 ! ring radii
        pintype(i)%pname=' '       ! pin type name
      enddo

      ibad=0

!--- read input

      write (*,'(2a)') ' Reading input file: ', trim(fname)
      inquire (file=fname, exist=ifxst)
      if (.not.ifxst) stop 'input file does not exist'

      open (linp, file=fname, status='old', action='read')

      do
        read (linp,'(a)',end=800) line
        k=index(line,'!')
        if (k.gt.0) line(k:)=' '    ! remove comment
        k=index(line,'#')
        if (k.gt.0) line(k:)=' '    ! remove comment
        if (line.eq.' ') cycle

        read (line,*) card

!> name: title
!> description: title of case
        if     (card.eq.'title') then
          read (line,*) card, title
        elseif (card.eq.'title2') then
          read (line,*) card, title2

!> name: square
!> description: obsolete card
        elseif (card.eq.'square') then
          read (line,*) card

!> name: nrow
!> description: number of rods along one side
        elseif (card.eq.'nrow') then
          if (nrow.gt.0) stop 'only one nrow card allowed'
          read (line,*) card, nrow
          allocate (irodmap(nrow,nrow))
          allocate (crodmap(nrow,nrow))
          irodmap=0
          crodmap=' '

!> name: ppitch
!> description: pin pitch
        elseif (card.eq.'ppitch') then
          read (line,*) card, ppitch

!> name: apitch
!> description: assembly pitch
        elseif (card.eq.'apitch') then
          read (line,*) card, apitch

!> name: matcool
!> description: coolant material number
        elseif (card.eq.'matcool') then
          read (line,*) card, matcool

!> name: pinrad
!> description: list of pin radii for a single pin type
        elseif (card.eq.'pinrad') then
          xtmp=-100.0d0
          read (line,*) card, ctmp, xtmp(:)
          call fill_pinrad(ctmp, xtmp)

!> name: pinmat
!> description: list of pin materials for a single pin type
        elseif (card.eq.'pinmat') then
          itmp=-100
          read (line,*) card, ctmp, itmp(:)
          call fill_pinmat(ctmp, itmp)

!> name: material
!> description: material name
!> notes: the material name is an edit printed to top of mesh file
        elseif (card.eq.'matname') then
          itmp=-100
          read (line,*) card, i, matname(i)

!> name: rodmap
!> description: array of rod types, start on next line
        elseif (card.eq.'rodmap') then
          ifrodmap=.true.
          read (line,*) card
          do j=1, nrow
            read (linp,*) (crodmap(i,j),i=1,nrow)
          enddo

        else
          write (0,   *) 'ERROR: invalid input card ', card
          stop 'invalid input card'
        endif

      enddo
  800 continue

      close (linp)

!--- write input

      write (*,'(2x,2a)') 'title:  ', trim(title)
      if (title2.ne.' ') then
        write (*,'(2x,2a)') 'title2: ', trim(title2)
      endif

      write (*,24) 'nrow      ', nrow
      write (*,24) 'matcool   ', matcool,' material coolant'
      do i=1, maxmatl
        if (matname(i).ne.' ') then
          write (*,24) 'material ', i, trim(matname(i))
        endif
      enddo
   20 format (2x,a,f12.5)
   24 format (2x,a,1x,i0,2x,a)

      if (apitch.lt.0.0d0) then
        apitch=nrow*ppitch
        write (*,*) 'setting default apitch for square problem'
      endif
      if (nrow.eq.1) then
        if (abs(apitch-ppitch).gt.0.000001d0) then
           write (*,*) 'invalid apitch for 1x1 problem'
           stop 'invalid apitch for 1x1 problem'
        endif
      endif

      write (*,20) 'ppitch ', ppitch
      write (*,20) 'apitch ', apitch

      if (ppitch.le.0.0d0) then
        write (0,*) 'no pin pitch specified'
        ibad=ibad+1
      endif

      xlen=nrow*ppitch
      if (apitch.lt.xlen) then
        write (0,*) 'assembly pitch is too small'
        write (0,*) 'apitch must be at least ', xlen
        ibad=ibad+1
      endif

      if (matcool.le.0) then
        write (0,*) 'invalid matcool'
        ibad=ibad+1
      endif

      if (ibad.gt.0) stop 'input errors encountered - check output'

!--- pin types

      if (npintype.le.0) then
        stop 'no pin descriptions have been entered in input'
      endif

      write (*,*)
      do i=1, npintype
        associate (pp => pintype(i))
          write (*,120) pp%pname, (pp%pinrad(j),j=1,pp%nring)
          write (*,125) pp%pname, (pp%pinmat(j),j=1,pp%nring)
        end associate
      enddo
  120 format ('  pinrad ', a, 50f9.5)
  125 format ('  pinmat ', a, 50(i6,3x))

!--- calculate number of rods in square geometry

      write (*,*)
      write (*,*) 'input rod map:'
      do j=1, nrow
         write (*,44) (crodmap(i,j),i=1,nrow)
         do i=1, nrow
           if (crodmap(i,j).eq.' ') stop 'blank found in pin map'
         enddo
      enddo
      write (*,*)
  44  format (4x,50a4)

!--- fill irodmap from character map

      if (ifrodmap) then
        do j=1, nrow
          do i=1, nrow
            irodmap(i,j)=find_pintype(crodmap(i,j),0)
          enddo
        enddo
      else
        irodmap=1    ! default to first pin type entered
      endif

!--- finished

      if (ibad.gt.0) stop 'input errors encountered - check output'

      deallocate(crodmap)

      return
      end subroutine readinput

!=======================================================================
!
!  add pinmat from input to data structure
!
      subroutine fill_pinmat(ctmp, itmp)
      implicit none
      character(len=4), intent(in) :: ctmp
      integer,          intent(in) :: itmp(maxring)

      integer :: itype
      integer :: ir
      integer :: i

!--- find pintype

      itype=find_pintype(ctmp, 1)

!---- find number of rings

      ir=0
      do i=maxring, 1, -1
        if (itmp(i).gt.0) then
          ir=i
          exit
        endif
      enddo
      if (ir.eq.0) stop 'no rings defined in input'
      if (pintype(itype)%nring.eq.0) then
        pintype(itype)%nring=ir
      else
        if (ir.ne.pintype(itype)%nring) stop 'unmatched rings in material'
      endif

!--- fill

      do i=1, ir
        pintype(itype)%pinmat(i)=itmp(i)
      enddo

      return
      end subroutine fill_pinmat

!=======================================================================
!
!  add pinrad from input to data structure
!
      subroutine fill_pinrad(ctmp, xtmp)
      implicit none
      character(len=4), intent(in) :: ctmp
      real(8),          intent(in) :: xtmp(maxring)

      integer :: itype
      integer :: ir
      integer :: i

!--- find pintype

      itype=find_pintype(ctmp, 1)

!---- find number of rings

      ir=0
      do i=maxring, 1, -1
        if (xtmp(i).gt.0.0d0) then
          ir=i
          exit
        endif
      enddo
      if (ir.eq.0) stop 'no rings defined in input'
      if (pintype(itype)%nring.eq.0) then
        pintype(itype)%nring=ir
      else
        if (ir.ne.pintype(itype)%nring) stop 'unmatched rings in material'
      endif

!--- fill

      do i=1, ir
        pintype(itype)%pinrad(i)=xtmp(i)
      enddo

      return
      end subroutine fill_pinrad

!=======================================================================
!
!  search list of pin types to convert character input to integer
!
      integer function find_pintype(ctmp, iopt)
      implicit none
      character(len=4), intent(in) :: ctmp
      integer,          intent(in) :: iopt   ! 1=allow new, 0 no

      integer :: itype
      integer :: i

!--- find pintype

      itype=0
      do i=1, npintype
        if (ctmp.eq.pintype(i)%pname) itype=i
      enddo
      if (itype.eq.0 .and. iopt.eq.1) then   ! create new type
        npintype=npintype+1
        itype=npintype
        if (npintype.gt.maxpintype) stop 'maxpintype exceeded'
        pintype(itype)%pname=ctmp
      endif
      if (itype.eq.0) then
        write (*,*) 'could not find pin type ', ctmp
        stop 'could not find pin type'
      endif

      find_pintype=itype

      return
      end function find_pintype

!=======================================================================
   end module mod_input
