module globals
  Implicit none
  !
  integer                         :: Seed = 0
  integer,parameter               :: maxSites = 5         ! maximun number of sites for a component
  integer,parameter               :: maxNpar = 10         ! maximun number of parameters for a model
  integer,parameter               :: maxdata = 1024       ! Maximun number of data in Compound file (pressure, loading, etc.)
  integer,parameter               :: maxcomponents = 10   ! maximun number of compounds in a mixture
  integer                         :: Ncomponents = 0
  !
  real                            :: MaxPressure = 1.0e+20
  real                            :: MinPressure = 1.0e-20
  real                            :: Chanels = 4096
  !
  ! Component structure
  type                            :: Components
    character(len=100)              :: IsothermFileName
    character(len=20)               :: Label 
    integer                         :: id
    integer                         :: sites              ! Number Of Sites
    character(len=20)               :: Models(maxSites)   ! Model of each site
    integer                         :: Npar(maxSites)     ! Number Of paramters per site
    integer                         :: TotalNpar
    real                            :: p(1:maxSites,1:maxNpar) ! Parameters
    integer                         :: NPressures              ! N Pressure Points in input isotherms
    real                            :: PurePressure(1:maxdata), PureLoading(1:maxdata)
    integer                         :: col1,col2
  end type
  type(Components), allocatable   :: Component(:)
  character(32*maxNpar),allocatable :: string_IEEE(:)
  !
  ! From Inputs:
  character(len=20)               :: SimulationType, DisplayName  ! Some system identification
  ! Flags for special actions:
  logical                         :: Fitting_Flag = .false.
  logical                         :: physical_constrains_flag = .false.
  logical                         :: seed_flag = .false.
  logical                         :: PressureRange_flag = .false.
  logical                         :: Refitting_Flag = .false.
  ! Timmings
  !real                            :: InitialisingTime
  !real                            :: FittingTime
  ! Isotherms
  !-----------
  Contains
! --------
  subroutine clean_screen                          ! Clean screen
   write(*,'(a1,a)',advance="no") achar(27), '[2J'
  end subroutine clean_screen
!
  pure character(len=32) function real2bin(a)      ! Convert real value in binary
   real,intent(in)                :: a
   write(real2bin(1:32) ,'(b32.32)') a
  end function real2bin
!
  pure real function bin2real(a)                   ! Convert binary to real
   character(len=32), intent(in)  :: a
   read(a(1:32),'(b32.32)') bin2real
  end function bin2real
!
  PURE INTEGER FUNCTION Clen(s)                    ! returns same result as LEN unless:
   CHARACTER(*),INTENT(IN)        :: s             ! last non-blank char is null
   INTEGER :: i
   Clen = LEN(s)
   i = LEN_TRIM(s)
   IF (s(i:i) == CHAR(0)) Clen = i-1
  END FUNCTION Clen
!
  PURE INTEGER FUNCTION Clen_trim(s)               ! returns same result as LEN_TRIM unless:
   CHARACTER(*),INTENT(IN)        :: s             ! last char non-blank is null, if true:
   INTEGER :: i                                    ! then len of C string is returned, note: "Ctrim is only user of this function"
   i = LEN_TRIM(s) ; Clen_trim = i
   IF (s(i:i) == CHAR(0)) Clen_trim = Clen(s)
  END FUNCTION Clen_trim
!
  subroutine ReadInputFile()
    ! Read input file
    ! Template for input.data
    !------------------------------------------------------------------------------------------------
    !    SimulationType           Fitting                           ! 
    !
    !    DisplayName              MAF-X8
    !    ModelConstrained
    !
    !    RestartFit               .true.
    !    Component 0 p-xylene
    !                CompoundPureIsothermFileName         Results.dat-MAF-x8-P1-Repeat-433K-p-xylene
    !                CompoundPureIsothermColumnPressure   3
    !                CompoundPureIsothermColumnLoading    8
    !                NumberOfIsothermSites      2
    !                 Langmuir-Freundlich        0.0  0.0     0.0
    !                 Langmuir-Freundlich        0.0  0.0     0.0
    !------------------------------------------------------------------------------------------------
    implicit none
    ! open files
    integer                       :: err_opening
    integer                       :: i,j,k,s
    character(len=100)            :: line
    logical                       :: FileNameCompoundFlag      = .false.
    logical                       :: NumberOfSitesCompoundFlag = .false.
    character(len=9)              :: spam
    read_input_do: do
      read(5,'(A)',iostat=err_opening) line
      ! End reading or ignore lines
      ! General stuff
      if(err_opening/=0)                      exit read_input_do
      if(line(1:1)=='#'.or.line(1:1)=='!')   cycle read_input_do
      if(line(1:3)=="end")                    exit read_input_do
      !
      if(line(1:14)=="SimulationType") then
        read(line(16:Clen_trim(line)),'(a)') SimulationType
        SimulationType = adjustl(SimulationType)
        if (SimulationType=='Fitting') Fitting_Flag = .true.
      end if
      ! Pressure range:
      if(line(1:13)=="PressureRange") then
        read(line(14:),*) MinPressure, MaxPressure
        PressureRange_flag = .true.
      end if
      ! Refitting flag
      if(line(1:10)=="RestartFit") then
        Refitting_Flag = .true.
      end if
      ! System name
      if(line(1:17)=="DisplayName") then
        read(line(16:Clen_trim(line)),'(a)') DisplayName
        DisplayName = adjustl(DisplayName)
      end if
      ! Model parameters have some constrains 
      if(line(1:16)=='ModelConstrained') physical_constrains_flag = .true.
      if(line(1:4)=='Seed') seed_flag = .true.
      ! Component section:
      ! Initialisation:
      if(line(1:9)=="Component") then
        Ncomponents = Ncomponents + 1
      end if
    end do read_input_do
    rewind(5)
    !
    ! Allocation Section:
    allocate(Component(1:Ncomponents))
    allocate(string_IEEE(1:Ncomponents))
    !
    i = 0
    ! Read Component Section:
    read_input_do_again: do 
      read(5,'(A)',iostat=err_opening) line
      if(err_opening/=0)                      exit read_input_do_again
      if(line(1:9)=="Component") then
        i = i + 1
        read(line,'(a,i2,a)') spam, Component(i)%id, Component(i)%Label
        write(6,'(a,i3,1x,a,1x,a)') 'Component',Component(i)%id,'is:',Component(i)%Label
        !
        FileNameCompoundFlag = .false.
        NumberOfSitesCompoundFlag = .false.
        !
        read_component: do
          read(5,'(A)',iostat=err_opening) line; line = trim(adjustl(line))
          if (line(1:28) == "CompoundPureIsothermFileName") then
            read(line(29:),'(a)') Component(i)%IsothermFileName 
            Component(i)%IsothermFileName = trim(adjustl(Component(i)%IsothermFileName))
            FileNameCompoundFlag = .true.
            write(6,'(a,1x,a)')' Filename is:', Component(i)%IsothermFileName
          end if
          !
          if (line(1:34)=="CompoundPureIsothermColumnPressure") then
            read(line(35:),*) Component(i)%col1
            write(6,'(a,1x,i2)')' Column Pressure: ', Component(i)%col1
          end if 
          !
          if (line(1:33)=="CompoundPureIsothermColumnLoading") then
            read(line(34:),*) Component(i)%col2
            write(6,'(a,1x,i2)')' Column Loading: ', Component(i)%col2
          end if
          !
          if (line(1:21) == "NumberOfIsothermSites") then
            read(line(22:),*) Component(i)%sites
            write(6,'(a,i3)') ' Number of sites:', Component(i)%sites
            s = 0
            do j = 1, Component(i)%sites
              read(5,'(A)',iostat=err_opening) line; line = trim(adjustl(line))
              read(line,*) Component(i)%Models(j), spam
              Component(i)%Models(j) = trim(adjustl(Component(i)%Models(j)))
              call CheckModelSite( Component(i)%Models(j), Component(i)%Npar(j))
              write(6,'(a,i1,a,i3)') ' Number of parameters for site(',j,') :',Component(i)%Npar(j)
              ! total number of parameters:
              Component(i)%TotalNpar = sum(Component(i)%Npar(1:Component(i)%sites))
              ! If Refitting_Flag is activated:
              !  Read parameter from input file:
              if( Refitting_Flag ) then
               read(line,*) Component(i)%Models(j), ( Component(i)%p(j,k), k =1, Component(i)%Npar(j) )
               write(6,'(a,1x,100(f25.12,1x))') Component(i)%Models(j), ( Component(i)%p(j,k), k =1, Component(i)%Npar(j) )
               ! dev: check binary numbers:
               do k = 1, Component(i)%Npar(j)
                s = s + 1
                string_IEEE(i)(32*(s-1)+1:32*s) = real2bin( Component(i)%p(j,k) )
               end do
              end if
            end do
            if( Refitting_Flag ) write(6,'(a,1x,i3)') string_IEEE(i)(1:32*Component(i)%TotalNpar), Component(i)%TotalNpar
            ! Dev: Check error in the convesion to real and binary and viceversa
            !s = 0
            !do j=1,Component(i)%sites
            !  do k = 1, Component(i)%Npar(j)
            !    s = s + 1
            !    write(6,'(f14.7)') 0.5*abs(Component(i)%p(j,k)-bin2real(string_IEEE(i)(32*(s-1)+1:32*s)))**2
            !  end do 
            !end do 
            !
            NumberOfSitesCompoundFlag = .true.
          end if
          if (FileNameCompoundFlag.and.NumberOfSitesCompoundFlag) exit read_component
        end do read_component
      end if
    end do read_input_do_again
    !
    ! Output the input file
    write(6,'(a)')'-----------------------------------------------------------'
    write(6,'(a,1x,a)') 'Simulation type:', SimulationType
    if ( physical_constrains_flag ) write(6,'(a)') 'The fitting is constrained'
    write(6,'(a,1x,a)') 'Display Name:', DisplayName
    write(6,'(a,1x,i3,a,i3,a)')'Number of components:', Ncomponents,' with ',&
      sum( Component(1:Ncomponents)%sites ),' total sites'
    !
    write(6,'(a)') 'Input file readed.'
    write(6,'(a)')'-----------------------------------------------------------'
    ! 
    ! Reading the isotherm files:
    write(6,'(a)') 'Reading isotherms...'
    ! Read reference data files for fitting:
    if(SimulationType=="Fitting") then
      do i = 1, Ncomponents
        call ReadIsotherms(i)
      end do 
    end if
    return
  end subroutine ReadInputFile
! 
  subroutine CheckModelSite(Model,Npar) 
    ! Read the number of parameters for allocation variables for each site and model
    ! Maybe move this to GA module?
    implicit none 
    character(len=20),intent(in) :: Model 
    integer,intent(out)          :: Npar 
    select case(Model)
    case("Langmuir","L")
      Npar = 2
    case("Langmuir-Freundlich","LF")
      Npar = 3
    case("Sips")
      Npar = 3
    case("Toth")
      Npar = 3
    case("Unilan")
      Npar = 3
    case("Asymptotic-Temkin")
      Npar = 3
    case default
      STOP 'Model Not Found'
    end select
  end subroutine CheckModelSite
!
  subroutine ReadIsotherms(ComponentTmp)
    ! Read all isotherm files
    implicit none
    integer,intent(in)            :: ComponentTmp
    character(len=500)            :: line
    integer                       :: u                   ! Unit File
    integer                       :: err_opening
    integer                       :: i
    real                          :: junk(4)
    integer                       :: col1,col2
    !
    col1 = Component(ComponentTmp)%col1
    col2 = Component(ComponentTmp)%col2
    !
    Component(ComponentTmp)%NPressures = 0
    write(6,'(a,a)') '[...] Reading isotherm file: ', Component(ComponentTmp)%IsothermFileName
    open( NewUnit=u, file=Component(ComponentTmp)%IsothermFileName, iostat=err_opening)
    if ( err_opening/=0 ) stop 'Error openning Isotherms'
    !
    read_isotherm: do 
      read(u,*,iostat=err_opening) line ; line=trim(adjustl(line))
      if( err_opening /= 0) exit read_isotherm
      if( line(1:1) == '#' ) cycle read_isotherm 
      Component(ComponentTmp)%NPressures = Component(ComponentTmp)%NPressures + 1
    end do read_isotherm
    rewind(u)
    !
    write(6,*)' Number of pressure lines:', Component(ComponentTmp)%NPressures
    write(6,*)' Columms:', col1, col2
    ! Read pure pressures and loadings:
    i = 0
    read_isotherm_again: do !i = 1, Component(ComponentTmp)%NPressures
      read(u,'(a)',iostat=err_opening) line
      if ( err_opening/=0 )  exit  read_isotherm_again
      if( line(1:1) == '#' ) cycle read_isotherm_again
      i = i + 1
      if( i > Component(ComponentTmp)%NPressures ) exit  read_isotherm_again
      read(line,*) junk(:col1-1), Component(ComponentTmp)%PurePressure(i), &
                junk(:col2-col1-1), Component(ComponentTmp)%PureLoading(i)
      write(6,'(2(f24.12,1x))') Component(ComponentTmp)%PurePressure(i), &
                                Component(ComponentTmp)%PureLoading(i)
    end do read_isotherm_again
    close(u)
    return
  end subroutine ReadIsotherms
!
  subroutine Bounds()
    Implicit None
    integer                 :: i, j
    real                    :: dx
    real                    :: LogMaxPressure
    real                    :: LogMinPressure
    !
    if(PressureRange_flag)then
      write(6,'(a)')'Pressure range from input file'
      write(6,*)'Min. and Max. pressures:', MinPressure, MaxPressure
    else
      do i = 1,Ncomponents
        do j = 1,Component(i)%NPressures
          if( Component(i)%PurePressure(j) >= MaxPressure ) MaxPressure = Component(i)%PurePressure(j)
          if( Component(i)%PurePressure(j) <= MinPressure ) MinPressure = Component(i)%PurePressure(j)
        end do
      end do
      write(6,'(a)')'Pressure range from extreme values in isotherms.'
      write(6,*)'Min. and Max. pressures:', MinPressure, MaxPressure
    end if
    !
    LogMaxPressure = log(MaxPressure)
    LogMinPressure = log(MinPressure)
    dx = real((LogMaxPressure-LogMinPressure)/real(Chanels))   ! log-space
    return
  end subroutine bounds
  !
  subroutine deallocateAll()
    deallocate(Component)
    deallocate(string_IEEE)
  end subroutine deallocateAll
end module globals
