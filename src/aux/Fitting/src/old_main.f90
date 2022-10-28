module general
    implicit none
    !
    integer                         :: Seed = 0
    integer,parameter               :: maxnp = 10
    integer                         :: Npar
    real,allocatable                :: param(:,:)
    real,allocatable                :: np(:)
    !
    integer,parameter               :: maxcompounds = 10
    integer                         :: compound
    integer                         :: Ncomponents
    !
    integer,parameter               :: maxdata = 1000      ! Maximun number of data in Compound file (pressure, loading, etc.)
    ! From Inputs:
    character(len=20)               :: SimulationType, DisplayName   ! Some system identification
    logical                         :: physical_constrains = .false. ! 
    ! Isotherms
    !-----------
    Contains
  ! --------
    subroutine clean_screen                      ! Clean screen
     write(*,'(a1,a)',advance="no") achar(27), '[2J'
    end subroutine clean_screen
  !
    pure character(len=32) function real2bin(a)  ! Convert real value in binary
     real,intent(in) :: a
     write(real2bin(1:32) ,'(b32.32)') a
    end function real2bin
  !
    pure real function bin2real(a)               ! Convert binary to real
     character(len=32), intent(in) :: a
     read(a(1:32),'(b32.32)') bin2real
    end function bin2real
  !
    PURE INTEGER FUNCTION Clen(s)                ! returns same result as LEN unless:
    CHARACTER(*),INTENT(IN) :: s                 ! last non-blank char is null
    INTEGER :: i
    Clen = LEN(s)
    i = LEN_TRIM(s)
    IF (s(i:i) == CHAR(0)) Clen = i-1
    END FUNCTION Clen
  !
    PURE INTEGER FUNCTION Clen_trim(s)           ! returns same result as LEN_TRIM unless:
    CHARACTER(*),INTENT(IN) :: s                 ! last char non-blank is null, if true:
    INTEGER :: i                                 ! then len of C string is returned, note:
                                                 !  "Ctrim is only user of this function"
    i = LEN_TRIM(s) ; Clen_trim = i
    IF (s(i:i) == CHAR(0)) Clen_trim = Clen(s)
    END FUNCTION Clen_trim
  !
  end module general
  
  module gaiast_globals
   use mod_random
   use general
   implicit none
  ! IAST
   integer,parameter               :: integration_points = 10000
   real,parameter                  :: precision_Newton = 1e-5  
   integer                         :: intervalos = 4096
   real                            :: tol = 0.001, tolfire = 0.25
   real,allocatable                :: concy(:),pi(:,:),iso(:,:)
   integer,allocatable             :: npress(:)
   real,target                     :: datas(2,maxcompounds,maxdata),f(maxdata)
   real,pointer                    :: x(:),y(:),alldat(:,:)
   character(100),allocatable      :: ajuste(:)
   character(100)                  :: String,intmethod,Equation
   character(32*maxnp),allocatable :: string_IEEE(:)
   character(5)                    :: inpt
   logical                         :: flag = .true., FlagFire = .false.,seed_flag=.true.
   logical                         :: physical_constrains = .false., range_flag =.true.
   logical                         :: Always_Use_Multicomponent_Solver=.false., IAST_flag = .true.
   logical                         :: refit_flag = .false.
   real,parameter                  :: R = 0.008314472 ! kJ/mol/K
   real                            :: T = 298.0000000 ! K
   real                            :: inferior
   real                            :: superior
   !
   Contains
  !--------
    subroutine ReadInputFile()
  ! Template for input.data
  !
  !    SimulationType           Fitting                           ! 
  !
  !    DisplayName              MAF-X8
  !    ModelConstrained
  !    Component 0 p-xylene
  !                CompoundPureIsothermFileName         Results.dat-MAF-x8-P1-Repeat-433K-p-xylene
  !                CompoundPureIsothermColumnPressure   3
  !                CompoundPureIsothermColumnLoading    8
  !                NumberOfIsothermSites      2                   ! NumberOfSites
  !                 Langmuir-Freundlich        0.0  0.0     0.0   ! SiteModel[1], SiteModel[1].Npar
  !                 Langmuir-Freundlich        0.0  0.0     0.0   ! SiteModel[2]
    implicit none
    ! open files
    integer                         :: err_opening
    character(len=100)              :: line
    !
    read_input_do: do
     read(5,'(A)',iostat=err_opening) line
     ! End reading or ignore lines
     ! General Stuff
     if(err_opening/=0)                                 exit read_input_do
     if(line(1:1)=='#'.or.line(1:1)=='!')              cycle read_input_do
     if(line(1:3)=="end".or.line(1:)=="end_iast_input") exit read_input_do
     !
     if(line(1:14)=="SimulationType") then
      read(line(16:Clen_trim(line)),'(a)') SimulationType
      SimulationType = adjustl(SimulationType)
     end if
     if(line(1:17)=="DisplayName") then
      read(line(16:Clen_trim(line)),'(a)') DisplayName
     end if
     ! Model parameters have some constrains 
     if(line(1:16)=='ModelConstrained') physical_constrains = .true.
  
     ! Component section:
     ! Initialisation:
     Ncomponents = 0
     !
     if(line(1:9)=="Component") then
      Ncomponents = Ncomponents + 1
     end if
    end do read_input_do
    ! Allocation Section:
    ! C
    allocate()
    ! Read Component Section:
  
     !--------------------------------------------------
     ! GAIAST version
     !if(line(1:5)=='ncomp')then
     ! read(line,*)inpt,ncomponents
     ! allocate(ajuste(ncomponents),concy(ncomponents))
     ! allocate(npress(ncomponents),np(ncomponents))
     ! allocate(string_IEEE(ncomponents))
     ! do ii=1,ncomponents
     !  read(5,*)concy(ii),ajuste(ii)
     !  string = ajuste(ii)
     !  call MakeInitPOP(string,npar)
     !  np(ii) = npar
     !  write(6,'(a,1x,i3,1x,a,1x,i3,1x,a)')'compound:',ii,'with',np(ii),'parameters'
     ! end do
     ! concy=concy/sum(concy)
     ! allocate(param(ncomponents,0:maxval(np)-1))
     ! param = 0.0
     ! cycle read_input_do
     !end if
     !if(line(1:5)=='ffit?')then
     ! read(line,*)inpt,flag,realorbinary
     ! if(flag.eqv..false.)then
     !  write(6,*)'[WARN] Fit is already done'
     !  select case(realorbinary)
     !   case('real')
     !    readparameter_real: do ii=1,ncomponents
     !     read(5,*)(param(ii,j),j=0,np(ii)-1)
     !    end do readparameter_real
     !   case default
     !    readparameter_bin: do ii=1,ncomponents
     !     do j=1,np(ii)
     !      read(5,'(a32)') string_IEEE(ii)(32*(j-1)+1:32*j)
     !      param( ii,j-1 ) = bin2real( string_IEEE(ii)(32*(j-1)+1:32*j) )
     !     end do
     !    end do readparameter_bin
     !   end select
     ! end if
      !cycle read_input_do
     !end if
     !if(line(1:5)=='refit') then
     ! write(6,*)'[WARN] Refitting parameters'
     ! read(line,*)inpt,refit_flag,realorbinary
     ! if(refit_flag)then
     !  select case(realorbinary)
     !  case('real')
     !   do ii=1,ncomponents
     !    read(5,*)(param(ii,j),j=0,np(ii)-1)
     !   end do
     !   do ii=1,ncomponents
     !    write(6,'(a,1x,i3,1x,a,1x,i3,1x,a)')'compound:',ii,'with',np(ii),'parameters'
     !    do j=1,np(ii)
     !     string_IEEE(ii)(32*(j-1)+1:32*j)=real2bin( param( ii,j-1 ) )
     !     write(6,'(i3,1x,i3,1x,f14.7,1x,a32)') ii,j,param( ii,j-1 ), string_IEEE(ii)(32*(j-1)+1:32*j)
     !    end do
     !   end do
     !  case default
     !   do ii=1,ncomponents
     !    write(6,'(a,1x,i3,1x,a,1x,i3,1x,a)')'compound:',ii,'with',np(ii),'parameters'
     !    do j=1,np(ii)
     !     read(5,'(a32)') string_IEEE(ii)(32*(j-1)+1:32*j) 
     !     param( ii,j-1 )=bin2real( string_IEEE(ii)(32*(j-1)+1:32*j) )
     !     write(6,'(i3,1x,i3,1x,f14.7,1x,a32)') ii,j,param( ii,j-1 ), string_IEEE(ii)(32*(j-1)+1:32*j)
     !    end do
     !   end do
     !  end select
     !  write(6,'(a)') string_IEEE(1:ncomponents)(1:32*np(ii))
     ! else
     !  string_IEEE(1:ncomponents)=' '
     ! end if
     !end if
     !if(line(1:5)=='inter') then
     ! read(line,*)inpt,intervalos
     !end if
     !if(line(1:5)=='toler') read(line,*)inpt, Tol
     !if(line(1:5)=='Range')then
     ! read(line,*)inpt, inferior, superior
     ! range_flag = .false.
     ! write(6,'(a)')'[WARN] Fugacity range specified'
     !end if
     !if(line(1:5)=='RSeed') then
     ! seed_flag=.false.
     ! read(line,*)inpt, seed
     !end if
    !end do read_input_do
    ! Allocate defaults options:
    !allocate(pi(ncomponents,0:intervalos),iso(ncomponents,0:intervalos))
   end subroutine ReadInputFile
  ! ...
   subroutine bounds()
    implicit none
    real               :: xxx(ncomponents),yyy(ncomponents),dx,x0,y0
    real               :: s1(ncomponents)
    real,allocatable   :: apar(:)
    integer            :: ib,jb
    character(len=100) :: funk
    s1 = 0.0
    if (range_flag) then
     superior = 0.0
     inferior = 1.0e30
     do ib = 1,ncomponents
      do jb = 1,npress(ib)
       if( datas(1,ib,jb) >= superior) superior = datas(1,ib,jb)
       if( datas(1,ib,jb) <= inferior) inferior = datas(1,ib,jb)
      end do
     end do
    end if
  !  ...
    write(6,*)'Inf:',inferior
    write(6,*)'Sup:',superior
  !  ...
    superior = log( superior )
    inferior = log( inferior )
    dx = real((superior-inferior)/real(intervalos))   ! log-space
  ! ...
    open(222,file='curves.txt')
    write(222,*)'#', inferior, superior
    do ib=0,intervalos-1
     do compound=1,ncomponents
      funk = ajuste( compound )
      allocate(apar(0:np(compound)-1))
      apar = 0.0
      do jb = 0, np(compound) - 1
       apar(jb) = param(compound,jb)
      end do
      xxx(compound) = exp( inferior + ib*dx )            ! real space
      x0=xxx(compound)
      yyy(compound) = model(apar, np(compound), x0, funk )
  ! integration of the spread pressure:
      y0 = exp( inferior + (ib+1)*dx )
      x0 = model(apar, np(compound), y0, funk )
      x0 = (yyy(compound) + x0)*dx/2.0
      s1(compound) = s1(compound) + x0
  ! spread pressure: definition
      pi(compound,ib+1) = s1(compound)
  ! adsorption value:
      iso(compound,ib+1) = yyy(compound)
      deallocate(apar)
     end do
     write(222,*)(xxx(compound) ,yyy(compound), pi(compound,ib+1), compound=1,ncomponents)
    end do
    close(222)
    return
   end subroutine bounds
  !
   pure real function integrate(x0,x1,a,n,funk)
    implicit none
    integer              ::  i
    integer,intent(in)   ::  n
    real,intent(in)      ::  x0,x1     ! intervalo (real)
    real,intent(in)      ::  a(0:n-1)  ! 
    real                 ::  delta,x
    real                 ::  factor
    character(100),intent(in):: funk
    delta=(x1-x0)/(integration_points)
    area: do i = 0,integration_points
     factor = 1.0
     if (i==0.or.i==integration_points-1) factor = 3.0/8.0
     if (i==1.or.i==integration_points-2) factor = 7.0/6.0
     if (i==2.or.i==integration_points-3) factor = 23.0/24.0
     x = x0 + delta*(0.5+i)
     Integrate = Integrate + factor*delta*model(a,n,x,funk)/x
    end do area
    return
   end function integrate
  
   subroutine ReadIsotherms()
    !character(len=2)                :: zlz
    !character(len=4)                :: extension=".gin"
    !GULPFilename=CIFFiles%filename(1:Clen_trim(CIFFiles%filename)-4)//extension
    !GULPFilename=adjustl(GULPfilename)
    !open(newunit=u,file=GULPFilename)
   implicit none
   integer        :: u, maxpress = 0
   integer        :: ii
   character(100) :: test_name
   npress = 0
   nisothermpure: do ii=1,ncomponents
    write (test_name, '( "isoterma", I1, ".dat" )' ) ii
    call system ("cp " // test_name // " isotermaN.dat")
    open(newunit=u,file="isotermaN.dat",status='old',iostat=err_apertura)
    if(err_apertura/=0) stop '[ERROR] isotermaN.dat :  Catastrophic failure'
    do0: do
     READ (u,'(A)',IOSTAT=err_apertura) line
     IF( err_apertura /= 0 ) EXIT do0
     npress(ii)=npress(ii)+1
    end do do0
    close(u)
   end do nisothermpure
   u = maxval(npress)
   nisothermpure1: do ii=1,ncomponents
    write (test_name, '( "isoterma", I1, ".dat" )' ) ii
    write(6,'(a15,1x,a15)') test_name,ajuste(ii)
    call system ("cp " // test_name // " isotermaN.dat")
    open(unit=456,file="isotermaN.dat",status='old',iostat=err_apertura)
    if(err_apertura/=0) stop '[ERROR] isotermaN.dat :  Catastrophic failure, check maxnp, maxcompounds or maxdata'
    do1: do i = 1,npress(ii)
     read (456,'(A)',iostat=err_apertura) line
     if(err_apertura/=0) exit do1
     read(line,*)datas(1,ii,i),datas(2,ii,i)
     write(6,'(2f25.10)') datas(1,ii,i),datas(2,ii,i)
    end do do1
    close(456)
   end do nisothermpure1
   end subroutine ReadIsotherms
  ! ...
   subroutine cite()
    implicit none
    integer           :: iii
    character(len=80) :: funk
    write(6,'(a)')     '====================================='
    write(6,'(a)')     'If you are using GAIAST and would like to cite it, then for GAIAST: '
    write(6,'(a)')     'GAIAST: 1. Salvador Rodríguez-Gómez Balestra, Rocio Bueno-Perez, Sofia Calero (2016) GAIAST, Zenodo'
    write(6,'(a)')     '        doi: 10.5281/zenodo.165844'
    write(6,'(a)')     'For the IAST theory:'
    write(6,'(a)')     '        A. L. Myers and J. M. Prausnitz, AIChE J., Thermodynamics of mixed-gas adsorption, 1965, 11, 121.'
    write(6,'(a)')     'For the isotherms models:'
    c01234: do iii=1,ncomponents
     funk = ajuste(iii)
     c01235: select case (funk)
      case ("langmuir")
       write(6,'(a,1x,i3,a,1x,a)')&
       "Model",iii,":","I. Langmuir, J. Am. Chem. Soc. 38 (11) (1916) 2221–2295."
      case ("langmuir_dualsite")
       write(6,'(a,1x,i3,a,1x,a)')&
       "Model",iii,":","Myers, A. L. (1983). AIChE journal, 29(4), 691-693, doi: 10.1002/aic.690290428"
      case ("toth")
       write(6,'(a,1x,i3,a,1x,a)')"Model",iii,":","Toth, Acta Chem. Acad. Hung. 69 (1971) 311–317."
      case ("jensen_seaton")
       write(6,'(a,1x,i3,a,1x,a)')"Model",iii,":","C. R. C. Jensen and N. A. Seaton, Langmuir, 1996, 12, 2866."
      case ("langmuir_freundlich")
   ! Also known as the Sips equation
       write(6,'(a,1x,i3)')'Model',iii
       write(6,*)'R. Sips, J. Chem. Phys., (1948); http://dx.doi.org/10.1063/1.1746922'
       write(6,*)'R. Sips, J. Chem. Phys. 18, 1024 (1950); http://dx.doi.org/10.1063/1.1747848'
       !write(6,*)'Turiel et al., 2003, DOI: 10.1039/B210712K'
       !write(6,*)'Umpleby, R. J., Baxter, S. C., Chen, Y., Shah, R. N., & Shimizu, K. D. (2001)., 73(19), 4584-4591.'
      case ("langmuir_sips")
       write(6,'(a,1x,i3)')'Model',iii
       write(6,'(8x,a)')"I. Langmuir, J. Am. Chem. Soc., 1916, 38(11), 2221–2295."
       write(6,'(8x,a)')"A. L. Myers, AIChE journal, 1983, 29(4), 691-693, doi: 10.1002/aic.690290428"
       write(6,'(8x,a)')'R. Sips, J. Chem. Phys., 1948 ; doi:http://dx.doi.org/10.1063/1.1746922'
      case("jovanovic")
       write(6,'(a,1x,i3,a,1x,a)')&
       "Model",iii,":","Jovanovic, D. S., Kolloid-Z.Z. Polym. 235, 1203 ( 1969 )"
      case("jovanovic_freundlich")
       write(6,'(a,1x,i3,a,1x,a)')&
        "Model",iii,":",'10.1006/jcis.1996.0518 and 10.1016/S0021-9673(97)01096-0'
     end select c01235
    end do c01234
   end subroutine cite
  !
   subroutine MakeInitPOP(funk,n)
   implicit none
   character(100),intent(in)::  funk
   integer,intent(out)::        n
   select case (funk)
    case("freundlich")
     n=2 ! numero de parametros del modelo.
    case ("langmuir")
     n=2
    case ("redlich_peterson")
     n=3
    case ("langmuir_freundlich","lf")
     n=3
    case ("toth")
     n=3
    case ("jensen_seaton")
     n=4
    case ("dubinin_raduschkevich")
     n=3
    case ("langmuir_dualsite")
     n=4
    case ("langmuir_freundlich_dualsite",'lf2')
     n=6
    case ("redlich_peterson_dualsite")
     n=6
    case ("langmuir_freundlich_3order","lf3")
     n=9
    case ("dubinin_astakhov")
     n=4
    case ("jovanovic_smoothed")
     n=3
    case ("jovanovic")
     n=2
    case ("jovanovic_freundlich")
     n=4
    case ("langmuir_sips")
     n=5
    case ("isobare")
     n = 4
    case ("asymtotic_temkin","cory")
     n = 3
    case ("Brunauer-Emmett-Teller","BET")
     n = 3
    case ("Quadratic")
     n = 3
   end select
   return
   end subroutine MakeInitPOP
  
   pure real function model(a,n,xx,funk,equation)
  ! Availabel models:
  ! ------------------------------------------------
  ! freundlich
  ! langmuir
  ! langmuir_freundlich
  ! redlich_peterson
  ! redlich_peterson_dualsite
  ! toth
  ! langmuir_dualsite
  ! langmuir_freundlich_dualsite
  ! langmuir_freundlich_3order
  ! dubinin_raduschkevich
  ! dubinin_astakhov
  ! jovanovic_smoothed
  ! jovanovic
  ! jovanovic_freundlich
  ! langmuir_sips
  ! jensen_seaton
  ! --------------------------------------------------
    implicit none
    integer,intent(in)                 :: n
    real,intent(in)                    :: a(0:n-1)
    real,intent(in)                    :: xx
    character(100),intent(in)          :: funk
    character(100),intent(in),optional :: equation
    select case (funk)
     case("freundlich","f")
      model = a(0)*xx**a(1)
     case("langmuir","l")
      model = a(0)*a(1)*xx/(1+a(1)*xx)
     case("langmuir_freundlich","lf")
      model = a(0)*a(1)*xx**a(2)/(1.0+a(1)*xx**a(2))
     case("redlich_peterson")
      model= a(0)*a(1)*xx/(1.0+a(1)*xx**a(2))
     case("redlich_peterson_dualsite")
      model = a(0)*a(1)*xx/(1+a(1)*xx**a(2))+a(3)*a(4)*xx/(1.0+a(4)*xx**a(5))
     case("toth","t") ! f(x)=Nmax*alfa*x/(1+(alfa*x)**c)**(1/c)
      model = a(0)*a(1)*xx/((1+(a(1)*xx)**(a(2)))**(1.0/a(2)))
     case("langmuir_dualsite","l2")
      model = a(0)*a(1)*xx/(1+a(1)*xx) + a(2)*a(3)*xx/(1.0+a(3)*xx)
     case("langmuir_freundlich_dualsite","lf2")
      model = a(0)*a(1)*xx**a(2)/(1+a(1)*xx**a(2))+a(3)*a(4)*xx**a(5)/(1.0+a(4)*xx**a(5))
     case("langmuir_freundlich_3order","lf3")
      model = a(0)*a(1)*xx**a(2)/(1+a(1)*xx**a(2))+a(3)*a(4)*xx**a(5)/(1.0+a(4)*xx**a(5)) + &
              a(6)*a(7)*xx**a(8)/(1+a(7)*xx**a(8))
     case ("dubinin_raduschkevich","dr") ! N=Nm*exp(-(RT/Eo ln(Po/P))^2)  #model
      model = a(0)*exp(-((R*T/a(1))*log(a(2)/xx) )**2)
     case ("dubinin_astakhov","da")       ! N=Nm*exp(-(RT/Eo ln(Po/P))^d) #model
      model = a(0)*exp(-((R*T/a(1))*log(a(2)/xx) )**a(3))
     case ("jovanovic_smoothed","js")
      model = a(0)*(1.0 - exp(-a(1)*xx))*exp(-a(2)*xx)
     case ("jovanovic","j")
      model = a(0)*(1.0 - exp(-a(1)*xx))
     case ("jovanovic_freundlich","jf")
      model = a(0)*(1.0 - exp(-(a(1)*xx)**a(2)))*exp(a(3)*xx**a(2))
     case ("langmuir_sips","ls")
      model = (a(0)*a(1)*xx)/(1+a(1)*xx) + (a(2)*a(3)*xx**a(4))/(1+a(3)*xx**a(4))
     case ("jensen_seaton","jsn")
      model = a(0)*xx*( 1.0 + ( a(0)*xx / (a(1)*( 1+a(2)*xx )) )**a(3) )**(-1.0/a(3))
     case ("asymtotic_temkin","cory")
      ! Asymptotic approximation to the Temkin Isotherm, (see DOI: 10.1039/C3CP55039G)
      model = a(0)*(a(1)*xx/(1+a(1)*xx) + a(2)*((a(1)*xx/(1+a(1)*xx))**2)*(a(1)*xx/(1+a(1)*xx)-1))
     case ("Brunauer-Emmett-Teller","BET")
      model = a(0)*( a(1)*xx/((1-a(2)*xx)*(1-a(2)*xx+a(1)*xx)))
     case ("Quadratic")
      model = a(0)*xx*(a(1)+2*a(2)*xx)/(1+a(1)*xx+a(2)*xx*xx)
     case ("UserDefined")
      model = 0.0
     case ("isobare")
      model = 1.0/(a(0)*exp(a(1) - a(2)/xx) + a(3))
    end select
    return
   end function model
  !
  end module gaiast_globals
  !
  program main
   use mod_random
   use general
   use gaiast_globals
   use mod_genetic
   use mod_simplex
   use, intrinsic :: iso_fortran_env
   call clear_screen()
   print '(4a)', 'This file was compiled by ', &
         compiler_version(), ' using the options ', &
         compiler_options()
   print '(a)',' ',&
  "#     ________    _____   .___    _____     ____________________  ",&
  "#    /  _____/   /  _  \  |   |  /  _  \   /   _____/\__    ___/  ",&
  "#   /   \  ___  /  /_\  \ |   | /  /_\  \  \_____  \   |    |     ",&
  "#   \    \_\  \/    |    \|   |/    |    \ /        \  |    |     ",&
  "#    \______  /\____|__  /|___|\____|__  //_______  /  |____|     ",&
  "#           \/         \/              \/         \/              ",&
  "# flama, cabesa.",&
  " "
   call read_input()
   if (seed_flag) then
    call init_random_seed( )
    seed_flag=.false.
   end if
   call ReadIsotherms()
   if(flag)then
    open(111,file="iso.dat")
    do compound = 1, ncomponents
     write(6,'(a)')' '
     write(6,'("Fitting compound:",1x,i2)') compound
     write(6,'(a14,1x,a24,1x,a24,10x,a14)')'Compound/Step:','Cromosome:','Parameters','Control'
     call fit(compound)
     call fit_simplex()
    end do
   end if
   if (IAST_flag) then
    call bounds()
    if(ncomponents>2.or.Always_Use_Multicomponent_Solver)then
     ! Always_Use_Multicomponent_Solver defined in the header of gaiast_globals module
     call IAST_multicomponent()
    else
     ! It works fine
     call IAST_binary()
    end if
   end if
   call cite()
   close(111)
   write(6,'(a)')'====================================='
   stop ':)'
  end program main