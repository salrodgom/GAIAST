! GAIAST
! Salvador R.G. Balestra, Rocio Bueno-Perez, and Sofia Calero
! Wednesday, 09, November, 2016
! Wednesday, 13, Jun, 2018 ( add SIMPLEX optimisation)
module mod_random
! module for pseudo random numbers
 implicit none
 private
 public init_random_seed, randint, r4_uniform
 contains
 subroutine init_random_seed( )
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer              :: i, n, un, istat, dt(8), pid
  integer(int64)       :: t
  call random_seed(size = n)
  allocate(seed(n))
  write(6,'(a)')"Random seed:"
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
   read(un) seed
   close(un)
   write(6,'(a)')"OS provides a random number generator"
  else
   ! Fallback to XOR:ing the current time and pid. The PID is
   ! useful in case one launches multiple instances of the same
   ! program in parallel.
   call system_clock(t)
   if (t == 0) then
      call date_and_time(values=dt)
      t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
           + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
           + dt(3) * 24_int64 * 60 * 60 * 1000 &
           + dt(5) * 60 * 60 * 1000 &
           + dt(6) * 60 * 1000 + dt(7) * 1000 &
           + dt(8)
   end if
   pid = getpid()
   t = ieor(t, int(pid, kind(t)))
   do i = 1, n
      seed(i) = lcg(t)
   end do
   write(6,'(a)')"Fallback to the current time and pid."
  end if
  call random_seed(put=seed)
  write(6,*) seed
 contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  pure integer function lcg(s)
    integer(int64),intent(in) :: s
    integer(int64)            :: ss
    if (s == 0) then
     ss = 104729
    else
     ss = mod(s, 4294967296_int64)
    end if
    ss = mod(ss * 279470273_int64, 4294967291_int64)
    lcg = int(mod(ss, int(huge(0), int64)), kind(0))
  end function lcg
 end subroutine init_random_seed
!
 integer function randint(i,j)
  integer,intent(in)    :: i,j
  real                  :: r
  call random_number(r)
  randint = i + floor((j+1-i)*r)
 end function randint
!
 real function r4_uniform(x,y)
  implicit none
  real,intent(in)       :: x,y
  real                  :: r
  call random_number(r)
  r4_uniform=(r*(y-x))+x
  return
 end function r4_uniform
end module mod_random
!
module qsort_c_module
! Recursive Fortran 95 quicksort routine sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms, 1997 printing
! Made F conformant by Walt Brainerd
 implicit none
 public  :: QsortC
 private :: Partition
 contains
!
 recursive subroutine QsortC(A)
  real(16), intent(in out), dimension(:) :: A
  integer                            :: iq
  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
 end subroutine QsortC
!
 subroutine Partition(A, marker)
  real(16), intent(in out), dimension(:) :: A
  integer, intent(out)               :: marker
  integer                            :: i, j
  real(16)                           :: temp
  real(16)                           :: x
  x = A(1)
  i= 0
  j= size(A) + 1
  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do
 end subroutine Partition
end module qsort_c_module

module gaiast_globals
 use mod_random
 implicit none
 integer                    :: npar,ncomponents,i
 integer(8)                 :: seed = 0
 integer,allocatable        :: np(:)
 integer                    :: err_apertura,ii,intervalos,j
 integer,parameter          :: integration_points = 10000
 real,parameter             :: precision_Newton = 1e-5
 real                       :: tol = 0.001, tolfire = 0.25
 real,allocatable           :: param(:,:),concy(:),pi(:,:),iso(:,:)
 integer,allocatable        :: npress(:)
 integer                    :: compound
 integer,parameter          :: maxnp = 10, maxcompounds = 10, maxdata = 1000
 real,target                :: datas(2,maxcompounds,maxdata),f(maxdata)
 real,pointer               :: x(:),y(:),alldat(:,:)
 character(100),allocatable      :: ajuste(:)
 character(32*maxnp),allocatable :: string_IEEE(:)
 character(100)             :: line,string,intmethod,Equation
 character(5)               :: inpt
 logical                    :: flag = .true., FlagFire = .false.,seed_flag=.true.
 logical                    :: physical_constrains = .false., range_flag =.true.
 logical                    :: Always_Use_Multicomponent_Solver=.true.
 logical                    :: refit_flag = .false.
 real,parameter             :: R = 0.008314472 ! kJ / mol / K
 real                       :: T = 298.0
 real                       :: inferior
 real                       :: superior
 contains
!
  pure character(len=32) function real2bin(x)
   implicit none
   real,intent(in) :: x
   write(real2bin(1:32) ,'(b32.32)') x
   return
  end function real2bin
!
  pure real function bin2real(a)
   implicit none
   character(len=32), intent(in) :: a
   read(a(1:32),'(b32.32)') bin2real
   return
  end function bin2real
!
  subroutine read_input()
  implicit none
  integer          :: ii,i,j
  character(len=4) :: realorbinary
  read_input_do: do
   read(5,'(A)',iostat=err_apertura)line
   if ( err_apertura /= 0 ) exit read_input_do
   if(line(1:1)=='#') cycle read_input_do
   if(line(1:3)=="end".or.line(1:)=="end_iast_input") exit read_input_do
   if(line(1:5)=='ncomp')then
    read(line,*)inpt,ncomponents
    allocate(ajuste(ncomponents),concy(ncomponents))
    allocate(npress(ncomponents),np(ncomponents))
    allocate(string_IEEE(ncomponents))
    do ii=1,ncomponents
     read(5,*)concy(ii),ajuste(ii)
     string = ajuste(ii)
     call MakeInitPOP(string,npar)
     np(ii) = npar
     write(6,'(a,1x,i3,1x,a,1x,i3,1x,a)')'compound:',ii,'with',np(ii),'parameters'
    end do
    concy=concy/sum(concy)
    allocate(param(ncomponents,0:maxval(np)-1))
    param = 0.0
    cycle read_input_do
   end if
   if(line(1:5)=='ffit?')then
    read(line,*)inpt,flag
    if(flag.eqv..false.)then
     write(6,*)'[WARN] Fit is already done'
     readparameter: do ii=1,ncomponents
      read(5,*)(param(ii,j),j=0,np(ii)-1)
     end do readparameter
    end if
    !cycle read_input_do
   end if
   if(line(1:5)=='refit') then
    write(6,*)'[WARN] Refitting parameters'
    read(line,*)inpt,refit_flag,realorbinary
    if(refit_flag)then
     select case(realorbinary)
     case('real')
      do ii=1,ncomponents
       read(5,*)(param(ii,j),j=0,np(ii)-1)
      end do
      do ii=1,ncomponents
       write(6,'(a,1x,i3,1x,a,1x,i3,1x,a)')'compound:',ii,'with',np(ii),'parameters'
       do j=1,np(ii)
        string_IEEE(ii)(32*(j-1)+1:32*j)=real2bin( param( ii,j-1 ) )
        write(6,'(i3,1x,i3,1x,f14.7,1x,a32)') ii,j,param( ii,j-1 ), string_IEEE(ii)(32*(j-1)+1:32*j)
       end do
      end do
     case default
      do ii=1,ncomponents
       write(6,'(a,1x,i3,1x,a,1x,i3,1x,a)')'compound:',ii,'with',np(ii),'parameters'
       do j=1,np(ii)
        read(5,'(a32)') string_IEEE(ii)(32*(j-1)+1:32*j) 
        param( ii,j-1 )=bin2real( string_IEEE(ii)(32*(j-1)+1:32*j) )
        write(6,'(i3,1x,i3,1x,f14.7,1x,a32)') ii,j,param( ii,j-1 ), string_IEEE(ii)(32*(j-1)+1:32*j)
       end do
      end do
     end select
     write(6,'(a)') string_IEEE(1:ncomponents)
    else
     do ii=1,ncomponents
      string_IEEE(ii)=' '
     end do
    end if
   end if
   if(line(1:5)=='inter') then
     read(line,*)inpt,intervalos
     allocate(pi(ncomponents,0:intervalos),iso(ncomponents,0:intervalos))
   end if
   if(line(1:5)=='toler') read(line,*)inpt, Tol
   if(line(1:5)=='Range')then
    read(line,*)inpt, inferior, superior
    range_flag = .false.
    write(6,'(a)')'[WARN] Fugacity range specified'
   end if
   if(line(1:5)=='RSeed') then
    seed_flag=.false.
    read(line,*)inpt, seed
   end if
   if(line(1:5)=='tempe') read(line,*)inpt, T
   if(line(1:5)=='InteM') read(line,*)inpt, IntMethod
   if(line(1:22)=='physically_constrained') then
    physical_constrains=.true.
    write(6,'(a)') '[WARN] The fits are physically constrained'
   end if
   if(line(1:5)=='Fire?') then
     read(line,*) inpt, FlagFire
     if(FlagFire) read(5,*) TolFire
   end if
   if(err_apertura/=0) exit read_input_do
  end do read_input_do
 end subroutine read_input
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
 subroutine IAST_binary()
  implicit none
  integer           :: i,j,k,ijk
  real              :: dx,p,concx(ncomponents),aux1,aux2
  real              :: min_tolerance
  integer           :: j_id
  real, allocatable :: j_tolerance(:)
  real              :: presion(ncomponents),n(0:ncomponents)
  dx = real((superior-inferior)/intervalos)
  allocate( j_tolerance(0:intervalos-1 ))
  open(unit=104,file='adsorcion.dat')
  do i=0,intervalos-1
   do j=0,intervalos-1
    if( abs( pi(1,i) - pi(2,j)) < tol ) then
   !j_tolerance = 0.0
   !j_tolerance(0:intervalos-1) = abs( pi(1,i) - pi(2,0:intervalos-1) )
   !min_tolerance = minval( j_tolerance )
   !j = minloc( j_tolerance,1,(j_tolerance<tol))
   !j = minloc( j_tolerance,1)
   !if( min_tolerance > tol ) STOP '[error] Tolerance value is bigger than the tolerance value by default]'
   ! {{
     presion(1) = exp(inferior+i*dx)
     presion(2) = exp(inferior+j*dx)
! {{ Calculate molar fraction of component 1 at spread pressure Pi from a general formula:
!    x1=1/(1+p1*y2/p2Y1 + w(o3) )
!    x2=1 - x2
! }}
     concx(1) = presion(2)*concy(1) / (concy(2)*presion(1) + concy(1)*presion(2))
     concx(2) = 1.0 - concx(1)   ! }}
     p  = presion(1)*concx(1)/concy(1)
     n(1) = iso(1,i)*concx(1)
     n(2) = iso(2,j)*concx(2)
     n(0) = n(1) + n(2)
     write(104,*)p,(n(ijk),ijk=1,ncomponents),n(0)
   ! }}
    end if
   end do
  end do
  close(104)
  deallocate(j_tolerance)
  return
 end subroutine IAST_binary
! ...
 subroutine IAST_multicomponent()
  implicit none
  integer          :: i,j,ijk,nn
  real,allocatable :: a(:)
  character(100)   :: funk
  real             :: comp,dx,p,concx(ncomponents),pureloading(ncomponents)
  real             :: presion(0:ncomponents,0:intervalos-1),n(0:ncomponents)
  real             :: auxP,aux1,aux2,lastPi(ncomponents),piValue=0.0,x,last=0.0
  logical          :: validPressure = .true., yIsZero = .true.,flag = .true.
  dx = real((superior-inferior)/intervalos)
  i = 0
  lastPi = 0.0
  presion = 0.0
  open(unit=104,file='adsorcion.dat')
  ! Solve IAST for every pressure point
  ScanPressures: do while (i<=intervalos-1)
! {{ initialisation:
    n    = 0.0
    aux1 = 0
    aux2 = 0
    auxP = datas(1,1,npress(1))
    presion(1,i) = exp(inferior+i*dx)
    last    = lastPi(1)                       !<---- works
    piValue = CalculatePi(1,presion,i,last)   !<---- works
    validPressure=CalculatePressures(presion,i,lastPi,piValue)   !<- WRONG! (why??, is it a mistake?)
    if(validPressure.and.i>0)then
     !write(6,*)i,(presion(ijk,i),ijk=1,ncomponents),piValue
     continue
    else
     i=i+1
     cycle ScanPressures
    end if
! }}
    aux1 = 0.0
! {{ Calculate molar fraction of component 1 at spread pressure Pi from a general formula: x1=1/(1+P1Y2/P2Y1+P1Y3/P3Y1+...)
    do ijk=1,ncomponents
     if (presion(ijk,i) /= 0.0) aux1=aux1+presion(1,i)*concy(ijk)/(presion(ijk,i)*concy(1))
    end do
! Case where there is not a Pressure i that makes Pi of component i = Pi of component 1: aux1 is zero by initialization
    if (aux1/=0.0) then
     concx(1) = 1.0/(1+aux1)
    else
     concx(1) = 0.0
    end if
! Special case when all y are zero apart from component 1
    yIsZero = .true.
    do ijk=2,ncomponents
     if (concy(ijk)/=0.0) yIsZero =.false.
    end do
    if (yIsZero) concx(1) = 1.0
! }}
!! {{ Partial Pressure from component 1:
    p = presion(1,i)*concx(1)/concy(1)
!! Calculate the molar fraction of the different components
    do ijk=2,ncomponents
     if (presion(ijk,i)/=0.0) concx(ijk)=(p*concy(ijk))/presion(ijk,i)
    end do
! {{ renormalise X:
    aux1 = sum( concx )
    concx = concx / aux1
! }}
! Calculate the total loading
    flag=.true.
    do ijk=1,ncomponents
!   Calculate pure loading:
     funk = ajuste(ijk)
     nn = np(ijk)
     allocate(a(0:nn-1))
     a = 0.0
     do j = 0, nn-1
      a(j) = param(ijk,j)
     end do
     x = presion(ijk,i)
     pureloading(ijk) = model(a,nn,x,funk)
     deallocate(a)
     if(pureloading(ijk)==0.0) flag=.false.
    end do
    if(flag)then
     do ijk=1,ncomponents
      if(pureloading(ijk)/=0.0)aux2=aux2+concx(ijk)/pureloading(ijk)
     end do
    else
     aux2 = 0.0
    end if
    if (aux2/=0.0) n(0) = 1.0/aux2
    do ijk=1,ncomponents
     n(ijk) = pureloading(ijk)*concx(ijk)
    end do
    write(104,*)p,(n(ijk),ijk=1,ncomponents),n(0) !,(pureloading(ijk),ijk=1,ncomponents),i
    i = i +1
  end do ScanPressures
  close(104)
  return
 end subroutine IAST_multicomponent
! ...
 logical function CalculatePressures(presion,point,lastPi,piValue)
  implicit none
  integer,intent(in) :: point
  real,intent(inout) :: presion(0:ncomponents,0:intervalos-1),piValue
  real,intent(inout) :: lastPi(ncomponents)
  real               :: oldPi, p,a1,a2,b1,b2,c,a,p1,p2
  integer            :: k,j,n
  real,allocatable   :: apar(:)
  character(100)     :: funk
  character(10)      :: solve_method = 'integral' ! OR "analytical"
  CalculatePressures = .true.
  k=1
  calc: do while (k<=ncomponents.and.CalculatePressures)
! {{ 
!  Interpolation and parameters from k-compound
   funk = ajuste( k )
   !if(funk=='langmuir') solve_method = 'analytical' !.or.funk=='langmuir_dualsite') solve_method = 'analytical'
   n = np(k)
   allocate(apar(0:n-1))
   do j = 0, n-1
    apar(j) = param(k,j)
   end do
! }}
   select case (solve_method)
    case ('analytical')
     select case (funk)
      case ('langmuir')
      presion(k,point) = (exp(piValue/apar(0))-1/apar(1))
      !case('langmuir_dualsite')
     end select
    case ('integral')
     oldPi = lastPi(k)
     apar = 0.0
     do j = 0, n-1
      apar(j) = param(k,j)
     end do
     if (point==0) then
      p = 0.0
      presion(k,point) = CalculatePressureIntegral(k,p,apar,n,funk,oldPi,piValue)
     else
      p = presion(k,point-1)
      presion(k,point) = CalculatePressureIntegral(k,p,apar,n,funk,oldPi,piValue)
     end if
    case default
     presion(k,point) = CalculatePressureLinear(k,piValue)
   end select
   lastPi(k) = piValue
   if ( point>0.and.presion(k,point) < presion(k,point-1)) then
    do j=1,ncomponents
     presion(j,point) = 0.0
    end do
    CalculatePressures = .false.
    deallocate(apar)
    cycle calc
   else
    CalculatePressures = .true.
   end if
   k = k +1
   deallocate(apar)
  end do calc
! }}
  return
 end function CalculatePressures
! ...
 real function CalculatePressureLinear(compound,piValue)
  implicit none
  integer        :: j,compound
  logical        :: lastPass = .false.
  real           :: piGuess=0.0,deltaPi = 0.0
  real           :: piValue
  CalculatePressureLinear = 0.0
  do j=0,npress(compound)
   if(j==0)then
     deltaPi=datas(2,compound,j+1) ! isother(k,j+1,1)
   else
     deltaPi=(datas(2,compound,j+1)-datas(2,compound,j))+(datas(1,compound,j+1)*datas(2,compound,j) &
            - datas(1,compound,j)*datas(2,compound,j+1))*log(datas(1,compound,j+1)/datas(1,compound,j))/&
              (datas(1,compound,j+1)-datas(1,compound,j))
   end if
   if(piGuess+deltaPi>=piValue)then
    if(j==0)then
     CalculatePressureLinear=piValue*datas(1,compound,j+1)/datas(2,compound,j+1)
    else
     CalculatePressureLinear=SolveNewtonLinear(compound,piValue-piGuess,j)
    end if
    exit
   else
    piGuess=piGuess+deltaPi
    if(j==npress(compound)-1) lastPass = .true.
   end if
  end do
  if(lastPass.and.piGuess<piValue) CalculatePressureLinear = 0.0
  return
 end function CalculatePressureLinear
!
 real function SolveNewtonLinear(compound,piValue,n)
  implicit none
  integer,intent(in) :: compound,n
  real,intent(in)    :: piValue
  real               :: X,X0,Y
  real               :: gradient
  real               :: a,b,c
  ! Calculate the coefficients of the fuction to find the zero
  a=(datas(2,compound,n+1)-datas(2,compound,n))/(datas(1,compound,n+1)-datas(1,compound,n))
  b=(datas(1,compound,n+1)*datas(2,compound,n)-datas(1,compound,n)*datas(2,compound,n+1))/ &
    (datas(1,compound,n+1)-datas(1,compound,n))
  c=-a*datas(1,compound,n)-b*log(datas(1,compound,n))-piValue
  ! Start looking for a solution at the middle of the interval [P1,P2]
  X=(datas(1,compound,n+1)+datas(1,compound,n))*0.5
  ! Calculate the value of the function at point X
 Y=a*X+b*log(X)+c
 do while(abs(Y)>precision_newton)
  ! Calculate the gradient in the point X
  gradient = 0.0
  gradient = a+b/X
  if (gradient/=0) then
   ! Calculate the equation of the tangent
   X0 = Y-gradient*X
   ! Calculate the 0 point of the tangent and make it our new X
   X = -X0/gradient
   ! Calculate the value of the function at point X
   Y = a*X+b*log(X)+c
  elseif (X > datas(1,compound,n+1)) then
   X=X-precision_Newton
  else
   X=X+precision_Newton
  end if
 end do
 if(X<datas(1,compound,n).or.X>datas(1,compound,n+1)) then
  write(6,*)"WARNING: A solution for the linear equation was not found between pressures",&
   datas(1,compound,n),'and',datas(1,compound,n+1),"."
  write(6,*)"No pressure will be considered in this interval"
  SolveNewtonLinear = 0.0
 else
  SolveNewtonLinear = X
 end if
 return
 end function SolveNewtonLinear

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

 real function CalculatePressureIntegral(k,x0,a,n,funk,oldPi,piValue)
  implicit none
  integer,intent(in)       :: k,n
  integer                  :: i,imax
  real                     :: x, delta, integral
  real,intent(in)          :: a(0:n-1), x0,oldPi,piValue
  character(100),intent(in):: funk
  character(100)           :: mode = 'version1'
  integral = 0.0
  x = 0.0
  i = 0
  imax = 10000000
  if ( x0 == 0.0 ) then
   delta = 1.0/real(integration_points)
  else
   delta = x0/real(integration_points)!+x0)
  end if
  pressint: do while ( oldPi+ integral <= piValue )
  select case (mode)
   case ('version1')
    !Integrate by the paralelogram rule
    x=x0+delta*(0.5+i)
    Integral = Integral + delta*model(a,n,x,funk)/x
    i = i+1
   case ('version2')
    x=x0+delta*0.5*(i*i+2.0*i+1.0)
    Integral = Integral + delta*model(a,n,x,funk)/x
  end select
  !
  if(oldPi + Integral >= piValue)then
   if((oldPi + Integral)/piValue >= 1.1) then
    write(6,*)'Integal greater than Pi!',(oldPi + Integral)/piValue
   end if
  end if
   if (i==imax) then
     write(6,*)'Overflow in the integral. Compound:',k,piValue,oldPi + Integral
     exit pressint
   end if
  end do pressint
  CalculatePressureIntegral = x
  !write(6,*)'[CalculatePressureIntegral]',i,k,CalculatePressureIntegral,oldPi,piValue,x0,delta
  !write(6,*)'[ParametersModel]',(a(i),i=0,n-1)
  return
 end function CalculatePressureIntegral

 real function CalculatePi(k,presion,i,lastPi)
  implicit none
  integer     :: i,k
  real        :: lastPi
  integer            :: n
  real,allocatable   :: apar(:)
  real               :: x0,x1,presion(0:ncomponents,0:intervalos-1)
  character(100)     :: funk,model
  funk = ajuste(k)
  n = np(k)
  allocate(apar(0:n-1))
  apar = 0.0
  do j = 0, n-1
   apar(j) = param(k,j)
  end do
  model = ajuste(k)
  select case (model)
  case ('langmuir')
   x1 = presion(k,i)
   lastPi = apar(0)*log(1+apar(1)*x1)
  case ('langmuir_dualsite')
   x1 = presion(k,i)
   lastPi = apar(0)*log(1+apar(1)*x1)+apar(2)*log(1+apar(3)*x1)
  case ('jovanovic')
   x1 = presion(k,i)
   lastPi = apar(0)*(x1+exp(-apar(1)*x1)/apar(1))
  case default
   if ( i==0 ) then
    x0=0.0
    x1=presion(k,i)
    lastPi = lastPi + integrate(x0,x1,apar,n,funk)
   else
    x0=presion(k,i-1)
    x1=presion(k,i)
    lastPi = lastPi + integrate(x0,x1,apar,n,funk)
   end if
  end select
  deallocate(apar)
  CalculatePi = lastPi
  return
 end function

 subroutine ReadIsotherms()
 implicit none
 integer        :: iso = 456, maxpress = 0
 character(100) :: test_name
 npress = 0
 nisothermpure: do ii=1,ncomponents
  write (test_name, '( "isoterma", I1, ".dat" )' ) ii
  call system ("cp " // test_name // " isotermaN.dat")
  open(unit=iso,file="isotermaN.dat",status='old',iostat=err_apertura)
  if(err_apertura/=0) stop '[ERROR] isotermaN.dat :  Catastrophic failure'
  do0: do
   READ (iso,'(A)',IOSTAT=err_apertura) line
   IF( err_apertura /= 0 ) EXIT do0
   npress(ii)=npress(ii)+1
  end do do0
  close(iso)
 end do nisothermpure
 iso = maxval(npress)
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
   write(6,'(2f20.5)') datas(1,ii,i),datas(2,ii,i)
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
  case ("langmuir_freundlich")
   n=3
  case ("toth")
   n=3
  case ("jensen_seaton")
   n=4
  case ("dubinin_raduschkevich")
   n=3
  case ("langmuir_dualsite")
   n=4
  case ("langmuir_freundlich_dualsite")
   n=6
  case ("redlich_peterson_dualsite")
   n=6
  case ("langmuir_freundlich_3order")
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
   n=4
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
   case ("UserDefined")
    model = 0.0
   case ("isobare")
    model = 1.0/(a(0)*exp(a(1) - a(2)/xx) + a(3))
  end select
  return
 end function model
!
 subroutine clear_screen
  write(*,'(a1,a)',advance="no") achar(27), '[2J'
 end subroutine clear_screen
end module gaiast_globals

module mod_genetic
 use mod_random
 use gaiast_globals
 use qsort_c_module
 implicit none
 !private
 !public :: fit
 integer,parameter             :: GA_size     = 2**12 ! numero de cromosomas
 real,parameter                :: GA_MutationRate = 0.3333 !2000/real(ga_size) ! ga_mutationrate=0.333
 real,parameter                :: GA_EliteRate= 0.25, GA_MotleyCrowdRate=0.25,GA_DisasterRate = 0.0000001
 integer,parameter             :: maxlinelength=maxnp*32
 integer,parameter             :: ga_Elitists = int( ga_size * ga_eliterate)
 integer,parameter             :: ga_Motleists = int( ga_size - ga_size*GA_MotleyCrowdRate)
 type                          :: typ_ga
  character(len=maxlinelength) :: genotype
  real                         :: phenotype(1:maxnp)
  real                         :: fitness
 end type
 type(typ_ga), pointer         :: parents(:)
 type(typ_ga), pointer         :: children(:)
 type(typ_ga), target          :: pop_alpha( ga_size )
 type(typ_ga), target          :: pop_beta( ga_size )
 contains
!
  type(typ_ga) function new_citizen(compound)
   implicit none
   integer     :: i,compound,j,k
   new_citizen%genotype = ' '
   do i = 1,32*np(compound)
    new_citizen%genotype(i:i) = achar(randint(48,49))
   end do
   do i = 1,np(compound)
    new_citizen%phenotype(i) = bin2real( new_citizen%genotype(32*(i-1)+1:32*i) )
    !read(new_citizen%genotype(32*(i-1)+1:32*i),'(b32.32)') new_citizen%phenotype(i)
   end do
   new_citizen%fitness = fitness( new_citizen%phenotype,compound)
   return
  end function new_citizen
!
  subroutine UpdateCitizen( axolotl ,compound )
   implicit none
   integer                     :: compound,i,GA_ELITISTS
   real                        :: infinite = 0.0
   type(typ_ga), intent(inout) :: axolotl
   do i = 1,np(compound)
    read(axolotl%genotype(32*(i-1)+1:32*i),'(b32.32)') axolotl%phenotype(i)
   end do
   axolotl%fitness = fitness( axolotl%phenotype,compound)
   return
  end subroutine UpdateCitizen
! 
  pure real function Fitness(phenotype,compound)
   implicit none
   real, intent(in)    :: phenotype(maxnp)
   integer,intent(in)  :: compound
   integer             :: i
   real                :: a(0:np(compound)-1),xx,yy,penalty
   character(len=100)  :: funk
   real,parameter      :: infinite = 50.0
   funk = ajuste(compound)
   do i = 0,np(compound)-1
    a(i) = phenotype(i+1)
   end do
   phys_constrains: if ( physical_constrains ) then
     select case (funk)
      case ('langmuir')
       if(a(0)<0.0.or.a(1)<=0.0)then
        ! constrains:
        ! a>0, b>0
        penalty = infinite
       else
        penalty = 0.0
       end if
      case ('langmuir_dualsite')
       if(a(0)<0.0.or.a(1)<0.0 .or.&
          a(2)<0.0.or.a(3)<0.0 ) then
        ! constrains:
        ! a>0, b>0, c>0, d>0
        penalty = infinite
       else
        penalty = 0.0
       end if
      case ('toth')
       if(a(0)<0.0.or.a(1)<0.0.or.a(2)<0.or.a(2)>1.0)then
        ! constrains:
        ! a>0 ; b>0 ; 0 < c < 1
        penalty = infinite
       else
        penalty = 0.0
       end if
      case ('langmuir_freundlich')
       if(a(0)<0.0.or.a(1)<0.0.or.a(2)<0.or.a(2)>1.0)then
        ! constrains:
        ! a>0 ; b>0 ; 0 < c < 1
        penalty = infinite
       else
        penalty = 0.0
       end if
      case ('jovanovic_freundlich')
       if(a(0)<0.0.or.a(2)<0.or.a(2)>1.0)then
        ! constrains:
        ! a>0 ; 0 < c < 1
        penalty = infinite
       else
        penalty = 0.0
       end if
      case ("redlich_peterson_dualsite")
       if(a(0)<0.0.or.a(1)<0.0.or.a(2)<0.or.&
          a(3)<0.0.or.a(4)<0.0.or.a(5)<0)then
        penalty = infinite
       else
        penalty = 0.0
       end if
      case ('langmuir_freundlich_dualsite')
       if(a(0)<0.0.or.a(1)<0.0.or.a(2)<0.or.a(2)>1.0 .or.&
          a(3)<0.0.or.a(4)<0.0.or.a(5)<0.or.a(5)>1.0 )then
       !if(a(0)<0.0.or.a(1)<0.0.or.a(2)<0.or.&
       !   a(3)<0.0.or.a(4)<0.0.or.a(5)<0 )then
        ! constrains:
        ! a>0 ; b>0 ; 0 < c < 1
        penalty = infinite
       else
        penalty = 0.0
       end if
      case ('langmuir_freundlich_3order')
       if(a(0)<0.0.or.a(1)<0.0.or.a(2)<0.or.&
          a(3)<0.0.or.a(4)<0.0.or.a(5)<0.or.&
          a(6)<0.0.or.a(7)<0.0.or.a(8)<0)then
        penalty = infinite
       else
        penalty = 0.0
       end if
      case ('jensen_seaton')
       if( a(0)<0 .or. a(1)<0 .or. a(2)<0.or.a(3)<0 ) then
        penalty = infinite
       else
        penalty = 0.0
       end if
      case ('langmuir_sips')
       if( a(0)<0.0.or.a(1)<0.0.or.a(4)<0.or.a(4)>1.0 .or.&
           a(2)<0.0.or.a(3)<0.0 )then
        ! constrains:
        ! a>0 ; b>0 ; 0 < c < 1
        penalty = infinite
       else
        penalty = 0.0
       end if
      case ("isobare")
       if( a(0)<0.0 ) then
        penalty = infinite
       else
        penalty = 0.0
       end if
      case default
       penalty = 0.0
     end select
   end if phys_constrains
   fitness = penalty
   do i = 1, npress(compound)
    xx = datas(1,compound,i)
    yy = datas(2,compound,i)
    fitness = fitness + 0.5*( yy - model(a,np(compound),xx,funk) )**2/&
     log(datas(1,Compound,npress(Compound))-datas(1,Compound,1)+1)
   end do
   if(isnan(fitness)) fitness      = 99999999.999999
   if(fitness==0.0000000000)fitness= 99999999.999999
   return
  end function Fitness

  subroutine WriteCitizen(k,kk,kkk,compound,lod,vgh)
   implicit none
   integer,intent(in)             :: k,kk,compound,lod,vgh,kkk
   integer                        :: i
   character(len=100)             :: fmt_
   character(len=32*np(compound)) :: wnowaste
   real                           :: wnowasteparam(1:32*np(compound)),wfitness
   do i=1,32*np(compound)
    wnowaste(i:i)=' '
   end do
   ! ...
   wnowaste = parents(k)%genotype
   do i=1,np(compound)
    wnowasteparam(i) = parents(k)%phenotype(i)
   end do
   wfitness = parents(k)%fitness
   write(6,'(i2,1x,i5,1x,a32,1x,e25.12,1x,f14.7,1x,a,1x,a,i8,a,i8,a,1x,a)')compound,kk,wnowaste(1:32),&
        wnowasteparam(1),wfitness,'[Fitness]','(',kkk,'/',vgh,')','[Similarity]' !,k
   do i=2,np(compound)
    if(lod>0.and.i==3)then
     write(6,'(9x,a32,1x,e25.12,10x,a,1x,i2,a)')wnowaste(1+32*(i-1):32*i),wnowasteparam(i),'Finishing:',lod,'/10'
    else
     write(6,'(9x,a32,1x,e25.12)')wnowaste(1+32*(i-1):32*i),wnowasteparam(i)
    end if
   end do
  end subroutine WriteCitizen

  subroutine piksrt(n,arr)
  ! Sort real array() of n elements
  implicit none
  integer :: n,j,i = 0
  REAL    :: a,arr(n)
  do j=2,n
    a = arr(j)
    do i=j-1,1,-1
      if (arr(i)<=a) goto 10
      arr(i+1)=arr(i)
    end do
    i=0
10  arr(i+1)=a
  end do
  return
  end subroutine piksrt

  subroutine SortByFitness()
   type(typ_ga)             ::  sorted(1:ga_size)
   integer                  ::  k,i
   real(16)                 ::  ftnss(1:ga_size)
   do k=1,ga_size
    ftnss(k)=dble(parents(k)%fitness)
    if(isnan(parents(k)%fitness)) ftnss(k) = 9999999999.d99
   end do
   call QsortC( ftnss )
   exter:do k=1,ga_size ! <- ordered
    inter:do i=1,ga_size
    if( dble(parents(i)%fitness) == ftnss(k))then
      sorted(k) = parents(i)
      cycle inter
    end if
    end do inter
   end do exter
   parents=sorted
   return
  end subroutine SortByFitness

  integer function Biodiversity( compound, animalito, suma  )
   implicit none
   integer,intent(in)             :: Compound
   type(typ_ga), intent(inout)    :: animalito(1:ga_size)
   integer,intent(out)            :: suma
   integer                        :: i,j,k,cont
   character(len=20)              :: mode = 'Normal'
   logical                        :: flag = .false.
   real                           :: error_d = 1e-3
   select case (mode)
    case('None')
     Biodiversity = 0
    case('Normal')
     suma=0
     Biodiversity = 0
     do k =1,ga_size
      do j=k+1,ga_size
       suma = suma + 1
       if( animalito(k)%genotype(1:32*np(Compound)) == animalito(j)%genotype(1:32*np(Compound)) )then
        Biodiversity = Biodiversity + 1
       end if
      end do
     end do
    case('Superficial')
     Biodiversity = 0.0
     suma=0
     do k = 1, ga_size
      do j = k+1,ga_size
       cont=0
       suma=suma+1
       dbio: do i = 1,np(compound)
        bio: if( abs(animalito(k)%phenotype(i) - animalito(j)%phenotype(i)) <= error_d )then
         cont=cont+1
        end if bio
       end do dbio
       if( cont == np(compound) ) Biodiversity = Biodiversity + 1
      end do
     end do
   end select
   return
  end function Biodiversity

  subroutine Mutate( macrophage , compound )
   implicit none
   type(typ_ga), intent(inout) :: macrophage
   integer                     :: ipos,compound
   do i = 1,np(compound)
    ipos = randint(32*(i-1)+1,32*i)
    macrophage%genotype(ipos:ipos) = achar(randint(48,49))
   end do
   return
  end subroutine Mutate

  subroutine NuclearDisaster(Compound)
   implicit none
   integer,intent(in) ::  Compound
   integer            :: k = 0, i, j
   real               :: rrr
   do i = 2, GA_Size
    children(i) = new_citizen(compound)
   end do
   return
  end subroutine NuclearDisaster

  subroutine Swap()
   if (associated(parents, target=pop_alpha)) then
       parents => pop_beta
       children => pop_alpha
   else
       parents => pop_alpha
       children => pop_beta
   end if
   return
  end subroutine Swap

  subroutine Elitism()
   children(:GA_ELITISTS) = parents(:GA_ELITISTS)
   return
  end subroutine
!
  subroutine crossover(s1,s2,i1,i2,j1,j2)
   implicit none
   integer,intent(in) :: s1,s2,i1,i2,j1,j2
   integer            :: k1,k2,spos,i
   do i=s1,s2
    call choose_randomly(i1,i2,j1,j2,k1,k2)
    spos = randint(0, 32*np(compound) )
    children(i)%genotype = parents(k1)%genotype(:spos) // parents(k2)%genotype(spos+1:)
   end do
   return
  end subroutine crossover
!
  subroutine choose_randomly(kk1,kk2,jj1,jj2,ii1,ii2)
   implicit none
   integer,intent(in)  :: kk1,kk2,jj1,jj2
   integer,intent(out) :: ii1,ii2
   ii1  = randint(kk1, kk2)
   ii2  = randint(jj1, jj2)
   do while ( ii1 == ii2 )
    ii2 = randint(jj1,jj2)
   end do
   return
  end subroutine choose_randomly
!
  subroutine Mate(compound)
   implicit none
   integer             :: i,i1,i2,spos
   integer, intent(in) :: compound
   real                :: rrr
   ! Retain the first 25% of the childrens
   call Elitism()
   ! Crossover:
   call crossover(ga_Motleists+1,ga_size-100,1,GA_elitists,GA_elitists+1,GA_size)
   call Crossover(GA_elitists+1,ga_Motleists,1,GA_elitists,1,GA_elitists)
   ! Replace the last 100 of the childrens by new childrens
   children(GA_size-100:GA_size) = new_citizen(compound)
   do i = 1, ga_size
    ! Mutation:
    ! I'll respect the first 10% of the childrens of the mutations.
    if(i>=GA_elitists*0.1)then
     rrr = r4_uniform(0.0,1.0)
     if ( rrr < GA_MutationRate) call Mutate(children(i),compound)
    end if
    call UpdateCitizen(children(i),compound)
   end do
   return
  end subroutine Mate
!
  subroutine Fit(compound)
   implicit none
   integer,intent(in)     :: compound
   integer,parameter  :: maxstep = 500, minstep = 10
   integer            :: kk, ii, i, k,vgh
   real               :: diff = 0.0, fit0 = 0.0
   integer            :: eps
   kk = 0
   ii = 0
   pop_alpha = [(new_citizen(compound), i = 1,ga_size)]
   parents =>  pop_alpha
   children => pop_beta
   if ( refit_flag ) then
    do i = 1,GA_elitists 
     parents(i)%genotype=string_IEEE(compound)
     children(i)%genotype=string_IEEE(compound)
     call UpdateCitizen(parents(i),compound)
     call UpdateCitizen(children(i),compound)
    end do
    call Mate(compound)
    call Swap()
    call SortByFitness()
   end if
   call WriteCitizen(1,ii,eps,compound,0, vgh )
   converge: do while ( .true. )
    ii=ii+1
    if ( ii == 1 ) write(6,'(a)') 'Go Down the Rabbit Hole >'
    call SortByFitness()
    call WriteCitizen(1,ii,eps,compound,kk, vgh )
    diff = eps
    eps = Biodiversity( compound, children, vgh )
    fire: if ( FlagFire ) then
     if ( ii >= minstep .and. parents(1)%fitness <= TolFire ) exit converge
    else
     if( ii>=minstep .and. parents(1)%fitness <= 0.1 .and. abs(parents(1)%fitness-fit0) <= 1e-5)then
      kk = kk + 1
     else
      kk = 0
     end if
     if ( ii >= maxstep .or. kk >= 20 ) exit converge
    end if fire
    call Mate(compound)
    call Swap()
    fit0 = parents(1)%fitness
   end do converge
   do i = 0, np( compound )-1
    param( compound,i ) = children(1)%phenotype(i+1)
   end do
   write(111,*)'#',(param(compound,i),i=0,np(compound )-1)
   write(111,*)'#','Fitness:',fit0,'Similarity:',eps
   return
  end subroutine fit
end module mod_genetic
!
module mod_simplex
 use mod_random
 use gaiast_globals
 use qsort_c_module
 use mod_genetic
 implicit none
 private
 public  :: fit_simplex
 contains
! ======================================================================
! Interface between program and module:
subroutine fit_simplex()
 implicit none
 real               :: e = 1.0e-4, scale = 1.0
 integer            :: iprint = 0
 integer            :: i,j
 real               :: sp(0:np(compound)-1)
 type(typ_ga)       :: axolotl
 write(6,'(a)')' '
 write(6,'(a)')'Minimising the cost function using the Nelder-Mead SIMPLEX method:'
 write(6,*)'# Compound:',compound, np(compound)
 do i = 0,np(compound)-1
  axolotl%phenotype(i+1) = param(compound,i)
  sp(i) = axolotl%phenotype(i+1)
 end do
 call simplex(sp,compound,np(compound),e,scale,iprint)
 do i= 0, np(compound)-1
  axolotl%phenotype(i+1)=sp(i)
  param(compound,i)=sp(i)
 end do
 write(111,*)'#',(param(compound,i),i=0,np(compound )-1)
 write(111,*)'#','Fitness:',func(np(compound),sp,compound)
end subroutine fit_simplex

! ======================================================================
! This is the function to be minimized
real function func(n,x,compound) result(rosen)
  implicit none
  integer,intent(in)   :: n
  real,   intent (in)  :: x(0:n-1)
  integer,intent(in)   :: compound
  real                 :: sp(1:np(compound))
  integer              :: i
  sp = 0.0
  do i=0,np(compound)-1
   sp(i+1) = x(i)
  end do
  rosen = fitness(sp,compound)
  return
end function func
! 
subroutine simplex(start, compound, n, EPSILON, scale, iprint)
! This is the simplex routine
! Michael F. Hutt
  implicit none
  integer, intent (in)                   :: n, iprint, compound
  real, intent (inout), dimension(0:n-1) :: start
  real, intent (in)                      :: EPSILON, scale
  integer, parameter                     :: MAX_IT = 1000000
  real, parameter                        :: ALPHA=1.0
  real, parameter                        :: BETA=0.5
  real, parameter                        :: GAMMA=2.0
! ======================================================================
! Variable Definitions
! vs = vertex with the smallest value
! vh = vertex with next smallest value
! vg = vertex with largest value
! i,j,m,row
! k = track the number of function evaluations
! itr = track the number of iterations
! v = holds vertices of simplex
! pn,qn = values used to create initial simplex
! f = value of function at each vertex
! fr = value of function at reflection point
! fe = value of function at expansion point
! fc = value of function at contraction point
! vr = reflection - coordinates
! ve = expansion - coordinates
! vc = contraction - coordinates
! vm = centroid - coordinates
! min
! fsum,favg,s,cent
! vtmp = temporary array passed to FUNC
! ======================================================================
  Integer :: vs,vh,vg
  Integer :: i,j,k,itr,m,row
  real, dimension(:,:), allocatable :: v
  real, dimension(:), allocatable  :: f
  real, dimension(:), allocatable :: vr
  real, dimension(:), allocatable :: ve
  real, dimension(:), allocatable :: vc
  real, dimension(:), allocatable :: vm
  real, dimension(:), allocatable :: vtmp
  real :: pn,qn
  real :: fr,fe,fc
  real :: min,fsum,favg,cent,s

  allocate (v(0:n,0:n-1))
  allocate (f(0:n))
  allocate (vr(0:n-1))
  allocate (ve(0:n-1))
  allocate (vc(0:n-1))
  allocate (vm(0:n-1))
  allocate (vtmp(0:n-1))

! create the initial simplex
! assume one of the vertices is 0.0

  pn = scale*(sqrt(n+1.)-1.+n)/(n*sqrt(2.))
  qn = scale*(sqrt(n+1.)-1.)/(n*sqrt(2.))

  DO i=0,n-1
    v(0,i) = start(i)
  END DO

  DO i=1,n
    DO j=0,n-1
      IF (i-1 == j) THEN
        v(i,j) = pn + start(j)
      ELSE
        v(i,j) = qn + start(j)
      END IF
    END DO
  END DO


! find the initial function values

  DO j=0,n
! put coordinates into single dimension array
! to pass it to FUNC
    DO m=0,n-1
      vtmp(m) = v(j,m)
    END DO
    f(j)=FUNC(n,vtmp,compound)
  END DO

! Print out the initial simplex
! Print out the initial function values
! find the index of the smallest value for printing
  IF (iprint == 0) THEN
   vs=0
   DO j=0,n
    If (f(j) .LT. f(vs)) Then
      vs = j
    END IF
   END DO
! print out the value at each iteration
   Write(6,'(a)') "Initial Values from genetic algorithm:"
   Write(6,*) (v(vs,j),j=0,n-1),'Fit:',f(vs)
  END IF
  k = n+1
! begin main loop of the minimization

DO itr=1,MAX_IT
! find the index of the largest value
  vg = 0
  DO j=0,n
    IF (f(j) .GT. f(vg)) THEN
      vg = j
    END IF
  END DO

! find the index of the smallest value
  vs = 0
  DO j=0,n
    If (f(j) .LT. f(vs)) Then
      vs = j
    END IF
  END DO

! find the index of the second largest value
  vh = vs
  Do j=0,n
    If ((f(j) .GT. f(vh)) .AND. (f(j) .LT. f(vg))) Then
      vh = j
    END IF
  END DO

! calculate the centroid
  DO j=0,n-1
  cent = 0.0
    DO m=0,n
      If (m .NE. vg) Then
        cent = cent + v(m,j)
      END IF
    END DO
    vm(j) = cent/n
  END DO

! reflect vg to new vertex vr
  DO j=0,n-1
    vr(j) = (1+ALPHA)*vm(j) - ALPHA*v(vg,j)
  END DO
  fr = FUNC(n,vr,compound)
  k = k+1

  If ((fr .LE. f(vh)) .AND. (fr .GT. f(vs))) Then
    DO j=0,n-1
      v(vg,j) = vr(j)
    END DO
    f(vg) = fr
  END IF

! investigate a step further in this direction
  If (fr .LE. f(vs)) Then
    DO j=0,n-1
      ve(j) = GAMMA*vr(j) + (1-GAMMA)*vm(j)
    END DO
    fe = FUNC(n,ve,compound)
    k = k+1

! by making fe < fr as opposed to fe < f(vs), Rosenbrocks function
! takes 62 iterations as opposed to 64.

    If (fe .LT. fr) Then
      DO j=0,n-1
        v(vg,j) = ve(j)
      END DO
      f(vg) = fe
    Else
      DO j=0,n-1
        v(vg,j) = vr(j)
      END DO
      f(vg) = fr
    END IF
  END IF

! check to see if a contraction is necessary
  If (fr .GT. f(vh)) Then
    DO j=0,n-1
      vc(j) = BETA*v(vg,j) + (1-BETA)*vm(j)
    END DO
    fc = FUNC(n,vc,compound)
    k = k+1
    If (fc .LT. f(vg)) Then
      DO j=0,n-1
        v(vg,j) = vc(j)
      END DO
    f(vg) = fc

! at this point the contraction is not successful,
! we must halve the distance from vs to all the
! vertices of the simplex and then continue.
! 10/31/97 - modified C program to account for
! all vertices.

  Else
    DO row=0,n
      If (row .NE. vs) Then
        DO j=0,n-1
          v(row,j) = v(vs,j)+(v(row,j)-v(vs,j))/2.0
        END DO
      END IF
    END DO
    DO m=0,n-1
      vtmp(m) = v(vg,m)
    END DO
    f(vg) = FUNC(n,vtmp,compound)
    k = k+1

    DO m=0,n-1
      vtmp(m) = v(vh,m)
    END DO
    f(vh) = FUNC(n,vtmp,compound)
    k = k+1
    END IF
  END IF
! find the index of the smallest value for printing
  !vs=0
  !DO j=0,n
  !  If (f(j) .LT. f(vs)) Then
  !    vs = j
  !  END IF
  !END DO
! print out the value at each iteration
  !IF (iprint == 0) THEN
  !  Write(6,*) "Iteration:",itr,(v(vs,j),j=0,n-1),'Value:',f(vs)
  !END IF
! test for convergence
  fsum = 0.0
  DO j=0,n
    fsum = fsum + f(j)
  END DO
  favg = fsum/(n+1.)
  !s = 0.0
  !DO j=0,n
  !  s = s + ((f(j)-favg)**2.)/n
  !END DO
  !s = sqrt(s)
  If (favg .LT. EPSILON.or.itr==MAX_IT) Then
! print out the value at each iteration
   Write(6,*) "Final Values:", itr
   Write(6,*) (v(vs,j),j=0,n-1),'Fit:',f(vs)
   IF(itr/=MAX_IT)then
    write(6,*)'Nelder-Mead has converged:',favg,'<',epsilon
   else
    write(6,*)'Maximun number of steps:',itr,'=',MAX_IT
   end if
   EXIT ! Nelder Mead has converged - exit main loop
  END IF
END DO
! end main loop of the minimization
! find the index of the smallest value
  vs = 0
  DO j=0,n
    If (f(j) .LT. f(vs)) Then
      vs = j
    END IF
  END DO
!  print out the minimum
  DO m=0,n-1
    vtmp(m) = v(vs,m)
  END DO
  min = FUNC(n,vtmp,compound)
  !write(6,*)'The minimum was found at ',(v(vs,k),k=0,n-1)
  !write(6,*)'The value at the minimum is ',min
  DO i=0,n-1
   start(i)=v(vs,i)
  END DO
250  FORMAT(A29,F7.4)
300  FORMAT(F11.6,F11.6,F11.6)
  return
  end subroutine simplex
! ======================================================================
end module mod_simplex
!
program main
 use mod_random
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
 call bounds()
 if(ncomponents>2.or.Always_Use_Multicomponent_Solver)then
  ! Always_Use_Multicomponent_Solver defined in the header of gaiast_globals module
  call IAST_multicomponent()
 else
  ! It works fine
  call IAST_binary()
 end if
 call cite()
 close(111)
 write(6,'(a)')'====================================='
 stop ':)'
end program main
