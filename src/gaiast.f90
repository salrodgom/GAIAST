! GAIAST
! Salvador Rodriguez-Gomez
! Rocio Bueno-Perez
! Sofia Calero-Diaz
! Wednesday, 09, November, 2016
module mod_random
! module for pseudo random numbers
 implicit none
 private
 public init_random_seed, randint, r4_uniform
 contains
 subroutine init_random_seed(seed)
  implicit none
  integer, intent(out) :: seed
! local
  integer   day,hour,i4_huge,milli,minute,month,second,year
  parameter (i4_huge=2147483647)
  double precision temp
  character*(10) time
  character*(8) date
  call date_and_time (date,time)
  read (date,'(i4,i2,i2)')year,month,day
  read (time,'(i2,i2,i2,1x,i3)')hour,minute,second,milli
  temp=0.0D+00
  temp=temp+dble(month-1)/11.0D+00
  temp=temp+dble(day-1)/30.0D+00
  temp=temp+dble(hour)/23.0D+00
  temp=temp+dble(minute)/59.0D+00
  temp=temp+dble(second)/59.0D+00
  temp=temp+dble(milli)/999.0D+00
  temp=temp/6.0D+00
  doext: do
    if(temp<=0.0D+00 )then
       temp=temp+1.0D+00
       cycle doext
    else
       exit doext
    end if
  enddo doext
  doext2: do
    if (1.0D+00<temp) then
       temp=temp-1.0D+00
       cycle doext2
    else
       exit doext2
    end if
  end do doext2
  seed=int(dble(i4_huge)*temp)
  if(seed == 0)       seed = 1
  if(seed == i4_huge) seed = seed-1
  return
 end subroutine init_random_seed
!
 integer function randint(i,j,seed)
  integer,intent(in) :: i,j,seed
  real               :: r
  CALL RANDOM_NUMBER(r)
  randint=int(r*(j+1-i))+i
 end function randint
! 
 real function r4_uniform(x,y,seed)
  implicit none
  real,intent(in)    :: x,y 
  real               :: r
  integer,intent(in) :: seed
  CALL RANDOM_NUMBER(r)
  r4_uniform=(r*(y-x))+x
  return
 end function r4_uniform
end module

module newton
    implicit none
    integer, parameter      :: maxiter = 20
    real(kind=8), parameter :: tol = 1.d-14
contains
subroutine solve(f, fp, x0, x, iters, debug)
    ! Estimate the zero of f(x) using Newton's method.
    ! Input:
    !   f:  the function to find a root of
    !   fp: function returning the derivative f'
    !   x0: the initial guess
    !   debug: logical, prints iterations if debug=.true.
    ! Returns:
    !   the estimate x satisfying f(x)=0 (assumes Newton converged!)
    !   the number of iterations iters
    implicit none
    real(kind=8), intent(in)    :: x0
    real(kind=8), external      :: f, fp
    logical, intent(in)         :: debug
    real(kind=8), intent(out)   :: x
    integer, intent(out)        :: iters
    ! Declare any local variables:
    real(kind=8)                :: deltax, fx, fxprime
    integer                     :: k
    ! initial guess
    x = x0
    if (debug) then
        print 11, x
 11     format('Initial guess: x = ', e22.15)
        endif
    ! Newton iteration to find a zero of f(x)
    do k=1,maxiter
        ! evaluate function and its derivative:
        fx = f(x)
        fxprime = fp(x)
        if (abs(fx) < tol) then
            exit  ! jump out of do loop
            endif
        ! compute Newton increment x:
        deltax = fx/fxprime
        ! update x:
        x = x - deltax
        if (debug) then
            print 12, k,x
 12         format('After', i3, ' iterations, x = ', e22.15)
            endif
        enddo
    if (k > maxiter) then
        ! might not have converged
        fx = f(x)
        if (abs(fx) > tol) then
            print *, '[Warning]: calculation has not converged yet'
            endif
        endif
    ! number of iterations taken:
    iters = k-1
end subroutine solve
end module newton

module qsort_c_module
! Recursive Fortran 95 quicksort routine sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms, 1997 printing
! Made F conformant by Walt Brainerd
 implicit none
 public  :: QsortC
 private :: Partition
 contains
 recursive subroutine QsortC(A)
  real(16), intent(in out), dimension(:) :: A
  integer                            :: iq
  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
 end subroutine QsortC
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
 integer                    :: npar,ncomponents,i,seed = 0
 integer,allocatable        :: np(:)
 integer                    :: err_apertura,ii,intervalos,j
 integer,parameter          :: integration_points = 10000
 real,parameter             :: precision_Newton = 1e-5
 real                       :: tol = 0.001, tolfire = 0.25
 real,allocatable           :: param(:,:),concy(:),pi(:,:),iso(:,:)
 integer,allocatable        :: npress(:)
 integer,parameter          :: maxnp = 10, maxcompounds = 10, maxdata = 1000
 real,target                :: datas(2,maxcompounds,maxdata),f(maxdata)
 real,pointer               :: x(:),y(:),alldat(:,:)
 character(100),allocatable      :: ajuste(:)
 character(32*maxnp),allocatable :: string_IEEE(:)
 character(100)             :: line,string,intmethod
 character(5)               :: inpt
 logical                    :: flag = .true., FlagFire = .false.,seed_flag=.true.
 logical                    :: physical_constrains = .false., range_flag =.true.
 logical                    :: refit_flag = .false.
 real,parameter             :: R = 0.008314472 ! kJ / mol / K
 real                       :: T = 298.0
 real                       :: inferior
 real                       :: superior
 contains
  subroutine read_input()
  implicit none
  read_input_do: do
   read(5,'(A)',iostat=err_apertura)line
   if ( err_apertura /= 0 ) exit read_input_do
   if(line(1:1)=='#') cycle read_input_do
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
    cycle read_input_do
   end if
   if(line(1:5)=='refit') then
    write(6,*)'[WARN] Refitting parameters'
    read(line,*)inpt,refit_flag
    if(refit_flag)then
    do ii=1,ncomponents
     do j=1,np(ii)
      read(5,'(a32)')  string
      read(string,'(a32)') string_IEEE(ii)(32*(j-1)+1:32*j)
     end do
     write(6,'(a)') string_IEEE(ii)
    end do
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
     read(line,*)inpt, FlagFire
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
  integer            :: ib,jb,compound
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
!   ...
    y0 = exp( inferior + (ib+1)*dx )
    x0 = model(apar, np(compound), y0, funk )
!   ...
    x0 = (yyy(compound) + x0)*dx/2.0
    s1(compound) = s1(compound) + x0
    pi(compound,ib+1) = s1(compound)
    iso(compound,ib+1) = yyy(compound)
    deallocate(apar)
   end do
   write(222,*)(xxx(compound) ,yyy(compound), pi(compound,ib+1), compound=1,ncomponents)
  end do
  close(222)
  return
 end subroutine bounds

 subroutine IAST_binary()
  implicit none
  integer        :: i,j,k,ijk
  real           :: comp1,comp2,dx,p,concx(ncomponents),aux1,aux2
  real           :: presion(ncomponents),n(0:ncomponents)
  dx = real((superior-inferior)/intervalos)
  open(unit=104,file='adsorcion.dat')
  do i=0,intervalos-1
   do j=0,intervalos-1
    comp1 = abs(pi(1,i)-pi(2,j))
    if (comp1 <= tol) then
     presion(1) = exp(inferior+i*dx)
     presion(2) = exp(inferior+j*dx)
! {{ Calculate molar fraction of component 1 at spread pressure Pi from a general formula: x1=1/(1+P1Y2/P2Y1 + w(o3) )
     concx(1) = presion(2)*concy(1) / (concy(2)*presion(1) + concy(1)*presion(2))
     concx(2) = 1.0 - concx(1)   ! }}
     p  = presion(1)*concx(1)/concy(1)
     n(1) = iso(1,i)*concx(1)
     n(2) = iso(2,j)*concx(2)
     n(0) = n(1) + n(2)
     write(104,*)p,(n(ijk),ijk=1,ncomponents),n(0)
    end if
   end do
  end do
  close(104)
  return
 end subroutine IAST_binary

 subroutine IAST_multicomponent()
  implicit none
  integer          :: i,j,ijk,nn
  real,allocatable :: a(:)
  character(100)   :: funk
  real           :: comp,dx,p,concx(ncomponents),pureloading(ncomponents)
  real           :: presion(0:ncomponents,0:intervalos-1),n(0:ncomponents)
  real           :: auxP,aux1,aux2,lastPi(ncomponents),piValue=0.0,x,last=0.0
  logical        :: validPressure = .true., yIsZero = .true.,flag = .true.
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
    !presion(1,i) = exp(log(datas(1,1,1))+(log(auxP)-log(datas(1,1,1)))*(i+1)/intervalos)
    last    = lastPi(1)                       !<---- works
    piValue = CalculatePi(1,presion,i,last)   !<---- works
    validPressure=CalculatePressures(presion,i,lastPi,piValue)   !<- WRONG!
    !write(6,*)i,(presion(ijk,i),ijk=1,ncomponents),piValue,validPressure
    if(validPressure.and.i>0)then
     !write(6,*)i,(presion(ijk,i),ijk=1,ncomponents),piValue
     continue
    else
     i=i+1
     cycle ScanPressures
    end if
    !lastPi = piValue
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

 logical function CalculatePressures(presion,point,lastPi,piValue)
  implicit none
  integer,intent(in) :: point
  real,intent(inout) :: presion(0:ncomponents,0:intervalos-1),piValue
  real,intent(inout) :: lastPi(ncomponents)
  real               :: oldPi, p,a1,a2,b1,b2,c,a,p1,p2
  integer            :: k,j,n
  real,allocatable   :: apar(:)
  character(100)     :: funk,mode= 'integral' ! 'Newton'
  CalculatePressures = .true.
  k=1
  calc: do while (k<=ncomponents.and.CalculatePressures)
! {{ interpolation and parameters from k-compound
   funk = ajuste( k )
   !if(funk=='langmuir') mode = 'analytical' !.or.funk=='langmuir_dualsite') mode = 'analytical'
   n = np(k)
   allocate(apar(0:n-1))
   do j = 0, n-1
    apar(j) = param(k,j)
   end do
! }}
   select case (mode)
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

 real function integrate(x0,x1,a,n,funk)
  implicit none
  integer              ::  i
  integer,intent(in)   ::  n
  real                 ::  delta,x
  real,intent(in)      ::  x0,x1
  real                 ::  factor
  real,intent(in)      ::  a(0:n-1)
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
 end select
 return
 end subroutine MakeInitPOP

 real function model(a,n,xx,funk)
  implicit none
  integer,intent(in)        :: n
  real,intent(in)           :: a(0:n-1)
  real,intent(in)           :: xx
  character(100),intent(in) :: funk
  select case (funk)
   case("freundlich","f")
    model = a(0)*xx**a(1)
   case("langmuir","l")
    model = a(0)*a(1)*xx/(1+a(1)*xx)
   case("langmuir_freundlich","lf")
    model = a(0)*a(1)*xx**a(2)/(1.0+a(1)*xx**a(2))
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
    model = (a(0)*a(1)*xx)/(1+a(1)*xx) + (a(2)*a(3)*xx**a(4))/(1+a(3)*xx*a(4))
   case ("jensen_seaton","jsn")
    model = a(0)*xx*( 1.0 + ( a(0)*xx / (a(1)*( 1+a(2)*xx )) )**a(3) )**(-1.0/a(3)) 
  end select
  return
 end function model
end module gaiast_globals

module mod_genetic
 use mod_random
 use gaiast_globals
 use qsort_c_module
 implicit none
 private
 public fit
 integer,parameter             :: ga_size     = 2**13 ! numero de cromosomas
 real,parameter                :: ga_mutationrate = 0.3333 !2000/real(ga_size) ! ga_mutationrate=0.333
 real,parameter                :: ga_eliterate= 0.25, GA_DisasterRate = 0.0000001
 integer,parameter             :: maxlinelength=maxnp*32
 integer,parameter             :: ga_elitists = int( ga_size * ga_eliterate)
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

  type(typ_ga) function new_citizen(compound,seed)
   implicit none
   integer :: i,compound,seed,j,k
   new_citizen%genotype = ' '
   do i = 1,32*np(compound)
    !j=randint(48,49,seed)
    !write(6,*)j,achar(j)
    new_citizen%genotype(i:i) = achar(randint(48,49,seed))
   end do
   !stop
   !write(6,'(a)') new_citizen%genotype
   do i = 1,np(compound)
    read(new_citizen%genotype(32*(i-1)+1:32*i),'(b32.32)') new_citizen%phenotype(i)
   end do
   new_citizen%fitness = fitness( new_citizen%phenotype,compound)
   return
  end function new_citizen

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

!  real function Fitness(phenotype,compound)
!   implicit none
!   real, intent(in)    :: phenotype(maxnp)
!   integer             :: i,compound,k = 0
!   real                :: a(0:np(compound)-1),xx,yy
!   character(len=100)  :: funk
!   logical             :: flagzero = .false.
!   real                :: infinite = HUGE(40843542)
!!   funk = ajuste(compound)
!   do i = 0,np(compound)-1
!    a(i) = phenotype(i+1)
!   end do
!   fitness = 0.0
!   do i = 1, npress(compound)
!    xx = datas(1,compound,i)
!    yy = datas(2,compound,i)
!    fitness = fitness + 0.5*( yy - model(a,np(compound),xx,funk) )**2
!   end do
!   return
!  end function Fitness
  real function Fitness(phenotype,compound)
   implicit none
   real, intent(in)    :: phenotype(maxnp)
   integer             :: i,compound,k = 0
   real                :: a(0:np(compound)-1),xx,yy
   character(len=100)  :: funk
   logical             :: flagzero = .false.
   !logical             :: physical_constrains = .false.
   real                :: infinite = 50.0,penalty = 0.0
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
      case default
       penalty = 0.0
     end select
   end if phys_constrains
   fitness = penalty
   do i = 1, npress(compound)
    xx = datas(1,compound,i)
    yy = datas(2,compound,i)
    fitness = fitness + 0.5*( yy - model(a,np(compound),xx,funk) )**2
   end do
   if(isnan(fitness)) fitness      = 99999999.999999
   if(fitness==0.0000000000)fitness= 99999999.999999
   return
  end function Fitness

  subroutine WriteCitizen(k,kk,kkk,compound)
   implicit none
   integer                        :: k,kk,i,compound
   character(len=100)             :: fmt_
   character(len=32*np(compound)) :: wnowaste
   real                           :: wnowasteparam(1:32*np(compound)),wfitness,kkk
   do i=1,32*np(compound)
    wnowaste(i:i)=' '
   end do
   ! ...
   wnowaste = parents(k)%genotype
   do i=1,np(compound)
    wnowasteparam(i) = parents(k)%phenotype(i)
   end do
   wfitness = parents(k)%fitness
   write(6,'(i2,1x,i5,1x,a32,1x,e25.12,1x,f20.10,1x,a,1x,f20.10,1x,a)')compound,kk,wnowaste(1:32),&
        wnowasteparam(1),wfitness,'[Fitness]',kkk,'[Similarity]' !,k
   do i=2,np(compound)
    write(6,'(9x,a32,1x,e25.12)')wnowaste(1+32*(i-1):32*i),wnowasteparam(i)
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

  real function Biodiversity( compound )
   implicit none
   integer,intent(in)             :: Compound
   character(len=32*np(compound)) :: str1,str2
   integer                        :: k
   character(len=20)              :: mode = 'Superficial'
   logical                        :: flag = .true.
   real                           :: error_d = 1e-3
   select case (mode)
    case('Superficial')
     do k = 1, ga_size
      do j = k+1,ga_size
      !if(children(k)%genotype == children(j)%genotype) then
      ! Biodiversity = Biodiversity + 1.0 /real(ga_size * ga_size)
      !end if
       dbio: do i = 1,np(compound)
        bio: if( abs(children(k)%phenotype(i) - children(j)%phenotype(i)) <= error_d )then
         flag = .true.
         exit dbio
        end if bio
       end do dbio
       if( flag ) then
        Biodiversity = Biodiversity + 1.0 /real(ga_size * ga_size)
       end if
      end do
     end do
     Biodiversity = 1.0/(Biodiversity)
    case('Deep')
     do k = 1, ga_size
      do j = k+1,ga_size
        Biodiversity = Biodiversity + &
         sum([(abs(iachar( children(k)%genotype(i:i)) - &
          iachar( children(j)%genotype(i:i))/real(32*np(compound))),i=1,32*np(compound))])
      end do
     end do
     Biodiversity = 1.0/(Biodiversity*ga_size*ga_size)
   end select
   return
  end function Biodiversity

  subroutine Mutate( macrophage , compound )
   implicit none
   type(typ_ga), intent(inout) :: macrophage
   integer                     :: ipos,compound
   do i = 1,np(compound)
    ipos = randint(32*(i-1)+1,32*i,seed)
    macrophage%genotype(ipos:ipos) = achar(randint(48,49,seed))
   end do
   return
  end subroutine Mutate

  subroutine NuclearDisaster(Compound)
   implicit none
   integer,intent(in) ::  Compound
   integer            :: k = 0, i, j
   real               :: rrr
   do i = GA_ELITISTS + 1, GA_Size
    do j=1,32*np(compound)
     Children%genotype(j:j) = achar(randint(48,49,seed))
    end do
    call UpdateCitizen(Children(i),Compound)
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

  subroutine Mate(compound)
   integer             :: i, i1, i2, spos
   integer, intent(in) :: compound
   real                :: rrr
   call Elitism()
   do i = GA_ELITISTS + 1, ga_size
    ! Crossover:
    ! {{ eleccion random del primer 50% de la tabla
    call choose_randomly(i1,i2)
    ! }}
    ! {{ eleccion proporcionalmente a su fitness
    !call choose_propto_fitness(i1,i2)
    !write(6,*)i1,i2
    ! }}
    spos = randint(0, 32*np(compound), seed )
    children(i)%genotype = parents(i1)%genotype(:spos) // parents(i2)%genotype(spos+1:)
    ! Mutate and NuclearDisaster:
    rrr = r4_uniform(0.0,1.0,seed)
    if ( rrr < GA_MUTATIONRATE) then
     call Mutate(children(i),compound)
    else if ( rrr >= GA_MutationRate .and. rrr <= GA_MutationRate + GA_DisasterRate ) then
     call NuclearDisaster(Compound)
    end if
   end do
   do i = 1, ga_size
    call UpdateCitizen(children(i),compound)
   end do
   return
  end subroutine Mate

  subroutine choose_randomly(j1,j2)
   implicit none
   integer,intent(out) :: j1,j2
   j1  = randint(1, int(ga_size/2),seed)
   j2  = randint(1, int(ga_size/2),seed)
   do while ( j1 == j2 )
    j2 = randint(1, int(ga_size/2),seed)
   end do
   return
  end subroutine choose_randomly

  subroutine choose_propto_fitness(j1,j2)
   implicit none
   integer,intent(out) :: j1,j2
   integer             :: i
   real                :: ftnss(ga_size),prop(0:ga_size)=0.0,rrr1,rrr2
   real                :: infinity = HUGE(2147483647)
   rrr1 = 0.0
   do i = 1, ga_size
    ftnss(i) = 1.0/parents(i)%fitness
    if ( isnan( parents(i)%fitness ) ) ftnss(i) = 0.0
    if ( parents(i)%fitness > infinity ) ftnss(i) = 0.0
    prop(i) = ftnss(i)
    if( ftnss(i) >= infinity ) then
      rrr1 = rrr1 + infinity
    else
      rrr1 = rrr1 + ftnss(i)
    end if
   end do
   prop = prop / rrr1
   ! select 1:
    rrr1 = r4_uniform(0.0,1.0,seed)
    slct1: do i=1,ga_size
     if(rrr1<=prop(i-1).and.rrr1>prop(i))then
      j1 = i
     end if
    end do slct1
    ! select 2:
    rrr2 = r4_uniform(0.0,1.0,seed)
    do while ( rrr1 == rrr2 )
     rrr2 = r4_uniform(0.0,1.0,seed)
    end do
    slct2: do i=1,ga_size
     if(rrr2<=prop(i-1).and.rrr2>prop(i))then
      j2 = i
     end if
    end do slct2
   return
  end subroutine choose_propto_fitness

  subroutine Fit(Compound,Seed)
   implicit none
   integer,intent(in) :: Compound, Seed
   integer,parameter  :: maxstep = 100, minstep = 10
   integer            :: kk, ii, i, k
   real               :: eps = 0.0, diff = 0.0, fit0 = 0.0
   kk = 0
   ii = 0
   pop_alpha = [(new_citizen(compound,seed), i = 1,ga_size)]
   parents =>  pop_alpha
   children => pop_beta
   if(refit_flag)then
    do i = 1, 2
     parents(i)%genotype=string_IEEE(compound)
     children(i)%genotype=string_IEEE(compound)
     call UpdateCitizen(parents(i),compound)
     call UpdateCitizen(children(i),compound)
    end do
    call Mate(compound)
    call Swap()
    call SortByFitness()
   end if
   call WriteCitizen(1,ii,eps,compound)
   converge: do while (.true.)
    ii=ii+1
    call SortByFitness()
    !do i=1,ga_size
     call WriteCitizen(1,ii,eps,compound)
    !end do
    !STOP
    diff = eps
    eps = Biodiversity( compound )
    fire: if ( FlagFire ) then
     if ( ii >= minstep .and. parents(1)%fitness <= TolFire ) exit converge
    else
     if( abs(diff - eps) <= 0.1 .and. ii >= minstep .and. &
      parents(1)%fitness - fit0 == 0 ) then
      kk = kk + 1
     else
      kk = 0
     end if
     if ( ii >= maxstep .or. kk >= 10 ) exit converge
    end if fire
    call Mate(compound)
    call Swap()
    !call SortByFitness()
    fit0 = parents(1)%fitness
   end do converge
   do i = 0, np( compound )-1
    param( compound,i ) = children(1)%phenotype(i+1)
   end do
   write(111,*)'#',(param(compound,i),i=0,np(compound )-1)
   write(111,*)'#','Fitness:',fit0,'Similarity:',eps,'Rseed',seed
   return
  end subroutine fit
end module mod_genetic

program main
 use mod_random
 use gaiast_globals
 use mod_genetic
 !use iso_fortran_env
 !print '(4a)', 'This file was compiled by ', &
 !      compiler_version(), ' using the options ', &
 !      compiler_options()
 print '(a)',' ',&
"#     ________    _____   .___    _____     ____________________  ",&
"#    /  _____/   /  _  \  |   |  /  _  \   /   _____/\__    ___/  ",&
"#   /   \  ___  /  /_\  \ |   | /  /_\  \  \_____  \   |    |     ",&
"#   \    \_\  \/    |    \|   |/    |    \ /        \  |    |     ",&
"#    \______  /\____|__  /|___|\____|__  //_______  /  |____|     ",&
"#           \/         \/              \/         \/              ",&
" "
 call read_input()
 if (seed_flag) then
  seed = 73709
  call init_random_seed(seed)
  seed_flag=.false.
 end if
 write(6,'("Random Seed:",1x,i10)') seed
 call ReadIsotherms()
 if(flag)then
  open(111,file="iso.dat")
  do ii = 1, ncomponents
   write(6,'("Fitting compound:",1x,i2)') ii
   write(6,'(a14,1x,a24,1x,a24,10x,a14)')'Compound/Step:','Cromosome:','Parameters','Control'
   call fit(ii,seed)
  end do
 end if
 call bounds()
 if(ncomponents>2)then
   call IAST_multicomponent()
 else
   call IAST_binary()
 end if
 call cite()
 close(111)
 write(6,'(a)')'====================================='
 stop ':)'
end program main
