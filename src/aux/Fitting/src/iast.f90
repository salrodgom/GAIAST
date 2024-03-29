module gaiast_globals
 use general
 use mod_random
 implicit none
 integer                         :: npar,ncomponents,i,ii,j
 integer(8)                      :: seed = 0
 integer,allocatable             :: np(:)
 integer                         :: err_apertura
 integer                         :: intervalos = 4096
 integer,parameter               :: integration_points = 10000
 real,parameter                  :: precision_Newton = 1e-5
 real                            :: tol = 0.001, tolfire = 0.25
 real,allocatable                :: param(:,:),concy(:),pi(:,:),iso(:,:)
 integer,allocatable             :: npress(:)
 integer                         :: compound
 integer,parameter               :: maxnp = 10, maxcompounds = 10, maxdata = 1000
 real,target                     :: datas(2,maxcompounds,maxdata),f(maxdata)
 real,pointer                    :: x(:),y(:),alldat(:,:)
 character(100),allocatable      :: ajuste(:)
 character(32*maxnp),allocatable :: string_IEEE(:)
 character(100)                  :: line,string,intmethod,Equation
 character(5)                    :: inpt
 logical                         :: flag = .true., FlagFire = .false.,seed_flag=.true.
 logical                         :: physical_constrains = .false., range_flag =.true.
 logical                         :: Always_Use_Multicomponent_Solver=.false., IAST_flag = .true.
 logical                         :: refit_flag = .false.
 real,parameter                  :: R = 0.008314472 ! kJ / mol / K
 real                            :: T = 298.0
 real                            :: inferior
 real                            :: superior
 contains
!
  subroutine read_input()
  implicit none
  integer                    :: ii,i,j
  character(len=4)           :: realorbinary
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
    read(line,*)inpt,flag,realorbinary
    if(flag.eqv..false.)then
     write(6,*)'[WARN] Fit is already done'
     select case(realorbinary)
      case('real')
       readparameter_real: do ii=1,ncomponents
        read(5,*)(param(ii,j),j=0,np(ii)-1)
       end do readparameter_real
      case default
       readparameter_bin: do ii=1,ncomponents
        do j=1,np(ii)
         read(5,'(a32)') string_IEEE(ii)(32*(j-1)+1:32*j)
         param( ii,j-1 ) = bin2real( string_IEEE(ii)(32*(j-1)+1:32*j) )
        end do
       end do readparameter_bin
      end select
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
     write(6,'(a)') string_IEEE(1:ncomponents)(1:32*np(ii))
    else
     string_IEEE(1:ncomponents)=' '
    end if
   end if
   if(line(1:5)=='inter') then
    read(line,*)inpt,intervalos
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
   if(line(1:5)=="IAST ") read(line,*)inpt, IAST_flag
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
  ! Allocate defaults options:
  allocate(pi(ncomponents,0:intervalos),iso(ncomponents,0:intervalos))
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
    presion(1,i) = exp(inferior+(i+1)*dx)
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
  integer                  :: i
  integer,parameter        :: imax = 100000000
  real                     :: x, delta, integral
  real,intent(in)          :: a(0:n-1), x0,oldPi,piValue
  character(100),intent(in):: funk
  character(100)           :: mode = 'version1'
  integral = 0.0
  x = 0.0
  i = 0
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
    i = i+1
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
 subroutine clear_screen
  write(*,'(a1,a)',advance="no") achar(27), '[2J'
 end subroutine clear_screen
end module gaiast_globals
!
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
