module mod_random
! module for pseudo random numbers
 implicit none
 private
 public init_random_seed, randreal, randint
contains
 subroutine init_random_seed()
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t
  call random_seed(size=n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
  form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
   read(un) seed
  close(un)
  else
  ! Fallback to XOR:ing the current time and pid. The PID is useful in case one launches multiple instances of the same program in parallel.
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
  end if
  call random_seed(put=seed)
 contains
 ! This simple PRNG might not be good enough for real work, but is sufficient for seeding a better PRNG.
  function lcg(s)
   integer :: lcg
   integer(int64) :: s
   if (s == 0) then
    s = 104729
   else
   s = mod(s, 4294967296_int64)
   end if
   s = mod(s * 279470273_int64, 4294967291_int64)
   lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function
 end subroutine
    ! Random real (0.0 <= r < 1.0)
 real function randreal()
  call random_number(randreal)
 end function
 ! Random int (a <= r <= b)
 integer function randint(a, b)
  integer, intent(in) :: a, b
  if (a > b) stop "a must be less than or equal to b"
   randint = a + int(randreal() * (b - a + 1))
  end function
end module

program iast_desarrollo_GA
 ! IAST program
 use mod_random
 implicit none
 real             :: comp,concx1,concx2,concy1,concy2
 real,parameter   :: tol = 0.001
 real             :: h1,h2,p,presion1,presion2,n,n1,n2
 real             :: inferior,superior,s1,s2,x0
 integer          :: err_apertura = 0,i,ii,k,l,j,intervalos,npar
 character(100)   :: line,ajuste,test_name
 real,allocatable :: x1(:),y1(:),x2(:),y2(:),coef1(:),coef2(:),PI1(:),PI2(:),area1(:),area2(:)
 real,allocatable :: e_coef1(:),e_coef2(:),param(:,:),eparam(:,:)
 real,allocatable :: puntos1(:),puntos2(:),funcion1(:),funcion2(:)
 
 call init_random_seed()

!========================================================================
! PARAMETROS
!========================================================================
 read(5,*)ajuste
 read(5,*)intervalos
 !read(5,*)tol
 read(5,*)concy1
 read(5,*)concy2
! ajuste = "toth"    !seleccion del ajuste 
! intervalos = 10000 !numero de puntos para calcular el area
! tol = 0.001        !tolerancia para igualar areas
! concy1 = 0.5       !concentracion fase gas componente1
! concy2 = 0.5       !concentracion fase gas componente2
!========================================================================
 nisothermpure: do ii=1,2
  write (test_name, '( "isoterma", I1, ".dat" )' ) ii
  call system ("cp " // test_name // " isotermaN.dat")
  open(unit=100,file="isotermaN.dat",status='old',iostat=err_apertura)
  if(err_apertura/=0) stop '[error] isotermaN.dat'
  k=0 
  do0: do
   READ (100,'(A)',IOSTAT=err_apertura) line 
   IF( err_apertura /= 0 ) EXIT do0
   k=k+1
  END DO do0
  if(ii==1) allocate(x1(1:k),y1(1:k))
  if(ii>=2) allocate(x2(1:k),y2(1:k))
  rewind(100)
  do1: do i=1,k
   read (100,'(A)',IOSTAT=err_apertura) line
   if(err_apertura/=0) EXIT do1
   if(ii==1)read(line,*) x1(i),y1(i)
   if(ii>=2)read(line,*) x2(i),y2(i)
  end do do1
  close(100)
 end do nisothermpure
!========================================================================
 select case (ajuste)
  case("freundlich")
   npar= 2
   ii  = 2**14 
  case ("langmuir")
   npar=2
   ii = 2**15
  case ("toth")
   npar=3
   ii = 2**18
  case ("jensen")
   npar=4
   ii = 2**20
  case ("dubinin_raduschkevich")
   npar=3
   ii = 2*20
  case ("langmuir_dualsite")
   npar=4
   ii = 2**18
  case ("dubinin_astakhov")
   npar=4
   ii = 2**20
 end select
 allocate(param(2,0:npar-1),eparam(2,0:npar-1))
 allocate(coef1(0:npar-1),coef2(0:npar-1),e_coef1(0:npar-1),e_coef2(0:npar-1))
 do i=1,2
  coef1=0
  coef2=0
  write (test_name, '( "isoterma", I1, ".dat" )' ) i
  write (6, '( "Fitting isoterma", I1, ".dat" )' ) i
  call system ("cp " // test_name // " isotermaN.dat")
  call fitgen(coef1,e_coef1,npar,ajuste,ii)
  do j=0,npar-1
   param(i,j)=coef1(j)
   eparam(i,j)=e_coef1(j)
  end do
 end do
 do j=0,npar-1
  coef1(j)=param(1,j)
  coef2(j)=param(2,j)
  e_coef1(j)=eparam(1,j)
  e_coef2(j)=eparam(2,j)
 end do
 deallocate(param,eparam)
 
!========================================================================
! CALCULO DE  AREAS
!========================================================================
 l=size(x2)
 k=size(x1)
 ALLOCATE (puntos1(0:intervalos),puntos2(0:intervalos))
   IF (x1(1)<x2(1)) THEN
    inferior =log10(x1(1))
   ELSE
    inferior =log10(x2(1))
   END IF
   IF (x1(k)<x2(l)) THEN
    superior =log10(x2(l))
   ELSE
    superior =log10(x1(k))
   END IF
   h1 = real((superior-inferior)/intervalos)
   DO i=0,intervalos
    puntos1(i)=10**(inferior+i*h1) 
   END DO
  IF (x1(1)<x2(1)) THEN
        inferior =log10(x1(1))
  ELSE
        inferior =log10(x2(1))
  END IF
  IF (x1(k)<x2(l)) THEN
        superior =log10(x2(l))
  ELSE
        superior =log10(x1(k))
  END IF
  h2 = real((superior-inferior)/intervalos)
  DO i=0,intervalos
     puntos2(i)=10**(inferior+i*h2)
  END DO
 ALLOCATE(funcion1(0:intervalos),funcion2(0:intervalos))
 do i=0,intervalos
  x0 = puntos1(i)
  funcion1(i) = model(coef1,x0,npar,ajuste)
  funcion2(i) = model(coef2,x0,npar,ajuste)
 end do
! calculamos el area bajo la curva de ajuste utilizando el metodo de los trapecios
 ALLOCATE (area1(0:intervalos),area2(0:intervalos),PI1(0:intervalos),PI2(0:intervalos))
 DO i=0,intervalos-1
  area1(i)= real((funcion1(i)+funcion1(i+1))*h1/2)!/puntos1(i)
  area2(i)= real((funcion2(i)+funcion2(i+1))*h2/2)!/puntos2(i)
  s1 = s1 + area1(i)
  s2 = s2 + area2(i)
  PI1(i+1)= s1
  PI2(i+1)= s2
  !
 END DO
 
 OPEN(221,file='iso1.dat')
 OPEN(222,file='iso2.dat')
 do i=1,intervalos
  WRITE(221,*)puntos1(i),funcion1(i),pi1(i)
  WRITE(222,*)puntos2(i),funcion2(i),pi2(i)
 end do
 close(221)
 close(222)
 
!========================================================================
! IAST (binary mixture)
!========================================================================
 open(unit=104,file='adsorcion.dat')
 do i=1,intervalos
  do j=1,intervalos
    comp = abs(PI1(i)-PI2(j))           ! las presiones se igualan
    if (comp <= tol) then              ! ...
     if (abs(x1(1)-x1(k))>100) then
      presion1 = 10**(inferior+i*h1)
     else
      presion1 = inferior+i+h1
     end if
     if (abs(x2(1)-x2(l))>100) then
      presion2 = 10**(inferior+j*h2)
     else
      presion2 = inferior+j*h2
     end if
     ! ... concentraciones ( 2 componentes )
      concx1 = presion2*concy1 / (concy2*presion1 + concy1*presion2)
      concx2 = 1.0 - concx1
     ! ... presion total y loading
      p = presion1*concx1/concy1
      n1 = funcion1(i)*concx1
      n2 = funcion2(j)*concx2
      n = n1+n2
      write(104,*)p,n1,n2,n
     end if
  end do
 end do
 close(104)
 call system ("awk '{print $1,$2,$3,$4}' adsorcion.dat | sort -gk1 > c")
 call system ("mv c adsorcion.dat")
 ! deallocate( todo )
 stop
 CONTAINS
!
! real function IntegrateIsotherm(x0,x1,x,integration_points)
!  implicit none
!  integer              ::  i,integration_points
!  real                 ::  delta
!  real                 ::  x0,x
!  real                 ::  factor
!  delta=(x1-x0)/(integration_points)
!  do i = 0,integration_points
!   factor = 1.0
!   if (i==0.or.i==integration_points-1) factor = 3.0/8.0
!   if (i==1.or.i==integration_points-2) factor = 7.0/6.0
!   if (i==2.or.i==integration_points-3) factor = 23.0/24.0
!   x = x0 + delta*(0.5+i)
!   IntegrateIsotherm=IntegrateIsotherm+factor*delta*/x
!   
!  end do
!
 subroutine fitgen(a,ea,n,funk,GA_POPSIZE)
! crea un vector de valores posibles de ajuste y va mezclando 'genes' para alcanzar
! un ajuste optimo minimizando el coste a una lectura de la isoterma
  implicit none
  integer              ::  i,j,k,l,h,err_apertura
  integer,intent(in)   ::  GA_POPSIZE
  integer,intent(in)   ::  n
  real,intent(out)     ::  a(0:n-1),ea(0:n-1)
  real                 ::  setparam(GA_POPSIZE,0:n-1),fit(GA_POPSIZE)
  character(100)       ::  line
  character(100),intent(in)::funk
  real,allocatable     ::  x(:),y(:)
  real,parameter       ::  mutationP0 = 0.25
  real                 ::  mutationP ! 10%
  real                 ::  renorm = 0.0, suma = 0.0,rrr
  real,parameter       ::  tol = 1.0
  real,parameter       ::  max_range = 100.0
  integer              ::  dimen
  fit = 99999999.99
  open(unit=123,file='isotermaN.dat',iostat=err_apertura)
  if(err_apertura/=0)stop '[error] open file'
  l=0
  dofit0: do
   read(123,'(A)',iostat=err_apertura) line
   if(err_apertura/=0) exit dofit0
   l=l+1
  end do dofit0
  allocate(x(l),y(l))
  rewind(123)
  dofit1: do i=1,l
   read(123,'(A)',iostat=err_apertura) line
   if(err_apertura/=0) exit dofit1
   read(line,*) x(i),y(i)
   write(6,*) x(i),y(i)
  end do dofit1
  close(123)
  dimen=l
  h = 0
  a  = 0.0
  ea = 0.0
  init_set: do k=1,GA_POPSIZE
   do j=0,n-1
    !setparam(k,j) = max_range*randreal()
    if(j==0) setparam(k,j) = max_range*randreal()
    if(j>=1) setparam(k,j) = randreal()
    if(j==3) setparam(k,j) = max_range*randreal()
   end do
  end do init_set
  mix: do
   mutationP=mutationP0
   i = randint(1,GA_POPSIZE)
   j = randint(1,GA_POPSIZE)
   do while (i==j)
    j = randint(1,GA_POPSIZE)
   end do
   do k=0,n-1
    a(k) = setparam(i,k)
   end do
   s1 = cost(a,x,y,n,dimen,funk)
   do k=0,n-1
    a(k) = setparam(j,k)
    ea(k)= 0.0
   end do
   s2 = cost(a,x,y,n,dimen,funk)
   fit(i)=s1
   fit(j)=s2
   if (abs(s1 - s2) <= 0.00001 .or. abs(s1)<=tol .or. abs(s2)<=tol ) then
      h = h + 1
      if ( h >= 1000 ) then
       exit mix
      end if
   end if
   renorm = s1*s2/(s1+s2)
   s1=renorm/s1
   s2=renorm/s2
   suma = 1.0/(s1+s2+mutationP)
   s1=s1*suma
   s2=s2*suma
   mutationP=mutationP*suma
   rrr = randreal()
   if(rrr<=s1)then
     do l=0,n-1 ! 2 <- 1
      setparam(j,l) = setparam(i,l)
     end do
   else if ( rrr > s1 .and. rrr<= s2 + s1 ) then
     do l=0,n-1 ! 1 <- 2
      setparam(i,l) = setparam(j,l)
     end do
   else
    k = randint(0,n-1) ! k-esima componente
    if(s1 <= s2) then  ! coste de i > coste de j
     forall (l=0:n-1)
      setparam(i,l) = setparam(j,l)
     end forall
     !setparam(i,k) = max_range*randreal()
     if(k==0)setparam(i,k) = max_range*randreal()
     if(k>=1)setparam(i,k) = randreal()
     if(k==3)setparam(i,k) = max_range*randreal()
    else
     forall (l=0:n-1)
      setparam(j,l) = setparam(i,l)
     end forall
     !setparam(i,k) = max_range*randreal()
     if(k==0)setparam(j,k) = max_range*randreal()
     if(k>=1)setparam(j,k) = randreal()
     if(k==3)setparam(j,k) = max_range*randreal()
    end if
   end if
  end do mix
  i = minloc(fit, dim=1)
  write(6,*)'Citizen Elitist:',i,'Fitness:',fit(i)
  !call sort_by_cost(q,n,setparam,fit,i)
  write(6,*)'Parameters:',(setparam(i,k),k=0,n-1), &
  'Deviations:',(sqrt( sum([((setparam(j,k)-a(k))**2,j=1,GA_POPSIZE)]) &
  / real(GA_POPSIZE-1) ),k=0,n-1)
  do k=0,n-1
   a(k) = setparam(i,k)
   ea(k)= sqrt( sum([((setparam(j,k)-a(k))**2,j=1,GA_POPSIZE)])/real(GA_POPSIZE-1) )
  end do
  return
 end subroutine fitgen
! ...
! subroutine sort_by_cost(q,n,setparam,fit,k)
! implicit none
! integer,intent(in) :: q,n
! integer            :: i
! integer,intent(out):: k
! real               :: fit(q)
! real,intent(in)    :: setparam(q,0:n-1)
! real               :: current_fit = 9999999.9
 ! ...
! do i=1,size(fit)
!  if(fit(i)<=current_fit)then
!   current_fit=fit(i)
!   k=i
!  end if
! end do
! end subroutine sort_by_cost
! ...
 real function cost(a,x,y,n,l,function_)
! calcula el coste
  implicit none
  integer,intent(in)        ::  n,l
  integer                   ::  i = 0
  real,intent(in)           ::  a(0:n-1),x(l),y(l)
  real                      ::  funk(l)
  real                      ::  x0 = 0.0
  character(100),intent(in) ::  function_
  cost = 0.0
  funk = 0.0
  do i=1,l
   x0 = x(i)
   funk(i) = model(a,x0,n,function_)
  end do
  cost = sum([(abs(y(i)-funk(i))**2,i=1,l)])
  return
 end function cost
!
 real function model(a,x,n,function_,T)
  implicit none
  integer,intent(in)        :: n
  real,intent(in)           :: a(0:n-1)
  real,intent(in)           :: x
  character(100),intent(in) :: function_
  real,intent(in),optional  :: T
  real,parameter            :: R = 0.008314472 ! kJ / mol / K
  select case (function_)
   case("freundlich")!n = a*x**b #function #model
    model = a(0)*x**a(1)
   case ("langmuir") !n = nmax*alfa*P/(1+alfa*P) #function #model
    model = a(0)*a(1)*x/(1+a(1)*x)
   case ("toth")     !n=f(x)=Nmax*alfa*x/(1+(alfa*x)**c)**(1/c) #function #model
    model = (a(0)*a(1)*x)/((1.0+(a(1)*x)**a(2))**(1/a(2)))
   case ("jensen")   !n = k1*x/(1+(k1*x/(alfa*(1+k2*x))**c))**(1/c) #function #model
    model = a(0)*x/(1+(a(0)*x/(a(1)*(1+a(3)*x))**a(2)))**(1/a(2))
   case ("dubinin_raduschkevich") ! N=Nm*exp(-(RT/Eo ln(Po/P))^2)  #model
    model = a(0)*exp(-((R*T/a(1))*log(a(2)/x) )**2)
   case ("langmuir_dualsite")      ! N=Nm*b*P/(1+b*P) + Nn*c*P/(1+c*P) #model
    model = a(0)*a(1)*x/(1+a(1)*x) + a(2)*a(3)*x/(1+a(3)*x)
   case ("dubinin_astakhov")       ! N=Nm*exp(-(RT/Eo ln(Po/P))^d) #model
    model = a(0)*exp(-((R*T/a(1))*log(a(2)/x) )**a(3))
  end select
  return
 end function model
end program iast_desarrollo_GA
