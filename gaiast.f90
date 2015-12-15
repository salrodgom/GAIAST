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

 integer function randint(i,j,seed)
  real               ::  a,b
  integer,intent(in) ::  i,j,seed 
  a = real(i)
  b = real(j)
  randint=int(r4_uniform(a,b+1.0,seed))
 end function randint

 REAL function r4_uniform2(b1,b2,seed)
  implicit none
  real b1,b2
  integer i4_huge,k,seed
  parameter (i4_huge=2147483647)
  if(seed == 0) then
   write(*,'(b1)')' '
   write(*,'(b1)')'R4_UNIFORM - Fatal error!'
   write(*,'(b1)')'Input value of SEED = 0.'
   stop '[ERROR] Chiquitan chiquititan tan tan &
    Que tun pan pan que tun pan que tepe tepe &
    Pan pan pan que tun pan que pin '
  end if
  k=seed/127773
  seed=16807*(seed-k*17773)-k*2836
  if(seed<0) then
    seed=seed+i4_huge
  endif
  r4_uniform=b1+(b2-b1)*real(dble(seed)* 4.656612875D-10)
  return
 end function r4_uniform2

 real function r4_uniform(b1,b2,seed)
  implicit none
  real,intent(in)    :: b1,b2
  integer,intent(in) :: seed
  r4_uniform = b1+(b2-b1)*iast_rand(seed,1)
  return
 end function r4_uniform

 real function iast_rand(iseed,imode)
!  imode = 1 => random number between 0 and 1
!        = 2 => random number between -1 and 1
  implicit none
  integer  :: iseed
  integer  :: imode
  integer       :: i
  integer, save :: ia1 = 625
  integer, save :: ia2 = 1741
  integer, save :: ia3 = 1541
  integer, save :: ic1 = 6571
  integer, save :: ic2 = 2731
  integer, save :: ic3 = 2957
  integer, save :: ix1
  integer, save :: ix2
  integer, save :: ix3
  integer       :: j
  integer, save :: m1  = 31104
  integer, save :: m2  = 12960
  integer, save :: m3  = 14000
  real          :: rm1
  real          :: rm2
  real,    save :: rr(97)
!
  rm1 = 1.0/m1
  rm2 = 1.0/m2
  if (iseed.lt.0) then
    ix1 = mod(ic1 - iseed,m1)
    ix1 = mod(ia1*ix1 + ic1,m1)
    ix2 = mod(ix1,m2)
    ix1 = mod(ia1*ix1 + ic1,m1)
    ix3 = mod(ix1,m3)
    do i = 1,97
      ix1 = mod(ia1*ix1 + ic1,m1)
      ix2 = mod(ia2*ix2 + ic2,m2)
      rr(i) = (float(ix1) + float(ix2)*rm2)*rm1
    enddo
    iseed = abs(iseed)
  endif
  ix3 = mod(ia3*ix3 + ic3,m3)
  j = 1 + (97*ix3)/m3
  if (j.gt.97.or.j.lt.1) then
    write(6,*)'array bounds exceeded in IAST_random'
    stop 'iast_random'
  endif
  if ( imode == 1 ) then
    iast_rand = rr(j)
  else
    iast_rand = 2.0*rr(j) - 1.0
  endif
  ix1 = mod(ia1*ix1 + ic1,m1)
  ix2 = mod(ia2*ix2 + ic2,m2)
  rr(j) = (float(ix1) + float(ix2)*rm2)*rm1
  return
  end function iast_rand
end module

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
  real, intent(in out), dimension(:) :: A
  integer                            :: iq
  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
 end subroutine QsortC
 subroutine Partition(A, marker)
  real, intent(in out), dimension(:) :: A
  integer, intent(out)               :: marker
  integer                            :: i, j
  real                               :: temp
  real                               :: x
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
 integer                    :: npar,ncomponents,i,seed
 integer,allocatable        :: np(:)
 integer                    :: err_apertura,ii,intervalos,j
 real                       :: tol = 0.001, tolfire = 0.25
 real,allocatable           :: param(:,:),concy(:),pi(:,:),iso(:,:)
 integer,allocatable        :: npress(:)
 integer,parameter          :: maxnp = 10, maxcompounds = 10, maxdata = 1000
 real,target                :: datas(2,maxcompounds,maxdata),f(maxdata)
 real,pointer               :: x(:),y(:),alldat(:,:)
 character(100),allocatable :: ajuste(:)
 character(100)             :: line,string,intmethod
 character(5)               :: inpt
 logical                    :: flag = .true., FlagFire = .false.
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
     write(6,*)'Fit is already done'
     readparameter: do ii=1,ncomponents
      if(ncomponents>=3) stop "[ERROR] Uh Oh: you're screwed"
      read(5,*)(param(ii,j),j=0,np(ii)-1)
     end do readparameter
    end if
    cycle read_input_do
   end if
   if(line(1:5)=='inter') then
     read(line,*)inpt,intervalos
     allocate(pi(ncomponents,0:intervalos),iso(ncomponents,0:intervalos))
   end if
   if(line(1:5)=='toler') read(line,*)inpt, Tol
   if(line(1:5)=='tempe') read(line,*)inpt, T
   if(line(1:5)=='InteM') read(line,*)inpt, IntMethod
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
  superior = 0.0
  inferior = 1.0e12
  do ib = 1,ncomponents
   do jb = 1,npress(ib)
    if( datas(1,ib,jb) >= superior) superior = datas(1,ib,jb)
    if( datas(1,ib,jb) <= inferior) inferior = datas(1,ib,jb)
   end do
  end do
! ...
  superior = log10( superior )
  inferior = log10( inferior )
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
    xxx(compound) = 10**( inferior + ib*dx )            ! real space
    x0=xxx(compound)
    yyy(compound) = model(apar, np(compound), x0, funk )
!   ...
    y0 = 10**( inferior + (ib+1)*dx )
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

 subroutine IAST()
  implicit none
  integer        :: i,j
  real           :: comp,dx,presion1,presion2,concx1,concx2,p,n1,n2,n

  dx = real((superior-inferior)/intervalos)

  open(unit=104,file='adsorcion.dat')
  do i=0,intervalos-1
   do j=0,intervalos-1
    comp = abs(pi(1,i)-pi(2,j))
    if (comp <= tol) then
     presion1 = 10**(inferior+i*dx)
     presion2 = 10**(inferior+j*dx)
     concx1 = presion2*concy(1) / (concy(2)*presion1 + concy(1)*presion2)
     concx2 = 1.0 - concx1
     p  = presion1*concx1/concy(1)
     n1 = iso(1,i)*concx1
     n2 = iso(2,j)*concx2
     n = n1 + n2
     write(104,*)p,n1,n2,n,concy
    end if
   end do
  end do
  close(104)
  return
 end subroutine IAST

 real function integrate(x0,x1,integration_points,a,n,funk)
  implicit none
  integer              ::  i
  integer,intent(in)   ::  n,integration_points
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
  if(err_apertura/=0) stop '[ERROR] isotermaN.dat :  Catastrophic failure'
  do1: do i = 1,npress(ii)
   read (456,'(A)',iostat=err_apertura) line
   if(err_apertura/=0) exit do1
   read(line,*)datas(1,ii,i),datas(2,ii,i)
   write(6,*) datas(1,ii,i),datas(2,ii,i)
  end do do1
  close(456)
 end do nisothermpure1
 end subroutine ReadIsotherms
! ...
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
  case ("jensen")
   n=4
  case ("dubinin_raduschkevich")
   n=3
  case ("langmuir_dualsite")
   n=4
  case ("langmuir_freundlich_dualsite")
   n=6
  case ("dubinin_astakhov")
   n=4
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
   case("freundlich")
    model = a(0)*xx**a(1)
   case("langmuir")
    model = a(0)*a(1)*xx/(1+a(1)*xx)
   case("langmuir_freundlich")
    model = a(0)*a(1)*xx**a(2)/(1+a(1)*xx**a(2))
   case("toth") ! f(x)=Nmax*alfa*x/(1+(alfa*x)**c)**(1/c)
    model = a(0)*a(1)*xx/((1.0+(a(1)*xx)**a(2))**(1.0/a(2)))
   case("langmuir_dualsite")
    model = a(0)*a(1)*xx/(1+a(1)*xx) + a(2)*a(3)*xx/(1.0+a(3)*xx)
   case("langmuir_freundlich_dualsite")
    model = a(0)*a(1)*xx**a(2)/(1+a(1)*xx**a(2))+a(3)*a(4)*xx**a(5)/(1.0+a(4)*xx**a(5))
   case("jensen")!f(x)=K1*x/(1+(K1*x/(alfa*(1+k2*x))**c))**(1/c)
    model = a(0)*xx/(1.0+(a(0)*xx/(a(1)*(1+a(2)*xx))**a(3)))**(1.0/a(3)) 
   !model = a(0)*xx/(1+(a(0)*xx/(a(1)*(1+a(2)*xx))**a(3)))**(1.0/a(3))
   case ("dubinin_raduschkevich") ! N=Nm*exp(-(RT/Eo ln(Po/P))^2)  #model
    model = a(0)*exp(-((R*T/a(1))*log(a(2)/xx) )**2)
   case ("dubinin_astakhov")       ! N=Nm*exp(-(RT/Eo ln(Po/P))^d) #model
    model = a(0)*exp(-((R*T/a(1))*log(a(2)/xx) )**a(3))
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
 integer,parameter        :: ga_size     = 2**13 ! numero de cromosomas
 real,parameter           :: ga_mutationrate = 0.3333 !2000/real(ga_size) ! ga_mutationrate=0.333
 real,parameter           :: ga_eliterate= 0.25, GA_DisasterRate = 0.0000001
 integer,parameter        :: ga_elitists = int( ga_size * ga_eliterate)
 type                     :: typ_ga
  character(len=32*maxnp) :: genotype
  real                    :: phenotype(1:maxnp)
  real                    :: fitness
 end type
 type(typ_ga), pointer    :: parents(:)
 type(typ_ga), pointer    :: children(:)
 type(typ_ga), target     :: pop_alpha( ga_size )
 type(typ_ga), target     :: pop_beta( ga_size )
 contains

  type(typ_ga) function new_citizen(compound,seed)
   implicit none
   integer :: i,compound,seed
   new_citizen%genotype = ' '
   do i = 1, 32*np(compound)
    new_citizen%genotype(i:i) = achar(randint(48,49,seed))
   end do
   do i = 0,np(compound)-1
    new_citizen%genotype(i*32+1:i*32+1) = '0'
   end do
   do i = 1,np(compound)
    read(new_citizen%genotype(32*(i-1)+1:32*i),'(b32.32)') new_citizen%phenotype(i)
   end do
   new_citizen%fitness = fitness( new_citizen%phenotype,compound)
  end function new_citizen

  subroutine UpdateCitizen( axolotl ,compound )
   implicit none
   integer                     :: compound,i
   type(typ_ga), intent(inout) :: axolotl
   do i = 1,np(compound)
    read(axolotl%genotype(32*(i-1)+1:32*i),'(b32.32)') axolotl%phenotype(i)
   end do
   axolotl%fitness = fitness( axolotl%phenotype,compound)
   return
  end subroutine UpdateCitizen
 
  real function Fitness(phenotype,compound)
   implicit none
   real, intent(in)    :: phenotype(maxnp)
   integer             :: i,compound,k = 0
   real                :: a(0:np(compound)-1),xx,yy
   character(len=100)  :: funk
   logical             :: flagzero = .false.
   real                :: infinite = HUGE(40843542)
   funk = ajuste(compound)
   do i = 0,np(compound)-1
    a(i) = phenotype(i+1)
    !if (a(i) == 0.0) k=k+1
   end do
   !if(k==np(compound)) STOP
   fitness = 0.0
   do i = 1, npress(compound)
    xx = datas(1,compound,i)
    yy = datas(2,compound,i)
    fitness = fitness + 0.5*( yy - model(a,np(compound),xx,funk) )**2
   end do
   !if( fitness == 0.0 ) fitness = infinite
   !else
   ! fitness = infinite
   !end if
   return
  end function Fitness

  subroutine WriteCitizen(k,kk,kkk,compound)
   implicit none
   integer                        :: k,kk,i,compound
   character(len=100)             :: fmt_
   character(len=32*np(compound)) :: wnowaste
   real                           :: wnowasteparam(1:32*np(compound)),wfitness,kkk
   !do i=1,32*np(compound)
   ! wnowaste(i:i)=' '
   !end do
   ! ...
   !wnowaste(1:32*np(compound))=parents(k)%genotype(1:32*np(compound))
   wnowaste = parents(k)%genotype
   wnowasteparam = parents(k)%phenotype
   !do i=1,np(compound)
   ! wnowasteparam(i)=parents(k)%phenotype(i)
   !end do
   wfitness = parents(k)%fitness
   write(6,'(i2,1x,i5,1x,a32,1x,e25.12,1x,f20.10,1x,f20.10)')compound,kk,wnowaste(1:32),wnowasteparam(1),wfitness,kkk
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
   real                     ::  ftnss(1:ga_size)
   do k=1,ga_size
    ftnss(k)=parents(k)%fitness
    if(isnan(parents(k)%fitness)) ftnss(k) = 9999999999.99
   end do
   !call PIKSRT(ga_size,ftnss)
   call QsortC( ftnss )
   exter:do k=1,ga_size
    inter:do i=1,ga_size
    if( parents(i)%fitness== ftnss(k))then
      sorted(k) = parents(i)
      cycle inter 
    end if
    end do inter
   end do exter
   parents=sorted
  end subroutine SortByFitness
  
  real function Biodiversity( compound )
   implicit none
   integer,intent(in)             :: Compound
   character(len=32*np(compound)) :: str1,str2
   integer                        :: k
   character(len=20)              :: mode = 'Superficial'
   select case (mode)
    case('Superficial')
     do k = 1, ga_size
      do j = k+1,ga_size
        if(children(k)%genotype == children(j)%genotype) then
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
   do i = 0,np(Compound)
    macrophage%genotype(i*32+1:i*32+1) = '0'
   end do
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
    do j = 0,np(compound)-1
     Children%genotype(j*32+1:j*32+1) = '0'
    end do
    call UpdateCitizen(Children(i),Compound)
   end do
  end subroutine NuclearDisaster

  subroutine Swap()
   if (associated(parents, target=pop_alpha)) then
       parents => pop_beta
       children => pop_alpha
   else
       parents => pop_alpha
       children => pop_beta
   end if
  end subroutine Swap

  subroutine Elitism()
   children(:GA_ELITISTS) = parents(:GA_ELITISTS)
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
  end subroutine Mate

  subroutine choose_randomly(j1,j2)
   implicit none
   integer,intent(out) :: j1,j2
   j1  = randint(1, int(ga_size/3),seed)
   j2  = randint(1, int(ga_size/3),seed)
   do while ( j1 == j2 ) 
    j2 = randint(1, int(ga_size/3),seed)
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
   ! write(6,*)prop(1:)
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
   integer,parameter  :: maxstep = 1000, minstep = 10
   integer            :: kk, ii, i, k
   real               :: eps = 0.0, diff = 0.0, fit0 = 0.0
   kk = 0
   ii = 0
   pop_alpha = [(new_citizen(compound,seed), i = 1,ga_size)]
   parents =>  pop_alpha
   children => pop_beta
   call WriteCitizen(1,ii,eps,compound)
   converge: do while (.true.)
    ii=ii+1
    call SortByFitness()
    diff = eps
    eps = Biodiversity( compound )
    call WriteCitizen(1,ii,eps,compound)
    fire: if ( FlagFire ) then
     if ( ii >= minstep .and. parents(1)%fitness <= TolFire ) exit converge
    else
     if( abs(diff - eps) <= 0.1 .and. ii >= minstep .and. &
      parents(1)%fitness - fit0 == 0 ) then
      kk = kk + 1 
     else
      kk = 0
     end if
     if ( ii >= maxstep .or. kk >= 5 ) exit converge
    end if fire
    call Mate(compound)
    call Swap()
    fit0 = parents(1)%fitness
   end do converge
   do i = 0, np( compound )-1
    param( compound,i ) = children(1)%phenotype(i+1)
   end do
   write(111,*)'#',(param(compound,i),i=0,np(compound )-1)
   write(111,*)'#','Fitness:',fit0,'Biodiversity:',eps
   return
  end subroutine fit 
  
end module mod_genetic

program main
 use mod_random
 use gaiast_globals
 use mod_genetic
 call init_random_seed(seed)
 call read_input()
 call ReadIsotherms()
 if(flag)then
  open(111,file="iso.dat")
  do ii = 1, ncomponents
   write(6,'("Fitting compound:",1x,i2)') ii
   write(6,'(a14,1x,a24,1x,a24,1x,a14,1x,a14)')'Compound/Step:','Cromosome:','Parameters','Fitness','Diversity'
   call fit(ii,seed)
  end do
 end if
 call bounds()
 call IAST()
 close(111)
 stop 'Hu ohgc!! SSSssss' 
end program main
