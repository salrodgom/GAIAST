module GA
 use globals
 use pseudorandoms
 use sorting
 implicit none
 private
 public                        :: FitGA, WriteModel, Fitness
 !
 integer,parameter             :: GA_Size     = 2**12           ! population size
 real,parameter                :: GA_MutationRate = 1.0/3.0     ! Mutation Rate
 real,parameter                :: GA_EliteRate= 0.1             ! Elitists population rate
 real,parameter                :: GA_MotleyCrowdRate=0.25       ! Pirates population rate
 real,parameter                :: GA_DisasterRate = 0.01      ! 
 integer,parameter             :: MaxLineLength = maxNpar*32      ! Max Lenght for a binary string
 integer,parameter             :: GA_Elitists = int(GA_size*GA_EliteRate) ! Number of elitists
 integer,parameter             :: GA_Motleists = int(GA_size*(1.0-GA_MotleyCrowdRate)) ! Number of pirates
 !
 type                          :: DNA
  character(len=MaxLineLength) :: genotype                      ! Codification in binary of the parameter
  real                         :: phenotype(1:maxNpar)          ! Float (real) representation of the genotype
  real                         :: fitness                       ! Fitness (cost function) of this DNA sequence
 end type
 type(DNA), pointer            :: parents(:)
 type(DNA), pointer            :: children(:)
 type(DNA), target             :: pop_alpha( ga_size )
 type(DNA), target             :: pop_beta( ga_size )
 contains
!
! type(DNA) function NewCitizen(ID)
! ! Done
! ! Create a new citizen in the Ensemble
!   implicit none
!   integer,intent(in) :: ID
!   integer            :: i
! ! 
!   NewCitizen%genotype = ' '
!   do i = 1,32*Component(id)%TotalNpar
!     NewCitizen%genotype(i:i) = achar(randint(48,49))
!   end do
!   do i = 1,Component(id)%TotalNpar
!     NewCitizen%phenotype(i) = bin2real( NewCitizen%genotype(32*(i-1)+1:32*i) )
!   end do
!   NewCitizen%fitness = fitness(NewCitizen%phenotype,id)
!   return
! end function NewCitizen
  type(DNA) function NewCitizen(ID)
  ! Done
  ! Create a new citizen in the Ensemble
    implicit none
    integer,intent(in) :: ID
    integer            :: i, j, Site, Param, NParam
  ! 
    NewCitizen%genotype = ' '
    i = 0
    do Site=1,Component(ID)%Sites
     j = 0
     do Param = 1, Component(ID)%Npar(Site)
      i = i + 1
      j = j + 1
      select case(Component(ID)%Models(Site))
       case("Langmuir")
        if(j==1) NewCitizen%phenotype(i) = r4_uniform(0.0,20.0)
        if(j==2) NewCitizen%phenotype(i) = r4_uniform(1E-20,1E+20)
       case("Langmuir-Freundlich","Sips")
        if(j==1) NewCitizen%phenotype(i) = r4_uniform(0.0,20.0) 
        if(j==2) NewCitizen%phenotype(i) = r4_uniform(1E-20,1E+20) 
        if(j==3) NewCitizen%phenotype(i) = r4_uniform(0.1,2.1)
       case("Toth")
        if(j==1) NewCitizen%phenotype(i) = r4_uniform(0.0,20.0)
        if(j==2) NewCitizen%phenotype(i) = r4_uniform(1E-20,1E+20)
        if(j==3) NewCitizen%phenotype(i) = r4_uniform(0.1,1.1)
       case default
        NewCitizen%phenotype(i) = r4_uniform(0.0,1.0E+20)
      end select 
     end do
    end do
    do i = 1,Component(ID)%TotalNpar
     NewCitizen%genotype(32*(i-1)+1:32*i) = real2bin(NewCitizen%phenotype(i))
     !read(NewCitizen%genotype(32*(i-1)+1:32*i),'(b32.32)') NewCitizen%phenotype(i)
    end do
    !write(6,'(6(f14.7,1x))') ( NewCitizen%phenotype(i), i =1, Component(ID)%TotalNpar )
    !write(6,'(a)') NewCitizen%genotype(1:32*Component(ID)%TotalNpar)
    NewCitizen%fitness = fitness(NewCitizen%phenotype,id)
    !STOP
    return
  end function NewCitizen
!
  subroutine UpdateCitizen( Citizen ,id)
  ! Done
  ! Updated fitness value for a citizen
    implicit none
    integer,intent(in)        :: id
    type(DNA), intent(in out) :: Citizen
    integer                   :: i
    do i = 1,Component(id)%TotalNpar
     read(Citizen%genotype(32*(i-1)+1:32*i),'(b32.32)') Citizen%phenotype(i)
    end do
    Citizen%fitness = fitness( Citizen%phenotype,ID)
    return
  end subroutine UpdateCitizen
!
  subroutine WriteCitizen(k,kk,kkk,id,lod,vgh)
  ! Done
  ! Output in screen the fitting evolution
  ! k -> Citizen
  ! kk -> Optimisation Cycle 
  ! kkk -> Diversity
  ! vgh -> Total Diversity
    implicit none
    integer,intent(in)             :: k,kk,kkk,id,lod,vgh
    integer                        :: i
    !
    write(6,'(i2,1x,i5,1x,i3,1x,a32,1x,e25.12,1x,f14.7,1x,a,1x,a,i8,a,i8,a,1x,a)') ID,kk,k,&
            parents(k)%genotype(1:32),&
            parents(k)%phenotype(1),&
            parents(k)%fitness,'[Fitness]','(',kkk,'/',vgh,')','[Similarity]' !,k
    do i=2,Component(id)%TotalNpar
      if(lod>0.and.i==3)then
        write(6,'(13x,a32,1x,e25.12,10x,a,1x,i2,a)')&
            parents(k)%genotype(1+32*(i-1):32*i),&
            parents(k)%phenotype(i),'Finishing:',lod,'/50'
      else
        write(6,'(13x,a32,1x,e25.12)')&
            parents(k)%genotype(1+32*(i-1):32*i),&
            parents(k)%phenotype(i)
      end if
    end do
  end subroutine WriteCitizen
  !
  real function Model(p,X,ID)
    Implicit None
    real,intent(in)    :: p(maxNpar)
    real,intent(in)    :: X
    integer,intent(in) :: ID
    integer            :: Site, Param, i
    real               :: tmp1, tmp2
    !
    Model = 0.0
    i = 0
    do Site=1,Component(ID)%Sites
      do Param =1, Component(ID)%Npar(Site)
        i = i + 1
        Component(ID)%p(Site,Param) = p(i)
      end do
      select case(Component(ID)%Models(Site))
      case("Langmuir")
        Model = Model + &
        Component(ID)%p(Site,1)*Component(ID)%p(Site,2)*X/&
         ( 1.0 + Component(ID)%p(Site,2)*X )
      case("Langmuir-Freundlich")
        Model = Model + &
        Component(ID)%p(Site,1)*Component(ID)%p(Site,2)*X**Component(ID)%p(Site,3)/&
         ( 1.0 + Component(ID)%p(Site,2)*X**Component(ID)%p(Site,3) )
      case("Sips")
        tmp1 = 1.0/Component(ID)%p(Site,3)
        Model = Model + &
        Component(ID)%p(Site,1)*(Component(ID)%p(Site,2)*X)**tmp1/&
         ( 1.0 + (Component(ID)%p(Site,2)*X)**tmp1 )
      case("Toth")
        tmp1 = 1.0/Component(ID)%p(Site,3)
        Model = Model + &
        Component(ID)%p(Site,1)*Component(ID)%p(Site,2)*X/&
         ( 1.0 + (Component(ID)%p(Site,2)*X)**Component(ID)%p(Site,3))**tmp1    
      case("Unilan")
        tmp1 = Component(ID)%p(Site,2)*exp(+Component(ID)%p(Site,3))
        tmp2 = Component(ID)%p(Site,2)*exp(-Component(ID)%p(Site,3))
        Model = Model + &
        (Component(ID)%p(Site,1)/(2*Component(ID)%p(Site,3)))*&
        log((1+tmp1*x)/(1+tmp2*x))
      case("Asymptotic-Temkin")
        Model = Model + Component(ID)%p(Site,1)*(Component(ID)%p(Site,2)*X/ &
         (1.0+Component(ID)%p(Site,2)*X)+Component(ID)%p(Site,3)*           &
         ((Component(ID)%p(Site,2)*X/(1.0+Component(ID)%p(Site,2)*X ))**2)* &
         (Component(ID)%p(Site,2)*X/(1.0+Component(ID)%p(Site,2)*X )-1.0) )
      case default
        stop 'Model Not Found'
      end select
    end do 
    return
  end function Model
  !
  subroutine WriteModel(p,ID)
    implicit none
    real,intent(in)     :: p(1:maxNpar)
    integer,intent(in)  :: ID
    integer             :: u, i, Site
    character(len=20)   :: FileName
    write(filename,'(a,i1,a)') 'model.',Component(ID)%ID,'.dat'
    filename = trim(adjustl(filename))
    open(newunit=u,file=filename)
    !
    write(u,'(a,1x,a,1x,a)')'# Parameters for',Component(ID)%Label,':'
    do Site=1, Component(ID)%Sites 
      write(u,'(a,i1,a,1x,100(e24.12,1x))') '# Site',Site,':',&
           (Component(ID)%p(Site,i),i=1,Component(ID)%Npar(Site))
    end do
    write(u,'(a)')'#'
    do i = 1, Component(ID)%NPressures
      write(u,*) Component(ID)%PurePressure(i), Component(ID)%PureLoading(i), &
                 Model(p,Component(ID)%PurePressure(i),ID)
    end do
    close(u)
  end subroutine WriteModel
  !
  real function Fitness(p,ID)
    ! Cost function
    implicit none
    real, intent(in)    :: p(maxNpar)
    integer,intent(in)  :: ID
    real, parameter     :: PenaltyCost = 50.0
    real                :: Penalty
    real                :: xx, yy 
    integer             :: i, Site, Param, NParam
    Penalty = 0.0
    !
    i = 0
    do Site=1,Component(ID)%Sites 
      do Param =1, Component(ID)%Npar(Site)
        i = i + 1
        Component(ID)%p(Site,Param) = p(i)
      end do
      ! 
      phys_constrains: if ( physical_constrains_flag ) then
        ! We check the parameters for each model in each Site
        select case(Component(ID)%Models(Site))
        case("Langmuir")
         ! constrains:
         ! a>0, b>0
         if(Component(ID)%p(Site,1) < 0.0 .or.&
            Component(ID)%p(Site,2) < 0.0 ) then
           penalty = PenaltyCost
           exit phys_constrains
         end if    
        case("Langmuir-Freundlich")
        ! a >0, b >0, 0<c<1
         if(Component(ID)%p(Site,1) < 0.0 .or.&
            Component(ID)%p(Site,2) < 0.0 .or.&
            Component(ID)%p(Site,3) < 0.0 ) then 
           !Component(ID)%p(Site,3) > 1.0 ) then ! <- Check! 
           penalty = PenaltyCost
           exit phys_constrains            
         end if
        case("Sips")
        ! a >0, b >0, 0<c<1
         if(Component(ID)%p(Site,1) < 0.0 .or.&
            Component(ID)%p(Site,2) < 0.0 .or.&
            Component(ID)%p(Site,3) < 0.0 ) then 
           !Component(ID)%p(Site,3) > 1.0 ) then ! <- Check! 
           penalty = PenaltyCost
           exit phys_constrains            
         end if
        case("Toth")
         ! a >0, b >0, 0<c<1
         if(Component(ID)%p(Site,1) < 0.0 .or.&
            Component(ID)%p(Site,2) < 0.0 .or.&
            Component(ID)%p(Site,3) < 0.0 ) then 
           !Component(ID)%p(Site,3) > 1.0 ) then ! <- Check! 
           penalty = PenaltyCost
           exit phys_constrains            
         end if 
        case("Unilan")
        ! 
         if(Component(ID)%p(Site,1) < 0.0 .or.&
            Component(ID)%p(Site,2) < 0.0 .or.&
            Component(ID)%p(Site,3) < 0.0 ) then 
           penalty = PenaltyCost
           exit phys_constrains            
         end if 
        case("Asymptotic-Temkin")
         if(Component(ID)%p(Site,1) < 0.0 .or.&
           Component(ID)%p(Site,2) < 0.0 .or.&
           Component(ID)%p(Site,3) < 0.0 ) then 
          penalty = PenaltyCost
          exit phys_constrains            
         end if 
        case default
          Stop 'Model Not Found in Fitness function'
        end select
      end if phys_constrains
    end do
    NParam  = i
    Fitness = Penalty
    !
    do i = 1, Component(ID)%NPressures
      xx = Component(ID)%PurePressure(i)
      yy = Component(ID)%PureLoading(i)
      fitness = fitness + (yy - model(P,xx,ID))**2
    end do
    fitness = sqrt( fitness / (Component(ID)%NPressures-NParam) )
    if(isnan(fitness)) fitness      = 99999999.999999
    if(fitness==0.0000000000)fitness= 99999999.999999
    return
  end function Fitness
!
  pure real function RCorrelation()
   implicit none
   RCorrelation = 0.0
   return
  end function RCorrelation
!
  pure real function IntegrateArea(ID)
    implicit none
    integer,intent(in)   :: ID
    integer              :: i
    real                 :: delta,x
    real                 :: factor
    !
    x = MinPressure
    delta=(MaxPressure-MinPressure)/(Component(ID)%NPressures)
    area: do i = 0,Component(ID)%NPressures-1
      factor = 1.0
      if (i==0.or.i==Component(ID)%NPressures-1) factor = 3.0/8.0
      if (i==1.or.i==Component(ID)%NPressures-2) factor = 7.0/6.0
      if (i==2.or.i==Component(ID)%NPressures-3) factor = 23.0/24.0
      x = x + delta*(0.5+i)
      IntegrateArea = IntegrateArea + factor*delta*Component(ID)%PurePressure(i+1)/x
    end do area
    return
  end function integrateArea
!
  subroutine SortByFitness()
   type(DNA)             ::  sorted(1:GA_size)
   integer               ::  k,i
   real(16)              ::  ftnss(1:ga_size)
   do k=1,GA_size
    ftnss(k)=dble(parents(k)%fitness)
    if(isnan(parents(k)%fitness)) ftnss(k) = 9999999999.d99
   end do
   call QsortC( ftnss )
   exter:do k=1,GA_size ! <- ordered
    inter:do i=1,GA_size
    if( dble(parents(i)%fitness) == ftnss(k))then
      sorted(k) = parents(i)
      cycle inter
    end if
    end do inter
   end do exter
   parents=sorted
   return
  end subroutine SortByFitness

  integer function Biodiversity( ID, Citizen, Total)
    implicit none
    integer,intent(in)             :: ID
    type(DNA), intent(in)          :: Citizen(1:GA_size)
    integer,intent(out)            :: Total
    integer                        :: i, j, k
    integer                        :: cont
    integer                        :: tmp1
    character(len=20)              :: mode = 'Normal'
    real                           :: error_d = 1e-3
    select case (mode)
    case('None')
      Biodiversity = 0
    case('Normal')
      Total = 0
      Biodiversity = 0
      ! Expensive
      do k =1,GA_size
        do j=k+1,GA_size
          Total = Total + 1
          tmp1 = Component(ID)%TotalNpar
          ! characters
          if( Citizen(k)%genotype(1:32*tmp1) == Citizen(j)%genotype(1:32*tmp1) ) &
            Biodiversity = Biodiversity + 1
        end do
      end do
    case('Superficial')
      Biodiversity = 0.0
      Total=0
      ! Very Expensive
      do k = 1, GA_size
        do j = k+1, GA_size
          cont=0
          Total = Total + 1
          tmp1 = Component(ID)%TotalNpar
          dbio: do i = 1,tmp1
            if(abs(Citizen(k)%phenotype(i)-Citizen(j)%phenotype(i)) <= error_d ) &
              cont = cont + 1
          end do dbio
        if( cont == tmp1 ) Biodiversity = Biodiversity + 1
        end do
      end do
    end select
    return
  end function Biodiversity
!
  subroutine NuclearDisaster(ID)
  ! "There can only be one".
    implicit none
    integer,intent(in) :: ID
    integer            :: i
    !
    do i = 2, GA_Size
      children(i) = NewCitizen(ID)
    end do
    !STOP 'Nuclear Disaster!'
    return
  end subroutine NuclearDisaster
!
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
!
  subroutine Elitism()
    children(:GA_ELITISTS) = parents(:GA_ELITISTS)
    return
  end subroutine
!
  subroutine Mutate( Mutant , ID )
  ! Cowabunga!
  !--------------------------------------------
  !     *
  ! 0000000000                 ->    0000100000
  !--------------------------------------------
    implicit none
    type(DNA), intent(inout)    :: Mutant
    integer,intent(in)          :: ID
    integer                     :: i, ipos
    do i = 1,Component(id)%TotalNpar
      ipos = randint(32*(i-1)+1,32*i)
      Mutant%genotype(ipos:ipos) = achar(randint(48,49))
    end do
    return
  end subroutine Mutate
!  
  subroutine Crossover(ID,s1,s2,i1,i2,j1,j2)
! [s1:s2] range of children
! [i1:i2] range of parent1
! [j1:j2] range of parent2
!----------------------------------------------
!  parent1    parent2                children
!    *          *                      *
!  00|000000  11|111111        ->    00|111111
!----------------------------------------------
   implicit none
   integer,intent(in) :: ID
   integer,intent(in) :: s1,s2,i1,i2,j1,j2
   integer            :: k1,k2,spos,i
   ! loop for each partner1
   do i=s1,s2
    call choose_randomly(i1,i2,j1,j2,k1,k2)
    spos = randint(0, 32*Component(id)%TotalNpar )
    children(i)%genotype = parents(k1)%genotype(:spos) // parents(k2)%genotype(spos+1:)
   end do
   return
  end subroutine Crossover
!
  subroutine choose_randomly(kk1,kk2,jj1,jj2,ii1,ii2)
  ! Choose ii1 in [kk1:kk2] range
  ! Choose ii2 in [jj1:jj2] range
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
  subroutine Mate(ID)
! Mate the genes
    implicit none
    integer, intent(in) :: ID
    integer             :: i
   ! Retain the first 25% of the childrens
    call Elitism()
    call Crossover(ID,GA_Elitists+1,GA_Motleists,1,GA_Elitists,1,GA_Elitists)
    call Crossover(ID,GA_Motleists+1,GA_Size-100,1,GA_Elitists,GA_Elitists+1,GA_size)
    ! Replace the last GA_Elitists of the childrens by new childrens
    children(GA_size-GA_Elitists:GA_size) = NewCitizen(ID)
    do i = GA_Elitists+1, GA_Size-GA_Elitists
    ! Mutation:
    ! I'll respect the first 10% of the childrens of the mutations.
      if( i >= GA_elitists*0.1 .and. r4_uniform(0.0,1.0) < GA_MutationRate) call Mutate(children(i),ID)
      call UpdateCitizen(children(i),ID)
    end do
    if(r4_uniform(0.0,1.0) < GA_DisasterRate ) call NuclearDisaster(ID)
   return
  end subroutine Mate
!
  subroutine FitGA(ID)
    implicit none
    integer,intent(in)     :: ID
    integer                :: Citizen
    !
    integer                :: OptimisationStep = 0
    integer,parameter      :: MaxOptimisationStep = 500
    integer                :: FullfilledConditionStep = 0
    integer,parameter      :: MaxFullfilledConditionStep = 50
    real, parameter        :: minimum_fitness = 5.0e-1
    real                   :: TempFitnessValue = 0.0
    integer                :: TempVarietyValue = 0.0
    integer                :: TotalBio = 0
    real,parameter         :: tolerance_equal_fitness = 1e-5
!
    integer                :: Site, Param, s
    integer,parameter      :: minstep = 10
    !
    TempFitnessValue = 999.000
    FullfilledConditionStep = 0; OptimisationStep = 0
    ! The 
    pop_alpha = [(NewCitizen(ID), Citizen = 1, GA_size)]
    parents =>  pop_alpha
    children => pop_beta
    !
    if ( Refitting_Flag ) then
      write(6,'(a)')'Refitting Activated'
      do Citizen = 1, GA_Elitists
        parents(Citizen)%genotype=string_IEEE(ID)
        children(Citizen)%genotype=string_IEEE(ID)
        call UpdateCitizen(parents(Citizen),ID)
        call UpdateCitizen(children(Citizen),ID)
      end do
      call Mate(ID)
      call Swap()
      call SortByFitness()
    end if
    ! Start the Fitting procedure
    converge: do while ( .true. )
      OptimisationStep = OptimisationStep + 1
      if ( OptimisationStep == 1 ) write(6,'(a)') 'Go Down the Rabbit Hole >'
      call SortByFitness()
      TempVarietyValue = Biodiversity( ID, children, TotalBio )
      do Citizen=1,1
        call WriteCitizen(Citizen,OptimisationStep,TempVarietyValue,ID,&
                          FullfilledConditionStep,TotalBio)
      end do
      if( OptimisationStep >= minstep .and.           &
          parents(1)%fitness <= minimum_fitness .and. &
          abs(parents(1)%fitness-TempFitnessValue) <= tolerance_equal_fitness ) then
        FullfilledConditionStep = FullfilledConditionStep + 1
      else
        FullfilledConditionStep = 0
      end if
    ! System converges!
      if ( OptimisationStep >= MaxOptimisationStep .or. &
         FullfilledConditionStep >= MaxFullfilledConditionStep ) exit converge
      ! Pairing 
      call Mate(ID)
      call Swap()  
      ! Take the best fitness value:
      TempFitnessValue = parents(1)%fitness
    end do converge
!   Output
    s = 0
    do Site =1, Component(ID)%Sites
      do Param = 1, Component(ID)%Npar(Site)
        s = s + 1
        Component(ID)%p(Site,Param) = Children(1)%phenotype(s)
      end do
    end do
    call WriteModel(Children(1)%phenotype,ID)
    return
  end subroutine FitGA
end module GA
