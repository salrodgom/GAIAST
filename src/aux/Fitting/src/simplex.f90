module simplex_
 use globals
 use ga
 implicit none
 private
 public             :: FitSIMPLEX
 contains
! ======================================================================
! Interface between program and module:
subroutine FitSIMPLEX(ID)
 implicit none
 integer,intent(in) :: ID
 real               :: e = 1.0e-4
 real               :: scale = 1.0
 integer            :: iprint = 0
 integer            :: Site, Param, i
 real               :: sp(0:Component(ID)%TotalNpar-1)
 real               :: p(1:maxNpar)
 !
 write(6,'(a)')' '
 write(6,'(a)')'Minimising the cost function using the Nelder-Mead SIMPLEX method:'
 write(6,*)'# Compound:',ID
 i = 0
 do Site = 1,Component(ID)%Sites
   do Param = 1,Component(ID)%Npar(Site)
     sp(i) = Component(ID)%p(Site,Param)
     i = i + 1
   end do 
 end do
 call SIMPLEX(sp,ID,Component(ID)%TotalNpar,e,scale,iprint)
 i = 0
 do Site = 1,Component(ID)%Sites
   do Param = 1,Component(ID)%Npar(Site)
     Component(ID)%p(Site,Param) = sp(i)
     p(i+1) = sp(i)
     i = i + 1
   end do 
 end do 
 call WriteModel(p,ID)
end subroutine FitSIMPLEX

! ======================================================================
! This is the function to be minimized
real function func(n,x,compound) result(rosen)
  implicit none
  integer,intent(in)   :: n
  real,   intent (in)  :: x(0:n-1)
  integer,intent(in)   :: compound
  real                 :: sp(1:Component(Compound)%TotalNpar)
  integer              :: i
  sp = 0.0
  do i=0,Component(Compound)%TotalNpar-1
   sp(i+1) = x(i)
  end do
  rosen = fitness(sp,compound)
  return
end function func
! 
subroutine SIMPLEX(start, compound, n, EPSILON, scale, iprint)
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
  real :: min,fsum,favg,cent

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
!250  FORMAT(A29,F7.4)
!300  FORMAT(F11.6,F11.6,F11.6)
  return
  end subroutine simplex
! ======================================================================
end module simplex_
