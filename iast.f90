PROGRAM iastdesarrollo
 IMPLICIT NONE
!========================================================================
! VARIABLES
!========================================================================
 REAL::comp,tol,concx1,concx2,concy1,concy2,h1,h2,p,presion1,presion2,presiontotal,conc1,conc2,conc3,n,n1,n2
 REAL::a1,a2,inferior,superior,punto,s1,s2
 INTEGER::err_apertura = 0,i,k,l,j,intervalos
 CHARACTER(100):: line,ajuste,divisiones1,divisiones2,archivo1,archivo2
 REAL,ALLOCATABLE::x1(:),y1(:),x2(:),y2(:),coef1(:),coef2(:),PI1(:),PI2(:),area1(:),area2(:)
 REAL,ALLOCATABLE::e_coef1(:),e_coef2(:)
 REAL,ALLOCATABLE::puntos1(:),puntos2(:),funcion1(:),funcion2(:)

!========================================================================
! PARAMETROS
!========================================================================
 !WRITE(6,*)"AJUSTE"
 READ(5,*)ajuste
 !WRITE(6,*)"INTERVALOS"
 READ(5,*)intervalos
 !WRITE(6,*)"TOLERANCIA"
 READ(5,*)tol
 !WRITE(6,*)"CONCENTRACION GAS COMPONENTE 1"
 READ(5,*)concy1
 !WRITE(6,*)"CONCENTRACION GAS COMPONENTE 2"
 READ(5,*)concy2
! ajuste = "toth"    !seleccion del ajuste 
! intervalos = 10000 !numero de puntos para calcular el area
! tol = 0.001        !tolerancia para igualar areas
! concy1 = 0.5       !concentracion fase gas componente1
! concy2 = 0.5       !concentracion fase gas componente2
!========================================================================
! LECTURA ISOTERMAS 
!========================================================================
 !abrimos el archivo con las isotermas
 OPEN(UNIT=100,FILE="isoterma1.dat",STATUS='old', IOSTAT=err_apertura)
      !inicializamos el contador
      k=0 
      ! esta foma de escribir los bucles nos permite dar nombre a cada bucle
      ! este bucle nos cuenta el numero de datos que tenemos en el archivo
      do0: DO  
           !con '(A)' le decimos que lo que va a leer es un caracter
           READ (100,'(A)',IOSTAT=err_apertura) line 
           IF( err_apertura /= 0 ) EXIT do0
           k=k+1
      END DO do0
      ! dimensionamos los vectores de datos con el numero de datos en el archivo
      ALLOCATE(x1(1:k),y1(1:k))
      ! rebobinamos el archivo para leer los datos
      REWIND(100)

      ! leemos los datos del archivo y los escribimos en los valores
      do1: DO i=1,k
           READ (100,'(A)',IOSTAT=err_apertura) line
           IF(err_apertura/=0) EXIT do1
           READ(line,*) x1(i),y1(i)
      END DO do1
 CLOSE(100)
!-----------------------------------------------------------------------
 OPEN(UNIT=100,FILE="isoterma2.dat",STATUS='old', IOSTAT=err_apertura)
     l=0
     do2: DO
          READ (100,'(A)',IOSTAT=err_apertura) line
          IF( err_apertura /= 0 ) EXIT do2
          l=l+1
     END DO do2
     ALLOCATE(x2(1:l),y2(1:l))
     REWIND(100)
     do3: DO i=1,l
          READ (100,'(A)',IOSTAT=err_apertura) line
          IF(err_apertura/=0) EXIT do3
          READ(line,*) x2(i),y2(i)
     END DO do3
 CLOSE(100)

!========================================================================
! AJUSTE ISOTERMAS
!========================================================================

! utilizamos gnuplot para hacer el ajuste y recuperamos los valores a traves del archivo coeficientes
!!!!mantener el archivo ajustelangmuir en la misma carpeta donde se corra el programa
 SELECT CASE (ajuste)
        case('freundlich')
         allocate(coef1(0:1),coef2(0:1))
         ALLOCATE(e_coef1(0:1),e_coef2(0:1))
          !-> isoterma1   
              CALL SYSTEM ('cp isoterma1.dat isotermaN.dat')
              CALL SYSTEM ('nohup ./ajustefreundlich.gnu')
              call system('mv ajuste.eps ajuste1.eps')
              call sleep(2)
!a               = 2.10888          +/- 0.06313      (2.993%)
!b               = 2.71054          +/- 0.2272       (8.381%)
              CALL SYSTEM ("grep 'a               =' fit.log  | tail -n1 | awk '{print $3,$5}' >  coeficientes")
              CALL SYSTEM ("grep 'b               =' fit.log  | tail -n1 | awk '{print $3,$5}' >> coeficientes")
              OPEN(unit=100,file="coeficientes",status='old')
                   READ(100,*)coef1(0),e_coef1(0)  !nmax
                   READ(100,*)coef1(1),e_coef1(1)  !alfa
              CLOSE(100)
              coef1(1)=1.0/coef1(1)
              !CALL SYSTEM ('rm coeficientes fit.log isotermaN.dat')
          !-> isoterma2
              CALL SYSTEM ('cp isoterma2.dat isotermaN.dat')
              CALL SYSTEM ('nohup ./ajustefreundlich.gnu')
              call system('mv ajuste.eps ajuste2.eps')
              call sleep(2)
              CALL SYSTEM ("grep 'a               =' fit.log  | tail -n1 | awk '{print $3,$5}' >  coeficientes")
              CALL SYSTEM ("grep 'b               =' fit.log  | tail -n1 | awk '{print $3,$5}' >> coeficientes")
              OPEN(unit=100,file="coeficientes",status='old')
                  READ(100,*)coef2(0),e_coef2(0)  !nmax
                  READ(100,*)coef2(1),e_coef2(1)  !alfa
              CLOSE(100)
              coef2(1)=1.0/coef2(1)
              !CALL SYSTEM ('rm coeficientes fit.log isotermaN.dat')
        CASE ("langmuir") !n = nmax*alfa*P/(1+alfa*P)______________________________________
              ALLOCATE(coef1(0:1),coef2(0:1))
              ALLOCATE(e_coef1(0:1),e_coef2(0:1))
          !-> isoterma1   
              CALL SYSTEM ('cp isoterma1.dat isotermaN.dat')
              CALL SYSTEM ('nohup ./ajustelangmuir.gnu')
              call system('mv ajuste.eps ajuste1.eps')
              call sleep(2)
              CALL SYSTEM ("grep 'Nmax            =' fit.log  | tail -n1 | awk '{print $3,$5}' >  coeficientes")
              CALL SYSTEM ("grep 'alfa            =' fit.log  | tail -n1 | awk '{print $3,$5}' >> coeficientes")
              OPEN(unit=100,file="coeficientes",status='old')
                   READ(100,*)coef1(0),e_coef1(0)  !nmax
                   READ(100,*)coef1(1),e_coef1(1)  !alfa
              CLOSE(100)
              CALL SYSTEM ('rm coeficientes fit.log isotermaN.dat')
          !-> isoterma2
              CALL SYSTEM ('cp isoterma2.dat isotermaN.dat')
              CALL SYSTEM ('nohup ./ajustelangmuir.gnu')
              call system('mv ajuste.eps ajuste2.eps')
              call sleep(2)
              CALL SYSTEM ("grep 'Nmax            =' fit.log  | tail -n1 | awk '{print $3,$5}' >  coeficientes")
              CALL SYSTEM ("grep 'alfa            =' fit.log  | tail -n1 | awk '{print $3,$5}' >> coeficientes")              
              OPEN(unit=100,file="coeficientes",status='old')
                  READ(100,*)coef2(0),e_coef2(0)  !nmax
                  READ(100,*)coef2(1),e_coef2(1)  !alfa
              CLOSE(100)
              CALL SYSTEM ('rm coeficientes fit.log isotermaN.dat')
        CASE ("toth") !n = nmax*alfa*P/(1+(alfa*P)^c)^(1/c)_________________________________
              ALLOCATE(coef1(0:2),coef2(0:2))
              ALLOCATE(e_coef1(0:1),e_coef2(0:1))
          !-> isoterma1
              CALL SYSTEM ('cp isoterma1.dat isotermaN.dat')
              CALL SYSTEM ('nohup ./ajustetoth.gnu')
              call system('mv ajuste.eps ajuste1.eps')
              call sleep(2)
              call system ("grep 'Nmax            =' fit.log  | tail -n1 | awk '{print $3,$5}' >  coeficientes")
              call system ("grep 'alfa            =' fit.log  | tail -n1 | awk '{print $3,$5}' >> coeficientes")
              call system ("grep 'c               =' fit.log  | tail -n1 | awk '{print $3,$5}' >> coeficientes")
              
              OPEN(unit=100,file="coeficientes",status='old')
                   READ(100,*)coef1(0),e_coef1(0)
                   READ(100,*)coef1(1),e_coef1(1)
                   READ(100,*)coef1(2),e_coef1(2)
              CLOSE(100)
             
              CALL SYSTEM ('rm coeficientes isotermaN.dat')
          !-> isoterma2
              CALL SYSTEM ('cp isoterma2.dat isotermaN.dat')
              CALL SYSTEM ('nohup ./ajustetoth.gnu')
              call system('mv ajuste.eps ajuste2.eps')
              call sleep(2)
              call system ("grep 'Nmax            =' fit.log  | tail -n1 | awk '{print $3,$5}' >  coeficientes")
              call system ("grep 'alfa            =' fit.log  | tail -n1 | awk '{print $3,$5}' >> coeficientes")
              call system ("grep 'c               =' fit.log  | tail -n1 | awk '{print $3,$5}' >> coeficientes")
              OPEN(unit=100,file="coeficientes",status='old')
                   READ(100,*)coef2(0),e_coef2(0)  
                   READ(100,*)coef2(1),e_coef2(1)
                   READ(100,*)coef2(2),e_coef2(2)
              CLOSE(100)
              CALL SYSTEM ('rm coeficientes isotermaN.dat')
        CASE ("jensen") !n = K*P/(1+(K*P/(alfa*(1+k*P)))**c)**(1/c)__________________________
              ALLOCATE(coef1(0:3),coef2(0:3))
              ALLOCATE(e_coef1(0:3),e_coef2(0:3))
          !-> isoterma1
              CALL SYSTEM ('cp isoterma1.dat isotermaN.dat')
              CALL SYSTEM ('nohup ./ajustejensen.gnu')
              call system('mv ajuste.eps ajuste1.eps')
              CALL SYSTEM ("grep 'K1              =' fit.log  | tail -n1 | awk '{print $3,$5}' >  coeficientes")
              CALL SYSTEM ("grep 'alfa            =' fit.log  | tail -n1 | awk '{print $3,$5}' >>  coeficientes")
              CALL SYSTEM ("grep 'k2              =' fit.log  | tail -n1 | awk '{print $3,$5}' >>  coeficientes")
              CALL SYSTEM ("grep 'c               =' fit.log  | tail -n1 | awk '{print $3,$5}' >>  coeficientes")
              OPEN(unit=100,file="coeficientes",status='old')
                   READ(100,*)coef1(0)  !K1
                   READ(100,*)coef1(1)  !alfa
                   READ(100,*)coef1(2)  !k2
                   READ(100,*)coef1(3)  !c
              write(6,*)(coef1(i),i=0,3)
              CLOSE(100)
              CALL SYSTEM ('rm coeficientes isotermaN.dat')
          !-> isoterma2
              CALL SYSTEM ('cp isoterma2.dat isotermaN.dat')
              CALL SYSTEM ('nohup ./ajustejensen.gnu')
              call system('mv ajuste.eps ajuste2.eps')
              CALL SYSTEM ("grep 'K1              =' fit.log  | tail -n1 | awk '{print $3,$5}' >  coeficientes")
              CALL SYSTEM ("grep 'alfa            =' fit.log  | tail -n1 | awk '{print $3,$5}' >>  coeficientes")
              CALL SYSTEM ("grep 'k2              =' fit.log  | tail -n1 | awk '{print $3,$5}' >>  coeficientes")
              CALL SYSTEM ("grep 'c               =' fit.log  | tail -n1 | awk '{print $3,$5}' >>  coeficientes")
              OPEN(unit=100,file="coeficientes",status='old')
                  READ(100,*)coef2(0)  !K
                  READ(100,*)coef2(1)  !alfa
                  READ(100,*)coef2(2)  !k
                  READ(100,*)coef2(3)  !c
              write(6,*)(coef1(i),i=0,3)
              CLOSE(100)
              CALL SYSTEM ('rm coeficientes fit.log isotermaN.dat')
 END SELECT
!========================================================================
! CALCULO DE  AREAS
!========================================================================

!vemos el rango de valores de presion de las isotermas para elegir las divisiones
 IF (abs(x1(1)-x1(k))>100)THEN
     divisiones1 = "logaritmica"
 ELSE
     divisiones1 = "linear"
 END IF

 IF (abs(x2(1)-x2(l))>100) THEN
     divisiones2 = "logaritmica"
 ELSE
     divisiones2 = "linear"
 END IF
 ALLOCATE (puntos1(0:intervalos),puntos2(0:intervalos))
 !SELECT CASE (divisiones1)
 ! CASE ("logaritmica")
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
  !CASE ("linear")
  ! IF (x1(1)<x2(1)) THEN
  !  inferior =x1(1)
  ! ELSE
  !  inferior =x2(1)
  ! END IF
  ! IF (x1(k)<x2(l)) THEN
  !  superior =x2(l)
  ! ELSE
  !  superior =x1(k)
  ! END IF
  ! h1 = real((superior-inferior)/intervalos)
  ! DO i=0,intervalos
  !      puntos1(i) = inferior+i*h1
  ! END DO
 !END SELECT
 !SELECT CASE (divisiones2)
  !      CASE ("logaritmica")
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
 !       CASE ("linear")
 !            IF (x1(1)<x2(1)) THEN
 !                  inferior =x1(1)
 !            ELSE
 !                  inferior =x2(1)
 !            END IF
 !            IF (x1(k)<x2(l)) THEN
 !                  superior =x2(l)
 !            ELSE
 !                  superior =x1(k)
 !            END IF
 !            h2 = real((superior-inferior)/intervalos)
 !            DO i=0,intervalos
 !                 puntos2(i) = inferior+i*h2
 !            END DO
 !END SELECT
 ALLOCATE(funcion1(0:intervalos),funcion2(0:intervalos))
 SELECT CASE (ajuste)
        case('freundlich')
         do i=0,intervalos
          funcion1(i) = coef1(0)*puntos1(i)**(1.0/coef1(1))
          funcion2(i) = coef2(0)*puntos2(i)**(1.0/coef2(1))
         end do
        CASE ("langmuir")
             DO i=0,intervalos          
                  funcion1(i) = coef1(0)*coef1(1)*puntos1(i)/(1+coef1(1)*puntos1(i))
                  funcion2(i) = coef2(0)*coef2(1)*puntos2(i)/(1+coef2(1)*puntos2(i))
             END DO
        CASE ("langmuiranalitico")
             !a1 y a2 son el extremo inferior de la integral definida para ambas funciones
             a1 = coef1(0)*(inferior-(log(1+coef1(1)*inferior)/coef1(1)))
             a2 = coef2(0)*(inferior-(log(1+coef2(1)*inferior)/coef2(1)))
             ALLOCATE(PI1(1:intervalos),PI2(1:intervalos))
             DO i=1,intervalos-1
                 ! punto = inferior +i*h
                  ! la integral de langmuir es : nmax*[p-ln(1+alfa*p)/alfa]
                  PI1(i)= coef1(0)*(punto-(log(1+coef1(1)*punto)/coef1(1)))-a1
                  PI2(i)= coef2(0)*(punto-(log(1+coef2(1)*punto)/coef2(1)))-a2
                  funcion1(i) = coef1(0)*coef1(1)*punto/(1+coef1(1)*punto)
                  funcion2(i) = coef2(0)*coef2(1)*punto/(1+coef2(1)*punto)
             END DO
             deallocate(PI1,PI2)
        CASE ("toth")
             DO i=0,intervalos
                  funcion1(i) = coef1(0)*coef1(1)*puntos1(i)/&
                                ((1+(coef1(1)*puntos1(i))**coef1(2))**(1/coef1(2)))
                  funcion2(i) = coef2(0)*coef2(1)*puntos2(i)/&
                                ((1+(coef2(1)*puntos2(i))**coef2(2))**(1/coef2(2)))
             END DO
        CASE ("jensen")
             DO i=0,intervalos
                  funcion1(i) = coef1(0)*puntos1(i)/(1+(coef1(0)*puntos1(i)/&
                               (coef1(1)*(1+coef1(2)*puntos1(i))))**coef1(3))**(1/coef1(3))
                  funcion2(i) = coef2(0)*puntos2(i)/(1+(coef2(0)*puntos2(i)/&
                               (coef2(1)*(1+coef2(2)*puntos2(i))))**coef2(3))**(1/coef2(3))
             END DO
 END SELECT
 OPEN(221,file='iso1.dat')
 OPEN(222,file='iso2.dat')
 do i=0,intervalos
  WRITE(221,*)puntos1(i),funcion1(i)
  WRITE(222,*)puntos2(i),funcion2(i)
 end do
 close(221)
 close(222)
! calculamos el area bajo la curva de ajuste utilizando el metodo de los trapecios
 ALLOCATE (area1(0:intervalos),area2(0:intervalos),PI1(0:intervalos),PI2(0:intervalos))
 DO i=0,intervalos-1
      area1(i)= real((funcion1(i)+funcion1(i+1))*h1/2)
      area2(i)= real((funcion2(i)+funcion2(i+1))*h2/2)
      s1 = s1 + area1(i)
      s2 = s2 + area2(i)
      PI1(i+1)= s1
      PI2(i+1)= s2
 END DO
!========================================================================
! IAST
!========================================================================
 OPEN(unit=104,file='adsorcion.dat')
 DO i=1,intervalos
      DO j=1,intervalos
          comp = abs(PI1(i)-PI2(j))
          IF (comp <= tol) THEN
               IF (abs(x1(1)-x1(k))>100)THEN
                  presion1 = 10**(inferior+i*h1)
               ELSE
                   presion1 = inferior+i+h1
               END IF
               IF (abs(x2(1)-x2(l))>100) THEN
                  presion2 = 10**(inferior+j*h2)
               ELSE
                   presion2 = inferior+j*h2
               END IF
               concx1 = presion2*concy1/(concy2*presion1+concy1*presion2)
               concx2 = 1-concx1
               
               p = presion1*concx1/concy1
               n1 = funcion1(i)*concx1
               n2 = funcion2(j)*concx2
               n = n1+n2
               WRITE(104,*)p,n1,n2,n
          END IF 
     END DO
 END DO
 CLOSE(104)
 call system ("awk '{print $1,$2,$3,$4}' adsorcion.dat | sort -gk1 > c")
 call system ("mv c adsorcion.dat")
 STOP 'Winter is coming'
 CONTAINS
! 
 subroutine get_seed(seed)
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
!  Force 0 < TEMP <= 1.
  DO
    IF(temp<=0.0D+00 )then
       temp=temp+1.0D+00
       CYCLE
    ELSE
       EXIT
    ENDIF
  ENDDO
  DO
    IF(1.0D+00<temp)then
       temp=temp-1.0D+00
       CYCLE
    else
       EXIT
    ENDIF
  ENDDO
  seed=int(dble(i4_huge)*temp)
  if(seed==0)       seed = 1
  if(seed==i4_huge) seed = seed-1
  RETURN
 END subroutine get_seed
!
 REAL function iast_rand(iseed,imode)
!  imode = 1 => random number between 0 and 1
!        = 2 => random number between -1 and 1
  implicit none
  integer  :: iseed
  integer  :: imode
!  Local variables
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
    WRITE(6,*)'array bounds exceeded in GULP_random'
    STOP 'iast_random'
  endif
  if (imode.eq.1) then
    iast_rand = rr(j)
  else
    iast_rand = 2.0*rr(j) - 1.0
  endif
  ix1 = mod(ia1*ix1 + ic1,m1)
  ix2 = mod(ia2*ix2 + ic2,m2)
  rr(j) = (float(ix1) + float(ix2)*rm2)*rm1
  return
  end function iast_rand
END PROGRAM
