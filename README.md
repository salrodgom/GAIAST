# gaiast
=======
This is the IAST code. Updated: Tue 29 Sep, 2015, 14:48:35

Input: 'input' file, isoterma1.dat, isoterma2.dat

$ cat input

jensen        <- isotherm model: langmuir, toth, 

5000

0.5           <- Molar fraction of compound 1 in reservoir

0.5           <- Molar fraction of compound 2 in reservoir

The total molar fraction in the mixture has to be equal to 1

Fitting:

 Genetic Algorithm, 

Models:

$ grep '#model' *.f90

grep '#model' *.f90

   case("freundlich")              !n = a*x**b                                    #model

   case("langmuir")               !n = nmax*alfa*P/(1+alfa*P)                    #model

   case("toth")                   !n=f(x)=Nmax*alfa*x/(1+(alfa*x)**c)**(1/c)     #model

   case("jensen")                 !n = k1*x/(1+(k1*x/(alfa*(1+k2*x))**c))**(1/c) #model

   case("dubinin_raduschkevich")  ! N=Nm*exp(-(RT/Eo ln(Po/P))^2)                #model
 
   case("langmuir_dualsite")      ! N=Nm*b*P/(1+b*P) + Nn*c*P/(1+c*P)            #model
 
   case("dubinin_astakhov")       ! N=Nm*exp(-(RT/Eo ln(Po/P))^d)                #model

Installation, Configuration, Usage:

$ make

make install
make execute
make clean
