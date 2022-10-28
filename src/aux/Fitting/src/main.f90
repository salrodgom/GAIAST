program main
  use, intrinsic                 :: iso_fortran_env
  use pseudorandoms
  use globals
  use GA
  use simplex_
  implicit none
  integer                        :: i
  call get_random_seed( )
  call clean_screen()
  print '(4a)', 'This program was compiled by ', &
       compiler_version(), ' using the options ', &
       compiler_options()
  print '(a)',' '
  call ReadInputFile()
  call Bounds()
  if(Fitting_Flag)then
    do i = 1, Ncomponents
      write(6,'(a)')' '
      write(6,'("Fitting compound:",1x,i2)') i
      write(6,'(a14,1x,a24,1x,a24,10x,a14)')'Compound/Step:','Cromosome:','Parameters','Control'
      call FitGA(i)
      !call WriteFit(i)
      call FitSimplex(i)
    end do
  end if
  call deallocateAll()
  stop ':)'
end program main
