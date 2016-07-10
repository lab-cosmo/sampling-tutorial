program mcnvt
  use random
  use routines

  implicit none

  ! * options to be given in input file *
  ! seed              initial seed for the random number generator
  ! temp              target temperature
  ! dataxyz           initial box
  ! nstep             number of steps to be performed
  ! stridetrj         stride outputting the trajectory
  ! stridelog         stride outputting the simulation logs
  ! mcstep            MC step
  ! outputf           output trajectory
  
  double precision temp, mcstep, beta, oldU, newU, pot
  double precision oldposition(3)
  integer argc, seed, nstep, stridetrj, stridelog
  character *256 inputf, dataxyz, outputf
  namelist /inp/ seed, temp, dataxyz, nstep, stridetrj, stridelog, mcstep, outputf

  double precision, allocatable :: q(:,:)

  integer istep, i , natoms
  integer rndindex, attempts, accepted

  ! reads command line
  argc = iargc()
  if ( argc .ne. 1 ) then
    write(6,*) ' Call me as: mcnvt <inputfile> '
    stop
  endif  

  ! reads input file
  call getarg (1,inputf)
  open(101,file=trim(inputf))
  read(101,inp)
  close(unit=101)
  call random_init(seed) ! initialize random number generator
  beta=1.d0/temp  
  
  ! read the box and the number of atoms
  call readdata(dataxyz,natoms,q)
  
  ! open the output files
  open(unit=12,file=trim(outputf),status='replace',action='write')
  ! open(unit=11,file="logs.dat",status='replace',action='write')
  ! write the header
  write(6,*) "## MC (NVT) code. "
  write(6,*) "## Natoms: ", natoms
  write(6,*) "## Temperature: ", temp
  write(6,*) "## MC step: ", mcstep
  write(6,*) "## Step , Pot. Energy "
  
  ! start mc code
  attempts= 0
  accepted= 0
  ! compute the potential energy of the initial configuration
  call get_U(natoms,q,pot)
  do istep=1,nstep
    attempts=attempts+1
    ! save the potential enrgy before the move
    oldU=pot
    ! select a particle randomly
    rndindex = int(dble(natoms)*random_uniform()) + 1
    ! save the original postion
    oldposition=q(:,rndindex)
    ! randomly displace the selected particles
    do i=1,3
      q(i,rndindex)=q(i,rndindex)+ &
                            mcstep*random_gaussian()
    enddo
    ! compute the potential enrgy after the move
    call get_U(natoms,q,newU)
    
    ! acceptance test
    if(random_uniform().lt.exp(-beta*(newU-oldU))) then
      ! accepted
      accepted = accepted+1
      pot=newU
    else
      q(:,rndindex)=oldposition
    endif
    ! write out positions and logs
    if(modulo(istep,stridelog)==0) &
      call write_log(6,istep,pot)
    if(modulo(istep,stridetrj)==0) then
      ! put the origin of the Cartesian coordinate system in the CM
      call adjustqcm(natoms,q)
      call write_xyz(12,istep,natoms,q)
    endif
  enddo
  ! print out the acceptance ratio : 
  write(6,*) "##"
  write(6,*) "## Total moves: ", attempts
  write(6,*) "## Accepteed moves: ", accepted
  write(6,*) "##"
  write(6,*) "## Ratio: ", dble(accepted)/dble(attempts)
  write(6,*) "##"  
  ! stop mc code 
  deallocate(q)
  ! close the files
  close(12)
!  close(11)
end program mcnvt
