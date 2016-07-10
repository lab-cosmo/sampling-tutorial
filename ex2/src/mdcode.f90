program mdcode
  use random
  use routines
  use glecn

  implicit none

  ! * options to be given in input file *
  ! seed              initial seed for the random number generator
  ! temp              target temperature
  ! dataxyz           initial box
  ! dt                timestep
  ! nstep             number of steps to be performed
  ! langevinWNtau     Langevin white noise thermostat relaxation time
  !                   (if 0 the WN thermostat won't be applied)
  ! gleafile          A matrix GLE CN thermostat 
  !                   (if 'NULL' the CN thermostat won't be applied) 
  ! mstep             apply the Langevin thermostat every m step
  ! andersentau       Andersen thermostat collision frequency
  !                   (if 0 the And. thermostat won't be applied)
  ! gleafile          A matrix GLE CN thermostat 
  ! stridetrj         stride outputting the trajectory
  ! stridelog         stride outputting the simulation logs
  ! outputf           output trajectory
  
  double precision temp, kin, pot, dt, dt2, dt2m
  double precision langevinWNtau, andersentau, andersenprob, friction, langham
  integer argc, seed, nstep, mstep, stridetrj, stridelog
  character*256 inputf, dataxyz, outputf, gleafile
  namelist /inp/ seed, temp, dataxyz, dt, nstep, langevinWNtau, gleafile, & 
                 mstep, andersentau, stridetrj, stridelog, outputf

  double precision, allocatable :: q(:,:),p(:,:),f(:,:)


  integer istep, i, natoms

  ! parse the command line
  argc = iargc()
  if ( argc .ne. 1 ) then
    write(6,*) ' Call me as: mdcode <inputfile> '
    stop
  endif  

  ! initialize defaults
  langevinWNtau = 0.0d0
  andersentau = 0.0d0
  gleafile='NULL'
  mstep=1
  ! read the input file
  call getarg (1,inputf)
  open(101,file=trim(inputf))
  read(101,inp)
  close(unit=101)
  call random_init(seed) ! initialize random number generator
  
  ! read the box and the number of atoms
  call readdata(dataxyz,natoms,q)
  ! positions is allocated inside the routine readdata
  ! now we have to allocate velocities and forces
  
  allocate(p(3,natoms),f(3,natoms))
  

  ! open the output files
  open(unit=12,file=trim(outputf),status='replace',action='write')
  !open(unit=13,file=trim(outputf)//".vxyz",status='replace',action='write')
  ! write the header
  write(6,*) "## MD code. "
  write(6,*) "## Natoms: ", natoms
  write(6,*) "## Timestep: ", dt
  write(6,*) "## Temperature: ", temp
  if(langevinWNtau.ne.0) then
    write(6,*) "## Use the Langevin white-noise thermostat with a tau ", langevinWNtau
    write(6,*) "## The thermostat is applied every ", mstep, " step "
  endif
  if(gleafile.ne."NULL") write(6,*) "## Use the Langevin colored-noise thermostat "
  if(andersentau.ne.0) write(6,*) "## Use the Andersen thermostat "
  write(6,*) "##"
  write(6,*) "## Step, Temp, Ekin, Epot, Etot, Conserved quantity"
  write(6,*) 
  
  ! generate the initial velocities according to the target T
  call generate_v(natoms,temp,p)
  ! if specified subtract the velocity of the CM
  call adjustpcm(natoms,p)
  ! compute the forces and the potential energy
  call get_FandU(natoms,q,f,pot)
  ! calculate the kinetic energy
  call get_Ek(natoms,p,kin)

  ! initialize to zero the accumulator for the 'conserved' quantity
  langham=0.d0    
  call write_log(6,0,natoms,pot,kin,langham)
  istep=0
  call write_xyz(12,istep,natoms,q)
  ! calculate once the follow
  dt2=0.5d0*dt
  dt2m=mstep*dt2
  if(andersentau.ne.0) andersenprob=dt/andersentau
  if(langevinWNtau.ne.0) friction=1.0d0/langevinWNtau
  !initialize the colored-noise thermostat
  if(gleafile.ne.'NULL') call gle_init(natoms,gleafile,temp,dt2m)

  ! start md loop
  do istep=1,nstep
    ! thermostatting
    langham = langham + sum(p*p)*0.5d0
    ! Andersen thermostat
    if((andersentau.ne.0).and.(random_uniform().lt.andersenprob)) then
      call generate_v(natoms,temp,p)
      call adjustpcm(natoms,p)
    endif
    ! Langevin thermostat is applied before and after the integrator
    if((langevinWNtau.ne.0).and.(modulo(istep,mstep).eq.0)) then
        call lwnthermostat(natoms,dt2m,friction,temp,p)
        call adjustpcm(natoms,p)
    endif
    if((gleafile.ne."NULL").and.(modulo(istep,mstep).eq.0))then
      call gle_step(natoms,p)
      call adjustpcm(natoms,p)
    endif
    ! conserved quantity
    langham = langham - sum(p*p)*0.5d0

    ! Velocity-Verlet integrator
    !   update velocities
    p=p+f*dt2
    !   update positions
    q=q+p*dt
    !   compute forces (and potential energy)
    call get_FandU(natoms,q,f,pot)
    !   update velocities
    p=p+f*dt2
    
    langham = langham + sum(p*p)*0.5d0
    
    ! Langevin thermostat
    if((langevinWNtau.ne.0).and.(modulo(istep,mstep).eq.0)) then  
        call lwnthermostat(natoms,dt2m,friction,temp,p)
        call adjustpcm(natoms,p)
    endif
    if((gleafile.ne."NULL").and.(modulo(istep,mstep).eq.0))then
      call gle_step(natoms,p)
      call adjustpcm(natoms,p)
    endif
    
    ! conserved quantity
    langham = langham - sum(p*p)*0.5d0
        
    ! calculate the kinetic energy
    call get_Ek(natoms,p,kin)
    
    ! write positions and logs
    if(modulo(istep,stridelog)==0) &
      call write_log(6,istep,natoms,pot,kin,langham)
    if(modulo(istep,stridetrj)==0) &
      call write_xyz(12,istep,natoms,q)
    !if(modulo(istep,stridetrj)==0) &
    !  call write_xyz(13,istep,natoms,px,py,pz)
  enddo
  ! stop the md loop
  deallocate(q,p,f)
  ! close the files
  close(12)
  !close(13)
end program mdcode
