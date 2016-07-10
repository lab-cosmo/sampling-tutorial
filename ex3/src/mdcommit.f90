program mdcommit
  use random
  use routines
  use structure

  implicit none

  ! * options to be given in input file *
  ! seed              initial seed for the random number generator
  ! temp              target temperature
  ! dataxyz           configuration files
  ! dt                timestep
  ! nstep             number of steps to be performed for each trajectory
  ! nshoot            number of trajectories to be shot from each point
  ! r0                parameters for the coordination numbers
  ! r1                
  ! sig
  ! cvoutf            file to output final coordinates of the trajectories
  ! sela(:)           range of na which defines product region [sela(1)..sela(2)]
  ! selb(:)           range of nb which defines product region [selb(1)..selb(2)]

  
  double precision temp, kin, pot, dt, dt2, dt2m
  integer argc, seed, nstep, nshoot, nprod
  character*256 inputf, dataxyz, cvoutf
  double precision r0,r1,sig,sig2
  double precision :: nc(2), phi(2), sela(2), selb(2), ncall(10), clist(2)    
  double precision, allocatable :: dnc(:,:,:),cn(:)
  namelist /inp/ seed, temp, dataxyz, dt, nstep, nshoot, &
                 r0, r1, sig, sela, selb, cvoutf, clist

  double precision, allocatable :: q0(:,:),q(:,:),p(:,:),f(:,:)


  integer istep, i, natoms, itraj, ishot, eof

  ! parse the command line
  argc = iargc()
  if ( argc .ne. 1 ) then
    write(6,*) ' Call me as: mdcommit <inputfile> '
    stop
  endif  

  ! initialize defaults
  sig = 0.5d0
  r0=1.5d0
  r1=1.2d0
  nshoot=1
  nstep=1000
  sela=-1.0d0
  selb=-1.0d0
  clist(1)=8
  clist(2)=9
  ! read the input file
  call getarg (1,inputf)
  open(101,file=trim(inputf))
  read(101,inp)
  close(unit=101)
  call random_init(seed) ! initialize random number generator

  ! calculate shorthands of different parameters
  dt2=0.5d0*dt
  sig2=sig*sig
  
  ! read the box and the number of atoms
  open(102,file=trim(dataxyz))
  call readdata(102,natoms,q0,eof)
  ! positions is allocated inside the routine readdata
  ! now we have to allocate velocities and forces
  
  allocate(q(3,natoms),p(3,natoms),f(3,natoms))
  allocate(dnc(2,3,natoms))
  allocate(cn(natoms))

  ! open the output files
  open(unit=12,file=trim(cvoutf),status='replace',action='write')

  ! write the header  
  write(6,*)  " ## Configuration Committor n4, ... n13 "
  write(12,*) " ## Configuration ShotNumber  nA   nB  "
  itraj=0
  do while (eof==0)
    itraj=itraj+1    
    nprod = 0
    do ishot=1,nshoot
        q=q0
        ! generate the initial velocities according to the target T
        call generate_v(natoms,temp,p)
        ! if specified subtract the velocity of the CM
        call adjustpcm(natoms,p)
        ! compute the forces and the potential energy
        call get_FandU(natoms,q,f,pot)
        ! calculate the kinetic energy
        call get_Ek(natoms,p,kin)
        ! start md loop
        do istep=1,nstep
          ! Velocity-Verlet integrator
          !   update velocities
          p=p+f*dt2
          !   update positions
          q=q+p*dt
          !   compute forces (and potential energy)
          call get_FandU(natoms,q,f,pot)
          !   update velocities
          p=p+f*dt2
        enddo   
    
        call get_ncdnc(natoms,q,clist,r0,r1,sig2,cn,nc,dnc)

        write(12,"(2I10,2(ES21.7E3))") itraj, ishot, nc(1), nc(2)

        if (sela(1)>=0.0d0 .or. selb(1)>0.0d0) then
         if ( (sela(1)<0 .or. (sela(1)<nc(1) .and. sela(2)>nc(1))  ) .and. &
              (selb(1)<0 .or. (selb(1)<nc(2) .and. selb(2)>nc(2))  ) ) then
            nprod=nprod+1
         endif   
        endif     
    enddo    
    call get_ncdnc(natoms,q0,clist,r0,r1,sig2,cn,nc,dnc)
    call get_ncall(natoms,cn,sig2,ncall)
    write(6,"(I10,11(ES21.7E3))") itraj, dble(nprod)/dble(nshoot), ncall
    call readdata(102,natoms,q0,eof)
  enddo
  ! stop the md loop
  deallocate(q0,q,p,f)
  deallocate(cn,dnc)
  
  ! close the files
  close(12)
  
end program mdcommit
