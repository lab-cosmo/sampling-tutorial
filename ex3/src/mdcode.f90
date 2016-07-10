program mdcode
  use random
  use routines
  use structure

  implicit none

  ! * options to be given in input file *
  ! seed              initial seed for the random number generator
  ! temp              target temperature
  ! dataxyz           initial box
  ! dt                timestep
  ! nstep             number of steps to be performed
  ! langevinWNtau     Langevin white noise thermostat relaxation time
  !                   (if 0 the WN thermostat won't be applied)
  !                   (if 'NULL' the CN thermostat won't be applied) 
  ! mstep             apply the Langevin thermostat every m step
  ! r0
  ! r1                
  ! sig
  ! stridetrj         stride outputting the trajectory
  ! stridelog         stride outputting the simulation logs
  ! stridestruct      stride outputting the structural number
  ! outputf           output trajectory
  ! structinfof       file that will contain the structure number 
  ! seloutf           file that will contain the selected structures
  ! clist(2)          coordination numbers for which to compute the nc CVs
  ! sela(2)           range of na for which structures are output [sela(1)..sela(2)]
  ! selb(2)           range of nb for which structures are output [selb(1)..selb(2)]
  
  double precision temp, kin, pot, dt, dt2, dt2m
  double precision langevinWNtau, friction, langham
  integer argc, seed, nstep, mstep, stridetrj, stridelog, stridecv
  character*256 inputf, dataxyz, outputf, outcvf, outallf, seloutf

  double precision r0,r1,sig,sig2
  double precision :: nc(2), phi(2), sela(2), selb(2), ncall(10), clist(2)    
  double precision, allocatable :: dnc(:,:,:),cn(:)
  namelist /inp/ seed, temp, dataxyz, dt, nstep, langevinWNtau, mstep, &
                 r0, r1, sig, stridetrj, stridelog, stridecv, &
                 outputf, outcvf, outallf, sela, selb, seloutf, clist

  double precision, allocatable :: q(:,:),p(:,:),f(:,:)


  integer istep, i, natoms, imove, eof

  ! parse the command line
  argc = iargc()
  if ( argc .ne. 1 ) then
    write(6,*) ' Call me as: mdcode <inputfile> '
    stop
  endif  

  ! initialize defaults
  langevinWNtau = 0.0d0
  sig = 0.5d0
  r0=1.5d0
  r1=1.2d0
  mstep=1
  outputf="out.xyz"
  outcvf="out.cv"
  outallf=""
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
  
  ! read the box and the number of atoms
  open(102,file=trim(dataxyz))
  call readdata(102,natoms,q,eof)
  ! positions is allocated inside the routine readdata
  ! now we have to allocate velocities and forces
  
  allocate(p(3,natoms),f(3,natoms))
  allocate(dnc(2,3,natoms))
  allocate(cn(natoms))
  ! open the output files
  open(unit=12,file=trim(outputf),status='replace',action='write')
  open(unit=13,file=trim(outcvf),status='replace',action='write')
  if (trim(outallf)/="") open(unit=15,file=trim(outallf),status='replace',action='write')

  if (sela(1)>=0.0d0 .or. selb(1)>0.0d0) open(unit=14,file=trim(seloutf),status='replace',action='write')
  ! write the header
  write(6,*) "## MD code. "
  write(6,*) "## Natoms: ", natoms
  write(6,*) "## Timestep: ", dt
  write(6,*) "## Temperature: ", temp
  if(langevinWNtau.ne.0) then
    write(6,*) "## Use the Langevin white-noise thermostat with a tau ", langevinWNtau
    write(6,*) "## The thermostat is applied every ", mstep, " step "
  endif
  write(6,*) "##"
  write(6,*) "## Step, Temp, Ekin, Epot, Etot, Conserved quantity"

  write(13,*) " ## Step   na   phia   nb   phib "
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
  ! calculate shorthands of different parameters
  dt2=0.5d0*dt
  dt2m=mstep*dt2
  sig2=sig*sig
  if(langevinWNtau.ne.0) friction=1.0d0/langevinWNtau

  ! start md loop
  do istep=1,nstep
    ! thermostatting
    langham = langham + sum(p*p)*0.5d0
    ! Langevin thermostat is applied before and after the integrator
    if((langevinWNtau.ne.0).and.(modulo(istep,mstep).eq.0)) then
        call lwnthermostat(natoms,dt2m,friction,temp,p)
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
    
    ! conserved quantity
    langham = langham - sum(p*p)*0.5d0
        
    ! calculate the kinetic energy
    call get_Ek(natoms,p,kin)
    
    
    !write(*,*) "=", nc8,nc9
    ! write positions and logs
    if(modulo(istep,stridelog)==0) &
      call write_log(6,istep,natoms,pot,kin,langham)
    if(modulo(istep,stridetrj)==0) &
      call write_xyz(12,istep,natoms,q)
    if(modulo(istep,stridecv)==0) then   

      call get_ncdnc(natoms,q,clist,r0,r1,sig2,cn,nc,dnc)

      phi(1)=dsqrt(sum(dnc(1,:,:)*dnc(1,:,:))*temp/(2*3.1415927))
      phi(2)=dsqrt(sum(dnc(2,:,:)*dnc(2,:,:))*temp/(2*3.1415927))

      write(13,"(I10,4(ES21.7E3))") istep, nc(1), phi(1), nc(2), phi(2)
      if (trim(outallf)/="") then
         call get_ncall(natoms,cn,sig2,ncall)
         write(15, '(10(ES21.7E3))') ncall
      endif
      if (sela(1)>=0.0d0 .or. selb(1)>0.0d0) then
         if ( (sela(1)<0 .or. (sela(1)<nc(1) .and. sela(2)>nc(1))  ) .and. &
              (selb(1)<0 .or. (selb(1)<nc(2) .and. selb(2)>nc(2))  ) ) then
           call write_xyz(14,istep,natoms,q)
         endif
      endif
    endif
  enddo
  ! stop the md loop
  deallocate(q,p,f)
  deallocate(cn,dnc)
  
  ! close the files
  close(12)
  !close(13)
end program mdcode
