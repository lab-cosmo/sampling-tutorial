module routines
  use random
  
  implicit none
  contains
  
  subroutine lwnthermostat(natoms,dt,friction,temp,p)
    ! Langevin white-noise thermostat
    implicit none
    integer, intent(in) :: natoms
    double precision, intent(in) :: dt
    double precision, intent(in) :: friction
    double precision, intent(in) :: temp
    double precision, intent(inout) :: p(3,natoms)
    double precision wnt,wns
    integer i,j
    
    wnt=exp(-friction*dt)
    wns=sqrt((1.0d0-wnt**2.0d0)*temp)
    do i=1,natoms
      do j=1,3
        ! Langevin thermostat contribution      
        p(j,i)=wnt*p(j,i)+wns*random_gaussian()
      enddo
    end do
  end subroutine lwnthermostat
  
  subroutine adjustpcm(natoms,p)
    ! subtract the momentum of the cm
    implicit none
    integer, intent(in) :: natoms
    double precision, intent(inout) :: p(3,natoms)
    double precision pcm(3)
    integer i
    
    pcm=0.0d0
    do i=1,3
      pcm(i)=sum(p(i,:))
    enddo
    pcm=pcm/natoms
    do i=1,natoms
      p(:,i)=p(:,i)-pcm
    enddo
  end subroutine adjustpcm

  subroutine generate_v(natoms,temp,p)
    ! initialize the velocities according to the T
    ! sampling them from a Gaussian distribution
    implicit none
    integer, intent(in)    :: natoms
    double precision, intent(in)    :: temp
    double precision, intent(inout) :: p(3,natoms)
    integer :: i, j
    do i=1,natoms
      do j=1,3
        p(j,i)=dsqrt(temp)*random_gaussian()
      enddo
    end do
  end subroutine generate_v
  
  subroutine get_FandU(natoms,q,f,pot)
    ! compute the LJ forces and potential energy given the positions
    implicit none
    integer, intent(in)  :: natoms
    double precision, intent(in) :: q(3,natoms)
    double precision, intent(out) :: f(3,natoms)
    double precision, intent(out) :: pot
    
    double precision :: dij(3), fij(3), fi(3)
    double precision :: idist2, idist6, idist12
    integer i,j
    f=0.d0
    pot=0.d0

    do i=1,natoms-1
      fi=0.0d0
      ! the neighbourgh list trick should be implemented here
      do j=i+1,natoms
        ! if you want to implement the pbc you should
        ! apply them here
        ! also the cutoff in the force calculation should
        ! be cosidered here

        dij=q(:,i)-q(:,j)
        idist2=1.0d0/dot_product(dij,dij)
        idist6=idist2*idist2*idist2
        idist12=idist6*idist6

        pot=pot+(idist12-idist6)
        fij=dij*(idist12+idist12-idist6)*idist2
        fi = fi+fij
        f(:,j)=f(:,j)-fij        
      enddo
      f(:,i)=f(:,i)+fi
    enddo
    pot=pot*4.0d0
    f=f*24.0d0
    
    ! the above code is a compuntational more efficient way
    ! to do the follow: 
!    f=0.d0
!    pot=0.d0
!    do i=1,natoms-1
!      ! the neighbourgh list trick should be implemented here
!      do j=i+1,natoms
!        ! if you want to implement the pbc you should
!        ! apply them here
!        ! also the cutoff in the force calculation should
!        ! be cosidered here
        
!        dij=q(:,i)-q(:,j)
!        dist2=dot_product(dij,dij)
!        pot=pot+4.0d0*(1.0d0/dist2**6-1.0d0/dist2**3)
!        fij=dij*8.0d0*(6.0d0/dist2**7-3.0d0/dist2**4)
!        f(:,i)=f(:,i)+fij
!        f(:,j)=f(:,j)-fij        
!      enddo
!    enddo
  end subroutine get_FandU
  
  subroutine get_Ek(natoms,p,kin)
    ! compute the kinetic energy
    implicit none
    integer, intent(in) :: natoms
    double precision, intent(in) :: p(3,natoms)
    double precision, intent(out) :: kin
    integer :: i
    
    kin=sum(p*p)*0.5d0
  end subroutine get_Ek
  
  subroutine write_log(fileid,istep,natoms,pot,kin,langham)
    ! write info to the logfile
    implicit none
    integer, intent(in) :: fileid
    integer, intent(in) :: istep
    integer, intent(in) :: natoms
    double precision, intent(in) :: pot
    double precision, intent(in) :: kin
    double precision, intent(in) :: langham
    write(fileid,"(I10,5(ES21.7E3))") istep,2.0*kin/(3.0*natoms-3),kin,pot,pot+kin,langham+pot+kin
  end subroutine write_log
  
  subroutine write_xyz(fileid,istep,natoms,q)
    ! read the positions and the box lenght from the stream fileid
    implicit none
    integer, intent(in) :: fileid
    integer, intent(in) :: istep
    integer, intent(in) :: natoms
    double precision, intent(in) :: q(3,natoms)
    integer i
    
    ! header
    ! wrie the atom numbers
    write(fileid,"(I4)") natoms
    WRITE(fileid,"(A7,I21)") "# Step ", istep
    ! write the atomic coordinates
    DO i=1,natoms
       WRITE(fileid,"(A3,ES21.7E3,A1,ES21.7E3,A1,ES21.7E3)") "Ar " , &
       q(1,i), " ", q(2,i), " ", q(3,i)
    ENDDO
  end subroutine write_xyz
   
  subroutine readdata(dataf,natoms,q)
    ! read the positions and the box lenght from the stream fileid
    implicit none
    character*256 dataf
    integer, intent(out) :: natoms
    double precision, allocatable, dimension(:,:), INTENT(INOUT) :: q
    integer i,endf
    character*256 dummy
    
    open(unit=11,file=dataf)
    ! read the atom numbers
    read(11,*,iostat=endf) natoms 
    ! check the ending of the file or other erros
    if(endf>0) error stop "*** Error occurred while reading file. ***"
    if(endf<0) return
    ! allocate the positions vector if necessesary 
    if (.not.(ALLOCATED(q))) allocate(q(3,natoms))
    !discard a line
    read(11,*,iostat=endf)
    if(endf>0) error stop "*** Error occurred while reading file. ***"
    if(endf<0) return
    !read the atomic coordinates
    do i=1,natoms
       ! indistinguishable particles : discard the atom name
       read(11,*) dummy,q(1,i),q(2,i),q(3,i)
       if(endf>0) error stop "*** Error occurred while reading file. ***"
       if(endf<0) return
    enddo
    close(unit=11) 
  end subroutine readdata

end module routines
