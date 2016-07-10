module routines
  implicit none
  contains

  subroutine get_U(natoms,q,pot)
    ! compute just the potential energy
    implicit none
    integer, intent(in)  :: natoms
    double precision, intent(in)  :: q(3,natoms)
    double precision, intent(out) :: pot
    
    double precision :: dij(3), dist
    integer i,j
    pot=0.d0
    do i=1,natoms-1
      do j=i+1,natoms
        dij=q(:,i)-q(:,j)
        dist=dsqrt(dot_product(dij,dij))
        pot=pot+(1.0d0/dist**12-1.0d0/dist**6)
      enddo
    enddo
    pot=pot*4.0d0
  end subroutine get_U
  
  subroutine adjustqcm(natoms,q)
    ! put the origin in the CM
    implicit none
    integer, intent(in)  :: natoms
    double precision, intent(inout)  :: q(3,natoms)
    double precision :: qcm(3)
    integer i
    qcm=0.0d0
    ! compute the position of CM
    do i=1,natoms
      qcm=qcm+q(:,i)
    enddo
    qcm=qcm/natoms
    ! traslate
    do i=1,natoms
      q(:,i)=q(:,i)-qcm
    enddo   
  end subroutine adjustqcm
  
  subroutine write_log(fileid,istep,pot)
    ! write info to the logfile
    implicit none
    integer, intent(in) :: fileid
    integer, intent(in) :: istep
    double precision, intent(in) :: pot
    write(fileid,"(I10,ES21.7E3)") istep,pot
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
    ! write a blank space
    write(fileid,"(A7,I21)") "# Step ", istep
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
    IF(endf>0) error STOP "*** Error occurred while reading file. ***"
    IF(endf<0) return
    ! allocate the positions vector if necessesary 
    IF (.not.(ALLOCATED(q))) ALLOCATE(q(3,natoms))
    !discard a line
    read(11,*,iostat=endf)
    if(endf>0) error stop "*** Error occurred while reading file. ***"
    if(endf<0) return
    !read the atomic coordinates
    DO i=1,natoms
       ! indistinguishable particles : discard the atom name
       READ(11,*) dummy,q(:,i)
       IF(endf>0) error STOP "*** Error occurred while reading file. ***"
       IF(endf<0) return
    ENDDO
    close(unit=11) 
  end subroutine readdata

end module routines
