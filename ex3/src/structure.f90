module structure
  implicit none
  contains
  
  subroutine smoothingfunc(y, cy, dcy)
     implicit none
     double precision, intent(in) :: y
     double precision, intent(inout) :: cy, dcy
     
     if(y.le.0) then
       cy=1
       dcy=0
     elseif(y.ge.1) then 
       cy=0
       dcy=0
     else
       cy=((y-1)**2)*(1+2*y)
       dcy=6*y*(y-1)
     endif
  end subroutine smoothingfunc
  
  subroutine get_ncall(natoms,cn,sig2,ncall)
    implicit none
    integer, intent(in)  :: natoms
    double precision, intent(in) :: sig2,cn(natoms)
    double precision, intent(out) :: ncall(10)
    integer c, i
    ncall = 0
    do i=1,natoms    
      do c = 1,10
        ncall(c)=ncall(c)+dexp(-(cn(i)-(c+3))**2/(2.0d0*sig2))
      enddo
    enddo    
  end subroutine

  subroutine get_ncdnc(natoms,q,c,r0,r1,sig2,cn,nc,dnc)
    ! comp nc, its first derivatives and coordination numbers
    implicit none
    double precision, intent(in) :: c(2)
    integer, intent(in)  :: natoms
    double precision, intent(in) :: q(3,natoms)
    double precision, intent(out) :: cn(natoms)
    double precision, intent(out) :: nc(2)
    double precision, intent(out) :: dnc(2,3,natoms)
    double precision, intent(in) :: r0,r1,sig2
    
    double precision :: dij(3),dist,cij,dcij,ai(2,natoms),bij
    double precision ::  nctmp
    integer i,j,k
    
    nc=0.0d0
    dnc=0.0d0
    cn=0.0d0

    do i=1,natoms
      do j=i+1,natoms
         ! get the distance
         dij=q(:,i)-q(:,j)
         dist=dsqrt(dot_product(dij,dij))
         ! get the coordination number
         call smoothingfunc((dist-r1)/(r0-r1), cij, dcij)
         cn(i)=cn(i)+cij
         cn(j)=cn(j)+cij
      enddo
      do k=1,2
        nctmp=dexp(-(cn(i)-c(k))**2/(2.0d0*sig2))
        ! structure number
        nc(k)=nc(k)+nctmp
        ! first derivative of the structure number
        ai(k,i) = (cn(i)-c(k)) * nctmp*(1.0d0/sig2) 
      enddo
    enddo

    do i=1,natoms
      do j=i+1,natoms
         ! get the distance
         dij=q(:,i)-q(:,j)
         dist=dsqrt(dot_product(dij,dij))
         ! get the coordination number
         call smoothingfunc((dist-r1)/(r0-r1), cij, dcij)
         if (dcij .ne. 0.0d0) then
            bij = -dcij/(r0-r1)/dist
            do k=1,2
              dnc(k,:,i)=dnc(k,:,i)+(ai(k,i)+ai(k,j))*dij*bij
              dnc(k,:,j)=dnc(k,:,j)-(ai(k,i)+ai(k,j))*dij*bij
            enddo
         endif         
      enddo      
    enddo
    ! write(*,*) cnp
    
    
  end subroutine get_ncdnc

end module structure
