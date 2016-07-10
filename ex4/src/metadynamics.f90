module metadynamics
  !use structure
  
  implicit none
  contains
  
  subroutine get_bias(natoms,mxmeta,ihill,hills,whills,sigmeta2,nc,dnc,fbias,bias)

    ! comp nc, its first derivatives and coordination numbers
    implicit none
    integer, intent(in)  :: natoms, mxmeta, ihill
    double precision, intent(in) :: hills(2,mxmeta)
    double precision, intent(in) :: whills(mxmeta), sigmeta2
    double precision, intent(in) :: nc(2)    
    double precision, intent(in) :: dnc(2,3,natoms)
    double precision, intent(out) :: bias
    double precision, intent(out) :: fbias(3,natoms)
    
    double precision :: diffs(2),tmpbias,dbias(2), d2ons2                         
    integer t,k
    
    bias=0.0d0
    fbias=0.0d0
    tmpbias=0.0d0
    dbias=0.0d0
    do t=1,ihill      
      diffs=nc-hills(:,t)
      d2ons2=sum(diffs*diffs)/sigmeta2
      if (d2ons2 .gt. 20.0d0) cycle  ! skip tiny terms
      ! compute the exponential
      tmpbias=dexp(-0.5d0*d2ons2)*whills(t)
      ! cumulate the contributions to dbias/ds
      dbias=dbias+diffs*tmpbias
      ! cumulate the gaussian to get the bias potential
      bias=bias+tmpbias
    end do
    ! mupltiply the weight
    bias=bias
    ! get the forces
    do k=1,2
      fbias=fbias+(dbias(k)*dnc(k,:,:))
    enddo
    fbias=fbias/sigmeta2
  end subroutine get_bias

end module metadynamics
