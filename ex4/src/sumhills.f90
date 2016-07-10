program sumhills
  implicit none
  
  double precision, allocatable :: hills(:,:), whills(:), shills(:) !,bias(:,:)
  double precision :: nci(2),ncf(2),bins(2)
  integer nhills
  character*1024 :: cmdbuffer   ! String used for reading text lines from files
  integer i,j,dim1,dim2,istep
  integer ccmd, endf, commas(4), par_count  ! stores the index of commas in the parameter string
	double precision vpar(5),tmpnc(2)
  
  istep=0
  nci=0.0d0
  ncf=12.0d0
  bins=0.01d0
  !!!!!!! Command line parser !!!!!!!!!!!!!
  do i = 1, iargc()
    call getarg(i, cmdbuffer)
    if (cmdbuffer == "-nci") then ! starting points
      ccmd = 3
    elseif (cmdbuffer == "-ncf") then ! final points
      ccmd = 4
    elseif (cmdbuffer == "-b") then ! bins
      ccmd = 5
    elseif (cmdbuffer == "-nhills") then ! max hills
      ccmd = 6
    elseif (cmdbuffer == "-h") then ! help flag
      call helpmessage
      call exit(-1)
    else
      if (ccmd == 0) then
        write(*,*) ""
        write(*,*) " Wrong usage. Insert the right parameters!"
        call helpmessage
        call exit(-1)
      elseif (ccmd == 3) then! nci1,nci2
        par_count = 1
        commas(1) = 0
        do while (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
          commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
          read(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) nci(par_count)
          par_count = par_count + 1
        enddo
        read(cmdbuffer(commas(par_count)+1:),*) nci(par_count)
      elseif (ccmd == 4) then! ncf1,ncf2
        par_count = 1
        commas(1) = 0
        do while (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
          commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
          read(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) ncf(par_count)
          par_count = par_count + 1
        enddo
        read(cmdbuffer(commas(par_count)+1:),*) ncf(par_count)
      elseif (ccmd == 5) then! b1,b2
        par_count = 1
        commas(1) = 0
        do while (index(cmdbuffer(commas(par_count)+1:), ',') > 0)
          commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
          read(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) bins(par_count)
          par_count = par_count + 1
        enddo
        read(cmdbuffer(commas(par_count)+1:),*) bins(par_count)
      elseif (ccmd == 6) then ! istep
        read(cmdbuffer,*) istep
      endif
    endif
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  dim1=int((ncf(1)-nci(1))/bins(1))+1
  dim2=int((ncf(2)-nci(2))/bins(2))+1
  
  call readinput(nhills,hills,whills,shills)
  shills = shills*shills ! getbias uses sigma^2
  if(istep.eq.0) istep=nhills

  !bias=0.0d0
  do i=1,dim1
    do j=1,dim2
      tmpnc(1)=nci(1)+bins(1)*(i-1)
      tmpnc(2)=nci(2)+bins(2)*(j-1)
      write(6,"(3(ES21.7E3))") tmpnc(1),tmpnc(2), &
               get_bias(tmpnc,istep,nhills,hills,whills,shills)
    enddo
    ! format the output for gnuplot
    write(6,*)
  enddo
  
  deallocate(hills, whills, shills)!,bias)

  contains

    subroutine helpmessage
      ! Banner to print out for helping purpose
    
      write(*,*) ""
      write(*,*) " USAGE: sumhills [-h] [-nci nci1,nci2] "
      write(*,*) "                  [-ncf ncf1,ncf2] [-b w1,w2] [-nhills step]"
      write(*,*) ""
      write(*,*) " This program require two columns in intput from the stdin! "
      write(*,*) ""
      write(*,*) " Example: "
      write(*,*) ""
      write(*,*) "   awk '{print $2,$3,$4,$5}' out.hills | sumhills \ "
      write(*,*) "       -nci 0,0 -ncf 16,16 -b 0.1,0.1 -nhills 1000 > bias-t1000.splot"
      write(*,*) ""
    end subroutine helpmessage
         
    double precision function get_bias(nc,istep,nhills,hills,whills,shills)
      ! comp nc, its first derivatives and coordination numbers
      implicit none
      integer, intent(in) :: istep
      integer, intent(in) :: nhills
      double precision, intent(in) :: hills(2,nhills), whills(nhills), shills(nhills)
      double precision, intent(in) :: nc(2) 
      
      double precision :: diffs(2),d2ons2                         
      integer t
      
      get_bias=0.0d0
      do t=1,istep      
        diffs=nc-hills(:,t)
        d2ons2=sum(diffs*diffs)/shills(t)
        if (d2ons2.gt.20.0d0) cycle  ! skip tiny terms
        get_bias=get_bias+dexp(-0.5d0*d2ons2) * whills(t)
      end do
    end function get_bias
    
    subroutine readinput(nhills,vout, wout, sout)
      integer, intent(out) :: nhills
      double precision, allocatable, intent(out) :: vout(:,:), wout(:), sout(:)
      
      integer, parameter :: nbuff = 1000
      double precision :: vbuff(2,nbuff), wbuff(nbuff), sbuff(nbuff)
      integer io_status,i
      
      nhills=0
      i=0
      allocate(vout(2,nbuff), wout(nbuff), sout(nbuff))
      vout=0.0d0
      do
        i=i+1
        read(5,*, iostat=io_status) vbuff(:,i), wbuff(i), sbuff(i)

        if(io_status<0) exit
        if(io_status>0) error stop "*** Error occurred while reading file. ***"
        if(i.EQ.nbuff) then
          call appendbuffer2(nhills,i,vout,vbuff(:,1:i))
          call appendbuffer(nhills,i,sout,sbuff(1:i))
          call appendbuffer(nhills,i,wout,wbuff(1:i))
  
          nhills = nhills+nbuff
          i=0
        endif
      enddo
      
      if(i>0) then
         call appendbuffer2(nhills,i,vout,vbuff(:,1:i))
         call appendbuffer(nhills,i,sout,sbuff(1:i))
         call appendbuffer(nhills,i,wout,wbuff(1:i))
      endif  
    end subroutine readinput
     
    subroutine appendbuffer2(n1,nbuff,v1,vbuff)
      ! ndimensional version
      ! collapse v1 and v2 into v1
      !
      ! Args:
      !    n1: size of the vector v1.
      !        n1 will be changed whend the routine will finish.
      !    nbuff: size of the vector vbuff
      !    v1: main vector
      !    vbuff: vector to append to v1
      integer, intent(inout) :: n1
      integer, intent(in) :: nbuff
      double precision, allocatable, dimension(:,:), INTENT(INOUT) :: v1
      double precision, dimension(:,:), INTENT(IN) :: vbuff
      
      double precision, allocatable, dimension(:,:) :: tmp_v
      
      if(n1.EQ.0)then
         deallocate(v1)
         allocate(v1(2,nbuff))
         v1=vbuff
      else
         allocate(tmp_v(2,n1+nbuff))
         tmp_v(:,1:n1)=v1
         tmp_v(:,n1+1:n1+nbuff)=vbuff
         deallocate(v1)
         allocate(v1(2,n1+nbuff))
         v1=tmp_v
         deallocate(tmp_v)
      endif
    end subroutine appendbuffer2

    subroutine appendbuffer(n1,nbuff,v1,vbuff)
      ! ndimensional version
      ! collapse v1 and v2 into v1
      !
      ! Args:
      !    n1: size of the vector v1.
      !        n1 will be changed whend the routine will finish.
      !    nbuff: size of the vector vbuff
      !    v1: main vector
      !    vbuff: vector to append to v1
      integer, intent(inout) :: n1
      integer, intent(in) :: nbuff
      double precision, allocatable, dimension(:), INTENT(INOUT) :: v1
      double precision, dimension(:), INTENT(IN) :: vbuff
      
      double precision, allocatable, dimension(:) :: tmp_v

      if(n1.EQ.0)then
         deallocate(v1)
         allocate(v1(nbuff))
         v1=vbuff
      else         
         allocate(tmp_v(n1+nbuff))         
         tmp_v(1:n1)=v1
         tmp_v(n1+1:n1+nbuff)=vbuff
         deallocate(v1)
         allocate(v1(n1+nbuff))
         v1=tmp_v
         deallocate(tmp_v)
      endif
    end subroutine 
    
end program sumhills
