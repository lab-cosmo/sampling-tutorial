module glecn
  use random
  
  implicit none
  real*8, allocatable, save ::  gS(:,:), gT(:,:), gp(:,:), ngp(:,:)
  integer ns
  
  contains
  
    ! matrix exponential by scale & square.
  ! one can diagonalize with lapack, but it's not worth it, as 
  ! we call this routine only once!      
  subroutine matrix_exp(M, n, j, k, EM)
    integer, intent(in)  :: n, j, k
    double precision, intent(in)   :: M(n,n)
    double precision, intent(out)   :: EM(n,n)
    
    double precision :: tc(j+1), tmp(n,n), SM(n,n)
    integer p, i
    tc(1)=1
    do i=1,j
       tc(i+1)=tc(i)/dble(i)
    enddo
    
    !scale
    SM=M*(1./2.**k)
    EM=0.
    do i=1,n
       EM(i,i)=tc(j+1)
    enddo
    
    !taylor exp of scaled matrix
    do p=j,1,-1
       EM=matmul(SM,EM);
       do i=1,n
          EM(i,i)=EM(i,i)+tc(p)
       enddo
    enddo
    
    !square
    do p=1,k
       EM=matmul(EM,EM)
    enddo
  end subroutine matrix_exp
  
  ! brute-force "stabilized" cholesky decomposition.
  ! in practice, we compute LDL^T decomposition, and force
  ! to zero negative eigenvalues.
  subroutine cholesky(SST, S, n)
    integer, intent(in)  :: n
    double precision, intent(in)   :: SST(n,n)
    double precision, intent(out)   :: S(n,n)
    double precision :: L(n,n), D(n,n) 
    integer i,j,k
    S=0.
    L=0.
    D=0.
    do i=1,n
       L(i,i)=1.0
       do j=1,i-1
          L(i,j)=SST(i,j);
          do k=1,j-1
             L(i,j)=L(i,j)-L(i,k)*L(j,k)*D(k,k)
          enddo
          if (D(j,j).ne. 0.0) then
            L(i,j)=L(i,j)/D(j,j) 
          else
            write(0,*) "Warning: zero eigenvalue in LDL^T decomposition."
            L(i,j)=0.
          end if
       enddo
       D(i,i)=SST(i,i)
       do k=1,i-1
          D(i,i)=D(i,i)-L(i,k)**2*D(k,k)
       end do
    enddo
    do i=1,n
       if ( D(i,i).ge. 0.0d0 ) then
         D(i,i)=sqrt(D(i,i))
       else
         write(0,*) "Warning: negative eigenvalue (",D(i,i),")in LDL^T decomposition."
         D(i,i)=0.0
       end if
    end do
    S=matmul(L,D)
  end subroutine cholesky
  
  ! initialize gle_init
  subroutine gle_init(natoms,gleafile,kt,dt)
    implicit none
    integer, intent(in) :: natoms
    character*256, intent(in) :: gleafile
    ! kt is the temp in reduced unit
    double precision, intent(in):: kt
    double precision, intent(in) :: dt!, wopt
    double precision, allocatable :: gA(:,:), gC(:,:), gr(:)
    integer i, j, k, h, cns, ios

    !reads in matrices
    !reads A (in units of the "optimal" frequency of the fitting)
    open(121,file=gleafile,status='OLD',iostat=ios)
    read(121,*) ns

    !allocate everything we need
    allocate(gA(ns+1,ns+1))
    allocate(gC(ns+1,ns+1))
    allocate(gS(ns+1,ns+1))
    allocate(gT(ns+1,ns+1))
    allocate(gp(3*natoms,ns+1))   
    allocate(ngp(3*natoms,ns+1))   
    allocate(gr(ns+1))
    
    if (ios.ne.0) write(0,*) "Error: could not read GLE-A file!"
    do i=1,ns+1
       read(121,*) gA(i,:)
    enddo
    close(121)

    ! gamma for a WN langevin will be 1/tau, which would make it optimal for w=1/(2tau) angular frequency. 
    ! we scale gA (which is expected to be fitted such that maximum fitted frequency is one) accordingly 
    !gA=gA*wopt
    
    ! init C to kT
    gC=0.0d0
    do i=1,ns+1
      gC(i,i)=kt
    enddo

    ! the deterministic part of the propagator is obtained in a second
    call matrix_exp(-dt*gA, ns+1,15,15,gT)
    
    ! the stochastic part is just as easy. we use gA as a temporary array
    gA=gC-matmul(gT,matmul(gC,transpose(gT)))
    call cholesky(gA, gS, ns+1)

    ! then, we must initialize the auxiliary vectors. we keep general - as we might be 
    ! using non-diagonal C to break detailed balance - and we use cholesky decomposition
    ! of C. again, since one would like to initialize correctly the velocities in 
    ! case of generic C, we use an extra slot for gp for the physical momentum, as we 
    ! could then use it to initialize the momentum in the calling code
    gA=gC   
    call cholesky(gA, gC, ns+1)
    
    do j=1,3*natoms
      do i=1,ns+1
        gr(i)=random_gaussian()
      enddo
      gp(j,:)=matmul(gC,gr)
    end do
    
    deallocate(gA)
    deallocate(gC)
    deallocate(gr)
  end subroutine gle_init

  ! the GLE propagator. 
  ! gT contains the deterministic (friction) part, and gS the stochastic (diffusion) part.
  ! gp(j,1) must be filled with the mass-scaled actual j-th momentum, and contains in
  ! gp(j,2:ns+1) the current values of additional momenta. 
  ! the matrix multiplies are performed on the whole array at once, and the new momentum
  ! is passed back to the caller, while the new s's are kept stored in gp.
  ! please note that one can avoid the double conversion between mass-scaled and actual
  ! momentum/velocity (as used by the calling code) by scaling with the mass the white-noise
  ! random numbers. 
  subroutine gle_step(natoms,p)
    implicit none
    integer, intent(in) :: natoms
    double precision, intent(inout) :: p(3,natoms)
    integer i, j
    double precision mfac, totm
    
    do i=1,natoms       
       gp(i,1)=p(1,i)   !<-- if m!= 1, here a proper scaling must be performed
       gp((natoms+i),1)=p(2,i)
       gp((2*natoms+i),1)=p(3,i)
    enddo

    ngp=transpose(matmul(gT,transpose(gp)))

    !now, must compute random part. 
    !first, fill up gp of random n
    do j=1,3*natoms
      do i=1,ns+1
        gp(j,i)=random_gaussian()     !<-- if m!= 1, alternatively one could perform the scaling here (check also init!)
      end do
    end do

    gp=ngp+transpose(matmul(gS,transpose(gp)))

    do i=1,natoms
      p(1,i)=gp(i,1)     !<-- if m!= 1, here a proper inverse scaling must be performed
      p(2,i)=gp((natoms+i),1)
      p(3,i)=gp((2*natoms+i),1)
    end do
    !write (*,*) sum(gp*gp)/(3*natoms*(ns+1))
  end subroutine gle_step
end module glecn
