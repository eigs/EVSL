!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!    MATRIX GENERATION ROUTINES  -- FINITE DIFFERENCE MATRICES         c
!----------------------------------------------------------------------c
! contents:                                                            c
!----------                                                            c
! gen57pt  : generates 5-point and 7-point matrices.                   c
! gen57bl  : generates block 5-point and 7-point matrices.             c
!                                                                      c
! supporting routines:                                                 c
!---------                                                             c
! gensten  : generate the stencil (point version)                      c
! bsten    : generate the stencil (block version)                      c
! fdaddbc  : finite difference add boundary conditions                 c
! fdreduce : reduce the system to eliminate node with known values     c
! clrow    : clear a row of a CSR matrix                               c
! lctcsr   : locate the position of A(i,j) in CSR format               c
!----------------------------------------------------------------------c
      subroutine gen57pt(nx,ny,nz,al,mode,n,a,ja,ia,iau,rhs)
      integer ja(*),ia(*),iau(*), nx, ny, nz, mode, n
      real*8 a(*), rhs(*), al(6)
!-----------------------------------------------------------------------
! On entry:
!
! nx      = number of grid points in x direction
! ny      = number of grid points in y direction
! nz      = number of grid points in z direction
! al      = array of size 6, carries the coefficient alpha of the
!           boundary conditions
! mode    = what to generate:
!           < 0 : generate the graph only,
!           = 0 : generate the matrix,
!           > 0 : generate the matrix and the right-hand side.
!
! On exit:
!
! n       = number of nodes with unknown values, ie number of rows
!           in the matrix
!
! a,ja,ia = resulting matrix in row-sparse format
!
! iau     = integer*n, containing the poisition of the diagonal element
!           in the a, ja, ia structure
!
! rhs     = the right-hand side
!
! External functions needed (must be supplied by caller)
!     afun, bfun, cfun, dfun, efun, ffun, gfun, hfun
!     betfun, gamfun
! They have the following prototype:
!     real*8 function xfun(x, y, z)
!     real*8 x, y, z
!-----------------------------------------------------------------------
! This subroutine computes the sparse matrix in compressed sparse row
! format for the elliptic equation:
!       d    du    d    du    d    du      du     du     du
! L u = --(A --) + --(B --) + --(C --) + D -- + E -- + F -- + G u = H u
!       dx   dx    dy   dy    dz   dz      dx     dy     dz
!
! with general Mixed Boundary conditions, on a rectangular 1-D,
! 2-D or 3-D grid using 2nd order centered difference schemes.
!
! The functions a, b, ..., g, h are known through the
! as afun, bfun, ..., gfun, hfun in this subroutine.
! NOTE: To obtain the correct matrix, any function that is not
! needed should be set to zero.  For example for two-dimensional
! problems, nz should be set to 1 and the functions cfun and ffun
! should be zero functions.
!
! The Boundary condition is specified in the following form:
!           du
!     alpha -- + beta u = gamma
!           dn
! Where alpha is constant at each side of the boundary surfaces.  Alpha
! is represented by parameter al.  It is expected to an array that
! contains enough elements to specify the boundaries for the problem,
! 1-D case needs two elements, 2-D needs 4 and 3-D needs 6.  The order
! of the boundaries in the array is left(west), right(east),
! bottom(south), top(north), front, rear.  Beta and gamma are functions
! of type real with three arguments x, y, z.  These two functions are
! known subroutine 'addbc' as betfun and gamfun.  They should following
! the same notion as afun ... hfun.  For more restriction on afun ...
! hfun, please read the documentation follows the subroutine 'getsten',
! and, for more on betfun and gamfun, please refer to the documentation
! under subroutine 'fdaddbc'.
!
! The nodes are ordered using natural ordering, first x direction, then
! y, then z.  The mesh size h is uniform and determined by grid points
! in the x-direction.
!
! The domain specified for the problem is [0 .ge. x .ge. 1],
! [0 .ge. y .ge. (ny-1)*h] and [0 .ge. z .ge. (nz-1)*h], where h is
! 1 / (nx-1).  Thus if non-Dirichlet boundary condition is specified,
! the mesh will have nx points along the x direction, ny along y and
! nz along z.  For 1-D case, both y and z value are assumed to zero
! when calling relavent functions that have three parameters.
! Similarly, for 2-D case, z is assumed to be zero.
!
! About the expectation of nx, ny and nz:
! nx is required to be .gt. 1 always;
! if the second dimension is present in the problem, then ny should be
! .gt. 1, else 1;
! if the third dimension is present in the problem, nz .gt. 1, else 1.
! when ny is 1, nz must be 1.
!-----------------------------------------------------------------------
!
!     stencil [1:7] has the following meaning:
!
!     center point = stencil(1)
!     west point = stencil(2)
!     east point = stencil(3)
!     south point = stencil(4)
!     north point = stencil(5)
!     front point = stencil(6)
!     back point = stencil(7)
!
!     al[1:6] carry the coefficient alpha in the similar order
!
!     west  side = al(1)
!     east  side = al(2)
!     south side = al(3)
!     north side = al(4)
!     front side = al(5)
!     back  side = al(6)
!
!                           al(4)
!                           st(5)
!                            |
!                            |
!                            |           al(6)
!                            |          .st(7)
!                            |     .
!         al(1)              | .             al(2)
!         st(2) ----------- st(1) ---------- st(3)
!                       .    |
!                   .        |
!               .            |
!            st(6)           |
!            al(5)           |
!                            |
!                           st(4)
!                           al(3)
!
!-------------------------------------------------------------------
!     some constants
!
      real*8 one
      parameter (one=1.0D0)
!
!     local variables
!
      integer ix, iy, iz, kx, ky, kz, node, iedge
      real*8  r, h, stencil(7)
      logical value, genrhs
!
!     nx has to be larger than 1
!
      if (nx.le.1) return
      h = one / dble(nx-1)
!
!     the mode
!
      value = (mode.ge.0)
      genrhs = (mode.gt.0)
!
!     first generate the whole matrix as if the boundary condition does
!     not exist
!
      kx = 1
      ky = nx
      kz = nx*ny
      iedge = 1
      node = 1
      do 100 iz = 1,nz
         do 90 iy = 1,ny
            do 80 ix = 1,nx
               ia(node) = iedge
!
!     compute the stencil at the current node
!
               if (value) then
                  call getsten(nx,ny,nz,mode,ix-1,iy-1,iz-1,stencil,h,r)
               end if
!     west
               if (ix.gt.1) then
                  ja(iedge)=node-kx
                  if (value) a(iedge) = stencil(2)
                  iedge=iedge + 1
               end if
!     south
               if (iy.gt.1) then
                  ja(iedge)=node-ky
                  if (value) a(iedge) = stencil(4)
                  iedge=iedge + 1
               end if
!     front plane
               if (iz.gt.1) then
                  ja(iedge)=node-kz
                  if (value) a(iedge) = stencil(6)
                  iedge=iedge + 1
               endif
!     center node
               ja(iedge) = node
               iau(node) = iedge
               if (value) a(iedge) = stencil(1)
               iedge = iedge + 1
!     east
               if (ix.lt.nx) then
                  ja(iedge)=node+kx
                  if (value) a(iedge) = stencil(3)
                  iedge=iedge + 1
               end if
!     north
               if (iy.lt.ny) then
                  ja(iedge)=node+ky
                  if (value) a(iedge) = stencil(5)
                  iedge=iedge + 1
               end if
!     back plane
               if (iz.lt.nz) then
                  ja(iedge)=node+kz
                  if (value) a(iedge) = stencil(7)
                  iedge=iedge + 1
               end if
!     the right-hand side
               if (genrhs) rhs(node) = r
               node=node+1
 80         continue
 90      continue
 100  continue
      ia(node)=iedge
!
!     Add in the boundary conditions
!
      call fdaddbc(nx,ny,nz,a,ja,ia,iau,rhs,al,h)
!
!     eliminate the boudary nodes from the matrix
!
      call fdreduce(nx,ny,nz,al,n,a,ja,ia,iau,rhs,stencil)
!
!     done
!
      return
!-----end-of-gen57pt----------------------------------------------------
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine getsten (nx,ny,nz,mode,kx,ky,kz,stencil,h,rhs)
      integer nx,ny,nz,mode,kx,ky,kz
      real*8 stencil(*),h,rhs,afun,bfun,cfun,dfun,efun,ffun,gfun,hfun
      external afun,bfun,cfun,dfun,efun,ffun,gfun,hfun
!-----------------------------------------------------------------------
!     This subroutine calculates the correct stencil values for
!     centered difference discretization of the elliptic operator
!     and the right-hand side
!
! L u = delx( A delx u ) + dely ( B dely u) + delz ( C delz u ) +
!       delx ( D u ) + dely (E u) + delz( F u ) + G u = H
!
!   For 2-D problems the discretization formula that is used is:
!
! h**2 * Lu == A(i+1/2,j)*{u(i+1,j) - u(i,j)} +
!              A(i-1/2,j)*{u(i-1,j) - u(i,j)} +
!              B(i,j+1/2)*{u(i,j+1) - u(i,j)} +
!              B(i,j-1/2)*{u(i,j-1) - u(i,j)} +
!              (h/2)*D(i,j)*{u(i+1,j) - u(i-1,j)} +
!              (h/2)*E(i,j)*{u(i,j+1) - u(i,j-1)} +
!              (h/2)*E(i,j)*{u(i,j+1) - u(i,j-1)} +
!              (h**2)*G(i,j)*u(i,j)
!-----------------------------------------------------------------------
!     some constants
!
      real*8 zero, half
      parameter (zero=0.0D0,half=0.5D0)
!
!     local variables
!
      integer k
      real*8 hhalf,cntr, x, y, z, coeff
!
!     if mode < 0, we shouldn't have come here
!
      if (mode .lt. 0) return
!
      do 200 k=1,7
         stencil(k) = zero
 200  continue
!
      hhalf = h*half
      x = h*dble(kx)
      y = h*dble(ky)
      z = h*dble(kz)
      cntr = zero
!     differentiation wrt x:
      coeff = afun(x+hhalf,y,z)
      stencil(3) = stencil(3) + coeff
      cntr = cntr + coeff
!
      coeff = afun(x-hhalf,y,z)
      stencil(2) = stencil(2) + coeff
      cntr = cntr + coeff
!
      coeff = dfun(x,y,z)*hhalf
      stencil(3) = stencil(3) + coeff
      stencil(2) = stencil(2) - coeff
      if (ny .le. 1) goto 99
!
!     differentiation wrt y:
!
      coeff = bfun(x,y+hhalf,z)
      stencil(5) = stencil(5) + coeff
      cntr = cntr + coeff
!
      coeff = bfun(x,y-hhalf,z)
      stencil(4) = stencil(4) + coeff
      cntr = cntr + coeff
!
      coeff = efun(x,y,z)*hhalf
      stencil(5) = stencil(5) + coeff
      stencil(4) = stencil(4) - coeff
      if (nz .le. 1) goto 99
!
! differentiation wrt z:
!
      coeff = cfun(x,y,z+hhalf)
      stencil(7) = stencil(7) + coeff
      cntr = cntr + coeff
!
      coeff = cfun(x,y,z-hhalf)
      stencil(6) = stencil(6) + coeff
      cntr = cntr + coeff
!
      coeff = ffun(x,y,z)*hhalf
      stencil(7) = stencil(7) + coeff
      stencil(6) = stencil(6) - coeff
!
! contribution from function G:
!
 99   coeff = gfun(x,y,z)
      stencil(1) = h*h*coeff - cntr
!
!     the right-hand side
!
      if (mode .gt. 0) rhs = h*h*hfun(x,y,z)
!
      return
!------end-of-getsten---------------------------------------------------
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine gen57bl (nx,ny,nz,nfree,na,n,a,ja,ia,iau,stencil)
!     implicit real*8 (a-h,o-z)
      integer ja(*),ia(*),iau(*),nx,ny,nz,nfree,na,n
      real*8 a(na,1), stencil(7,1)
!--------------------------------------------------------------------
! This subroutine computes the sparse matrix in compressed
! format for the elliptic operator
!
! L u = delx( a . delx u ) + dely ( b . dely u) + delz ( c . delz u ) +
!       delx ( d . u ) + dely (e . u) + delz( f . u ) + g . u
!
! Here u is a vector of nfree componebts and each of the functions
! a, b, c, d, e, f, g   is an (nfree x nfree) matrix depending of
! the coordinate (x,y,z).
! with Dirichlet Boundary conditions, on a rectangular 1-D,
! 2-D or 3-D grid using centered difference schemes.
!
! The functions a, b, ..., g are known through the
! subroutines  afunbl, bfunbl, ..., gfunbl. (user supplied) .
!
! uses natural ordering, first x direction, then y, then z
! mesh size h is uniform and determined by grid points
! in the x-direction.
! 
! The output matrix is in Block -- Sparse Row format. 
!
!--------------------------------------------------------------------
! parameters:
!-------------
! Input:
! ------
! nx      = number of points in x direction
! ny      = number of points in y direction
! nz      = number of points in z direction
! nfree   = number of degrees of freedom per point
! na      = first dimension of array a as declared in calling
!           program. Must be .ge. nfree**2
!
! Output: 
! ------ 
! n       = dimension of matrix (output)
!
! a, ja, ia = resulting matrix in  Block Sparse Row format
!           a(1:nfree**2, j ) contains a nonzero block and ja(j) 
!           contains the (block) column number of this block.
!           the block dimension of the matrix is n (output) and 
!           therefore the total number of (scalar) rows is n x nfree.
!     
! iau     = integer*n containing the position of the diagonal element
!           in the a, ja, ia structure
!
! Work space:
!------------ 
! stencil =  work array of size (7,nfree**2) [stores local stencils]
!
!--------------------------------------------------------------------
!
!     stencil (1:7,*) has the following meaning:
!
!     center point = stencil(1)
!     west point   = stencil(2)
!     east point   = stencil(3)
!     south point  = stencil(4)
!     north point  = stencil(5)
!     front point  = stencil(6)
!     back point   = stencil(7)
!
!
!                           st(5)
!                            |
!                            |
!                            |
!                            |          .st(7)
!                            |     .
!                            | .
!         st(2) ----------- st(1) ---------- st(3)
!                       .    |
!                   .        |
!               .            |
!            st(6)           |
!                            |
!                            |
!                           st(4)
!
!-------------------------------------------------------------------
!     some constants
!
      real*8 one
      parameter (one=1.0D0)
!
!     local variables
!
      integer iedge,ix,iy,iz,k,kx,ky,kz,nfree2,node
      real*8  h
!
      h = one/dble(nx+1)
      kx = 1
      ky = nx
      kz = nx*ny
      nfree2 = nfree*nfree
      iedge = 1
      node = 1
      do 100 iz = 1,nz
         do 90 iy = 1,ny
            do 80 ix = 1,nx
               ia(node) = iedge
               call bsten(nx,ny,nz,ix,iy,iz,nfree,stencil,h)
!     west
               if (ix.gt.1) then
                  ja(iedge)=node-kx
                  do 4 k=1,nfree2
                  a(k,iedge) = stencil(2,k)
 4                continue
                  iedge=iedge + 1
               end if
!     south
               if (iy.gt.1) then
                  ja(iedge)=node-ky
                  do 5 k=1,nfree2
                  a(k,iedge) = stencil(4,k)
 5                continue
                  iedge=iedge + 1
               end if
!     front plane
               if (iz.gt.1) then
                  ja(iedge)=node-kz
                  do 6 k=1,nfree2
                  a(k,iedge) = stencil(6,k)
 6                continue
                  iedge=iedge + 1
               endif
!     center node
               ja(iedge) = node
               iau(node) = iedge
               do 7 k=1,nfree2
                  a(k,iedge) = stencil(1,k)
 7             continue
               iedge = iedge + 1
!     -- upper part
!     east
               if (ix.lt.nx) then
                  ja(iedge)=node+kx
                  do 8 k=1,nfree2
                  a(k,iedge) = stencil(3,k)
 8                continue
                  iedge=iedge + 1
               end if
!     north
               if (iy.lt.ny) then
                  ja(iedge)=node+ky
                  do 9 k=1,nfree2
                  a(k,iedge) = stencil(5,k)
 9                continue
                  iedge=iedge + 1
               end if
!     back plane
               if (iz.lt.nz) then
                  ja(iedge)=node+kz
                  do 10 k=1,nfree2
                     a(k,iedge) = stencil(7,k)
 10               continue
                  iedge=iedge + 1
               end if
!------next node -------------------------
               node=node+1
 80         continue
 90      continue
 100  continue
!     
!     -- new version of BSR -- renumbering removed. 
!     change numbering of nodes so that each ja(k) will contain the
!     actual column number in the original matrix of entry (1,1) of each
!     block (k).
!      do 101 k=1,iedge-1
!         ja(k) = (ja(k)-1)*nfree+1
! 101  continue
!
!      n = (node-1)*nfree
      n = node-1 
      ia(node)=iedge
      return
!--------------end-of-gen57bl-------------------------------------------
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
      subroutine bsten (nx,ny,nz,kx,ky,kz,nfree,stencil,h)
!-----------------------------------------------------------------------
!     This subroutine calcultes the correct block-stencil values for
!     centered difference discretization of the elliptic operator
!     (block version of stencil)
!
! L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) +
!       d delx ( u ) + e dely (u) + f delz( u ) + g u
!
!   For 2-D problems the discretization formula that is used is:
!
! h**2 * Lu == a(i+1/2,j)*{u(i+1,j) - u(i,j)} +
!              a(i-1/2,j)*{u(i-1,j) - u(i,j)} +
!              b(i,j+1/2)*{u(i,j+1) - u(i,j)} +
!              b(i,j-1/2)*{u(i,j-1) - u(i,j)} +
!              (h/2)*d(i,j)*{u(i+1,j) - u(i-1,j)} +
!              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} +
!              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} +
!              (h**2)*g(i,j)*u(i,j)
!-----------------------------------------------------------------------
!     some constants
!
      real*8  zero,half
      parameter(zero=0.0D0,half=0.5D0)
!
!     local variables
!
      integer i,k,kx,ky,kz,nfree,nfree2,nx,ny,nz
      real*8 stencil(7,*), cntr(225), coeff(225),h,h2,hhalf,x,y,z
!------------
      if (nfree .gt. 15) then
         print *, ' ERROR ** nfree too large '
         stop
      endif
!
      nfree2 = nfree*nfree
      do 200 k=1, nfree2
         cntr(k) = zero
         do 199 i=1,7
            stencil(i,k) = zero
 199     continue
 200  continue
!------------
      hhalf = h*half
      h2 = h*h
      x = h*dble(kx)
      y = h*dble(ky)
      z = h*dble(kz)
! differentiation wrt x:
      call afunbl(nfree,x+hhalf,y,z,coeff)
      do 1 k=1, nfree2
      stencil(3,k) = stencil(3,k) + coeff(k)
      cntr(k) = cntr(k) + coeff(k)
 1    continue
!
      call afunbl(nfree,x-hhalf,y,z,coeff)
      do 2 k=1, nfree2
         stencil(2,k) = stencil(2,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 2    continue
!
      call dfunbl(nfree,x,y,z,coeff)
      do 3 k=1, nfree2
         stencil(3,k) = stencil(3,k) + coeff(k)*hhalf
         stencil(2,k) = stencil(2,k) - coeff(k)*hhalf
 3    continue
      if (ny .le. 1) goto 99
!
! differentiation wrt y:
!
      call bfunbl(nfree,x,y+hhalf,z,coeff)
      do 4 k=1,nfree2
         stencil(5,k) = stencil(5,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 4    continue
!
      call bfunbl(nfree,x,y-hhalf,z,coeff)
      do 5 k=1, nfree2
         stencil(4,k) = stencil(4,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 5    continue
!
      call efunbl(nfree,x,y,z,coeff)
      do 6 k=1, nfree2
         stencil(5,k) = stencil(5,k) + coeff(k)*hhalf
         stencil(4,k) = stencil(4,k) - coeff(k)*hhalf
 6    continue
      if (nz .le. 1) goto 99
!
! differentiation wrt z:
!
      call cfunbl(nfree,x,y,z+hhalf,coeff)
      do 7 k=1, nfree2
         stencil(7,k) = stencil(7,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 7    continue
!
      call cfunbl(nfree,x,y,z-hhalf,coeff)
      do 8 k=1, nfree2
         stencil(6,k) = stencil(6,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 8    continue
!
      call ffunbl(nfree,x,y,z,coeff)
      do 9 k=1, nfree2
         stencil(7,k) = stencil(7,k) + coeff(k)*hhalf
         stencil(6,k) = stencil(6,k) - coeff(k)*hhalf
 9    continue
!
! discretization of  product by g:
!
 99   call gfunbl(nfree,x,y,z,coeff)
      do 10 k=1, nfree2
         stencil(1,k) = h2*coeff(k) - cntr(k)
 10   continue
!
      return
!------------end of bsten-----------------------------------------------
!-----------------------------------------------------------------------
      end
      subroutine fdreduce(nx,ny,nz,alpha,n,a,ja,ia,iau,rhs,stencil)
      implicit none
      integer nx,ny, nz, n, ia(*), ja(*), iau(*)
      real*8  alpha(*), a(*), rhs(*), stencil(*)
!-----------------------------------------------------------------------
! This subroutine tries to reduce the size of the matrix by looking
! for Dirichlet boundary conditions at each surface and solve the boundary
! value and modify the right-hand side of related nodes, then clapse all
! the boundary nodes.
!-----------------------------------------------------------------------
!     parameters
!
      real*8   zero
      parameter(zero=0.0D0)
!
!     local variables
!
      integer  i,j,k,kx,ky,kz,lx,ux,ly,uy,lz,uz,node,nbnode,lk,ld,iedge
      real*8   val
      integer  lctcsr
      external lctcsr
!
!     The first half of this subroutine will try to change the right-hand
!     side of all the nodes that has a neighbor with Dirichlet boundary
!     condition, since in this case the value of the boundary point is
!     known.
!     Then in the second half, we will try to eliminate the boundary
!     points with known values (with Dirichlet boundary condition).
!
      kx = 1
      ky = nx
      kz = nx*ny
      lx = 1
      ux = nx
      ly = 1
      uy = ny
      lz = 1
      uz = nz
!
!     Here goes the first part. ----------------------------------------
!
!     the left (west) side
!
      if (alpha(1) .eq. zero) then
         lx = 2
         do 10 k = 1, nz
            do 11 j = 1, ny
               node = (k-1)*kz + (j-1)*ky + 1
               nbnode = node + kx
               lk = lctcsr(nbnode, node, ja, ia)
               ld = iau(node)
               val = rhs(node)/a(ld)
!     modify the rhs
               rhs(nbnode) = rhs(nbnode) - a(lk)*val
 11         continue
 10      continue
      endif
!
!     right (east) side
!
      if (alpha(2) .eq. zero) then
         ux = nx - 1
         do 20 k = 1, nz
            do 21 j = 1, ny
               node = (k-1)*kz + (j-1)*ky + nx
               nbnode = node - kx
               lk = lctcsr(nbnode, node, ja, ia)
               ld = iau(node)
               val = rhs(node)/a(ld)
!     modify the rhs
               rhs(nbnode) = rhs(nbnode) - a(lk)*val
 21         continue
 20      continue
      endif
!
!     if it's only 1-D, skip the following part
!
      if (ny .le. 1) goto 100
!
!     the bottom (south) side
!
      if (alpha(3) .eq. zero) then
         ly = 2
         do 30 k = 1, nz
            do 31 i = lx, ux
               node = (k-1)*kz + i
               nbnode = node + ky
               lk = lctcsr(nbnode, node, ja, ia)
               ld = iau(node)
               val = rhs(node)/a(ld)
!     modify the rhs
               rhs(nbnode) = rhs(nbnode) - a(lk)*val
 31         continue
 30      continue
      endif
!
!     top (north) side
!
      if (alpha(4) .eq. zero) then
         uy = ny - 1
         do 40 k = 1, nz
            do 41 i = lx, ux
               node = (k-1)*kz + i + (ny-1)*ky
               nbnode = node - ky
               lk = lctcsr(nbnode, node, ja, ia)
               ld = iau(node)
               val = rhs(node)/a(ld)
!     modify the rhs
               rhs(nbnode) = rhs(nbnode) - a(lk)*val
 41         continue
 40      continue
      endif
!
!     if only 2-D skip the following section on z
!
      if (nz .le. 1) goto 100
!
!     the front surface
!
      if (alpha(5) .eq. zero) then
         lz = 2
         do 50 j = ly, uy
            do 51 i = lx,  ux
               node = (j-1)*ky + i
               nbnode = node + kz
               lk = lctcsr(nbnode, node, ja, ia)
               ld = iau(node)
               val = rhs(node)/a(ld)
!     modify the rhs
               rhs(nbnode) = rhs(nbnode) - a(lk)*val
 51         continue
 50      continue
      endif
!
!     rear surface
!
      if (alpha(6) .eq. zero) then
         uz = nz - 1
         do 60 j = ly, uy
            do 61 i = lx, ux
               node = (nz-1)*kz + (j-1)*ky + i
               nbnode = node - kz
               lk = lctcsr(nbnode, node, ja, ia)
               ld = iau(node)
               val = rhs(node)/a(ld)
!     modify the rhs
               rhs(nbnode) = rhs(nbnode) - a(lk)*val
 61         continue
 60      continue
      endif
!
!     now the second part ----------------------------------------------
!
!     go through all the actual nodes with unknown values, collect all
!     of them to form a new matrix in compressed sparse row format.
!
 100  kx = 1
      ky = ux - lx + 1
      kz = (uy - ly + 1) * ky
      node = 1
      iedge = 1
      do 80 k = lz, uz
         do 81 j = ly, uy
            do 82 i = lx, ux
!
!     the corresponding old node number
               nbnode = ((k-1)*ny + j-1)*nx + i
!
!     copy the row into local stencil, copy is done is the exact
!     same order as the stencil is written into array a
               lk = ia(nbnode)
               if (i.gt.1) then
                  stencil(2) = a(lk)
                  lk = lk + 1
               end if
               if (j.gt.1) then
                  stencil(4) = a(lk)
                  lk = lk + 1
               end if
               if (k.gt.1) then
                  stencil(6) = a(lk)
                  lk = lk + 1
               end if
               stencil(1) = a(lk)
               lk = lk + 1
               if (i.lt.nx) then
                  stencil(3) = a(lk)
                  lk = lk + 1
               endif
               if (j.lt.ny) then
                  stencil(5) = a(lk)
                  lk = lk + 1
               end if
               if (k.lt.nz) stencil(7) = a(lk)
!
!     first the ia pointer -- points to the beginning of each row
               ia(node) = iedge
!
!     move the values from the local stencil to the new matrix
!
!     the neighbor on the left (west)
               if (i.gt.lx) then
                  ja(iedge)=node-kx
                  a(iedge) =stencil(2)
                  iedge=iedge + 1
               end if
!     the neighbor below (south)
               if (j.gt.ly) then
                  ja(iedge)=node-ky
                  a(iedge)=stencil(4)
                  iedge=iedge + 1
               end if
!     the neighbor in the front
               if (k.gt.lz) then
                  ja(iedge)=node-kz
                  a(iedge)=stencil(6)
                  iedge=iedge + 1
               endif
!     center node (itself)
               ja(iedge) = node
               iau(node) = iedge
               a(iedge) = stencil(1)
               iedge = iedge + 1
!     the neighbor to the right (east)
               if (i.lt.ux) then
                  ja(iedge)=node+kx
                  a(iedge)=stencil(3)
                  iedge=iedge + 1
               end if
!     the neighbor above (north)
               if (j.lt.uy) then
                  ja(iedge)=node+ky
                  a(iedge)=stencil(5)
                  iedge=iedge + 1
               end if
!     the neighbor at the back
               if (k.lt.uz) then
                  ja(iedge)=node+kz
                  a(iedge)=stencil(7)
                  iedge=iedge + 1
               end if
!     the right-hand side
               rhs(node) = rhs(nbnode)
!------next node -------------------------
               node=node+1
!
 82         continue
 81      continue
 80   continue
!
      ia(node) = iedge
!
!     the number of nodes in the final matrix is stored in n
!
      n = node - 1
      return
!-----------------------------------------------------------------------
      end
!-----end of fdreduce-----------------------------------------------------
!-----------------------------------------------------------------------
      subroutine fdaddbc(nx,ny,nz,a,ja,ia,iau,rhs,al,h)
      integer nx, ny, nz, ia(nx*ny*nz), ja(7*nx*ny*nz), iau(nx*ny*nz)
      real*8  h, al(6), a(7*nx*ny*nz), rhs(nx*ny*nz)
!-----------------------------------------------------------------------
! This subroutine will add the boundary condition to the linear system
! consutructed without considering the boundary conditions
!
! The Boundary condition is specified in the following form:
!           du
!     alpha -- + beta u = gamma
!           dn
! Alpha is stored in array AL.  The six side of the boundary appares
! in AL in the following order: left(west), right(east), bottom(south),
! top(north), front, back(rear). (see also the illustration in gen57pt)
! Beta and gamma appears as the functions, betfun and gamfun.
! They have the following prototype
!
! real*8 function xxxfun(x, y, z)
! real*8 x, y, z
!
! where x, y, z are vales in the range of [0, 1][0, (ny-1)*h]
! [0, (nz-1)*h]
!
! At the corners or boundary lines, the boundary conditions are applied
! in the follow order:
! 1) if one side is Dirichlet boundary condition, the Dirichlet boundary
!    condition is used;
! 2) if more than one sides are Dirichlet, the Direichlet condition
!    specified for X direction boundary will overwrite the one specified
!    for Y direction boundary which in turn has priority over Z
!     direction boundaries.
! 3) when all sides are non-Dirichlet, the average values are used.
!-----------------------------------------------------------------------
!     some constants
!
      real*8   half,zero,one,two
      parameter(half=0.5D0,zero=0.0D0,one=1.0D0,two=2.0D0)
!
!     local variables
!
      character*2 side
      integer  i,j,k,kx,ky,kz,node,nbr,ly,uy,lx,ux
      real*8   coeff, ctr, hhalf, x, y, z
      real*8   afun, bfun, cfun, dfun, efun, ffun, gfun, hfun
      external afun, bfun, cfun, dfun, efun, ffun, gfun, hfun
      real*8   betfun, gamfun
      integer  lctcsr
      external lctcsr, betfun, gamfun
!
      hhalf = half * h
      kx = 1
      ky = nx
      kz = nx*ny
!
!     In 3-D case, we need to go through all 6 faces one by one. If
!     the actual dimension is lower, test on ny is performed first.
!     If ny is less or equals to 1, then the value of nz is not
!     checked.
!-----
!     the surface on the left (west) side
!     Concentrate on the contribution from the derivatives related to x,
!     The terms with derivative of x was assumed to be:
!
!     a(3/2,j,k)*[u(2,j,k)-u(1,j,k)] + a(1/2,j,k)*[u(0,j,k)-u(1,j,k)] +
!     h*d(1,j,k)*[u(2,j,k)-u(0,j,k)]/2
!
!     But they actually are:
!
!     2*{a(3/2,j,k)*[u(2,j,k)-u(1,j,k)] -
!     h*a(1,j,k)*[beta*u(1,j,k)-gamma]/alpha]} +
!     h*h*d(1,j,k)*[beta*u(1,j,k)-gamma]/alpha
!
!     Therefore, in terms of local stencil the right neighbor of a node
!     should be changed to 2*a(3/2,j,k),
!     The matrix never contains the left neighbor on this border, nothing
!     needs to be done about it.
!     The following terms should be added to the center stencil:
!     -a(3/2,j,k) + a(1/2,j,k) + [h*d(1,j,k)-2*a(1,j,k)]*h*beta/alpha
!
!     And these terms should be added to the corresponding right-hand side
!     [h*d(1,j,k)-2*a(1,j,k)]*h*gamma/alpha
!
!     Obviously, the formula do not apply for the Dirichlet Boundary
!     Condition, where alpha will be zero. In that case, we simply set
!     all the elements in the corresponding row to zero(0), then let
!     the diagonal element be beta, and the right-hand side be gamma.
!     Thus the value of u at that point will be set. Later on point
!     like this will be removed from the matrix, since they are of
!     know value before solving the system.(not done in this subroutine)
!
      x = zero
      side = 'x1'
      do 20 k = 1, nz
         z = (k-1)*h
         do 21 j = 1, ny
            y = (j-1)*h
            node = 1+(j-1)*ky+(k-1)*kz
!
!     check to see if it's Dirichlet Boundary condition here
!
            if (al(1) .eq. zero) then
               call clrow(node, a, ja, ia)
               a(iau(node)) = betfun(side,x,y,z)
               rhs(node) = gamfun(side,x,y,z)
            else
!
!     compute the terms formulated above to modify the matrix.
!
!     the right neighbor is stroed in nbr'th posiiton in the a
               nbr = lctcsr(node, node+kx, ja, ia)
!
               coeff = two*afun(x,y,z)
               ctr = (h*dfun(x,y,z) - coeff)*h/al(1)
               rhs(node) = rhs(node) + ctr * gamfun(side,x,y,z)
               ctr = afun(x-hhalf,y,z) + ctr * betfun(side,x,y,z)
               coeff = afun(x+hhalf,y,z)
               a(iau(node)) = a(iau(node)) - coeff + ctr
               a(nbr) = two*coeff
            end if
 21      continue
 20   continue
!
!     the right (east) side boudary, similarly, the contirbution from
!     the terms containing the derivatives of x were assumed to be
!
!     a(nx+1/2,j,k)*[u(nx+1,j,k)-u(nx,j,k)] +
!     a(nx-1/2,j,k)*[u(nx-1,j,k)-u(nx,j,k)] +
!     d(nx,j,k)*[u(nx+1,j,k)-u(nx-1,j,k)]*h/2
!
!     Actualy they are:
!
!     2*{h*a(nx,j,k)*[gamma-beta*u(nx,j,k)]/alpha +
!     a(nx-1/2,j,k)*[u(nx-1,j,k)-u(nx,j,k)]} +
!     h*h*d(nx,j,k)*[gamma-beta*u(nx,j,k)]/alpha
!
!     The left stencil has to be set to 2*a(nx-1/2,j,k)
!
!     The following terms have to be added to the center stencil:
!
!     -a(nx-1/2,j,k)+a(nx+1/2,j,k)-[2*a(nx,j,k)+h*d(nx,j,k)]*beta/alpha
!
!     The following terms have to be added to the right-hand side:
!
!     -[2*a(nx,j,k)+h*d(nx,j,k)]*h*gamma/alpha
!
      x = one
      side = 'x2'
      do 22 k = 1, nz
         z = (k-1)*h
         do 23 j = 1, ny
            y = (j-1)*h
            node = (k-1)*kz + j*ky
!
            if (al(2) .eq. zero) then
               call clrow(node, a, ja, ia)
               a(iau(node)) = betfun(side,x,y,z)
               rhs(node) = gamfun(side,x,y,z)
            else
               nbr = lctcsr(node, node-kx, ja, ia)
!
               coeff = two*afun(x,y,z)
               ctr = (coeff + h*dfun(x,y,z))*h/al(2)
               rhs(node) = rhs(node) - ctr * gamfun(side,x,y,z)
               ctr = afun(x+hhalf,y,z) - ctr * betfun(side,x,y,z)
               coeff = afun(x-hhalf,y,z)
               a(iau(node)) = a(iau(node)) - coeff + ctr
               a(nbr) = two*coeff
            end if
 23      continue
 22   continue
!
!     If only one dimension, return now
!
      if (ny .le. 1) return
!
!     the bottom (south) side suface, This similar to the situation
!     with the left side, except all the function and realted variation
!     should be on the y.
!
!     These two block if statment here is to resolve the possible conflict
!     of assign the boundary value differently by different side of the
!     Dirichlet Boundary Conditions. They ensure that the edges that have
!     be assigned a specific value will not be reassigned.
!
      if (al(1) .eq. zero) then
         lx = 2
      else
         lx = 1
      end if
      if (al(2) .eq. zero) then
         ux = nx-1
      else
         ux = nx
      end if
      y = zero
      side = 'y1'
      do 24 k = 1, nz
         z = (k-1)*h
         do 25 i = lx, ux
            x = (i-1)*h
            node = i + (k-1)*kz
!
            if (al(3) .eq. zero) then
               call clrow(node, a, ja, ia)
               a(iau(node)) = betfun(side,x,y,z)
               rhs(node) = gamfun(side,x,y,z)
            else
               nbr = lctcsr(node, node+ky, ja, ia)
!
               coeff = two*bfun(x,y,z)
               ctr = (h*efun(x,y,z) - coeff)*h/al(3)
               rhs(node) = rhs(node) + ctr * gamfun(side,x,y,z)
               ctr = bfun(x,y-hhalf,z) + ctr * betfun(side,x,y,z)
               coeff = bfun(x,y+hhalf,z)
               a(iau(node)) = a(iau(node)) - coeff + ctr
               a(nbr) = two*coeff
            end if
 25      continue
 24   continue
!
!     The top (north) side, similar to the right side
!
      y = (ny-1) * h
      side = 'y2'
      do 26 k = 1, nz
         z = (k-1)*h
         do 27 i = lx, ux
            x = (i-1)*h
            node = (k-1)*kz+(ny-1)*ky + i
!
            if (al(4) .eq. zero) then
               call clrow(node, a, ja, ia)
               a(iau(node)) = betfun(side,x,y,z)
               rhs(node) = gamfun(side,x,y,z)
            else
               nbr = lctcsr(node, node-ky, ja, ia)
!
               coeff = two*bfun(x,y,z)
               ctr = (coeff + h*efun(x,y,z))*h/al(4)
               rhs(node) = rhs(node) - ctr * gamfun(side,x,y,z)
               ctr = bfun(x,y+hhalf,z) - ctr * betfun(side,x,y,z)
               coeff = bfun(x,y-hhalf,z)
               a(iau(node)) = a(iau(node)) - coeff + ctr
               a(nbr) = two*coeff
            end if
 27      continue
 26   continue
!
!     If only has two dimesion to work on, return now
!
      if (nz .le. 1) return
!
!     The front side boundary
!
!     If the edges of the surface has been decided by Dirichlet Boundary
!     Condition, then leave them alone.
!
      if (al(3) .eq. zero) then
         ly = 2
      else
         ly = 1
      end if
      if (al(4) .eq. zero) then
         uy = ny-1
      else
         uy = ny
      end if
!
      z = zero
      side = 'z1'
      do 28 j = ly, uy
         y = (j-1)*h
         do 29 i = lx, ux
            x = (i-1)*h
            node = i + (j-1)*ky
!
            if (al(5) .eq. zero) then
               call clrow(node, a, ja, ia)
               a(iau(node)) = betfun(side,x,y,z)
               rhs(node) = gamfun(side,x,y,z)
            else
               nbr = lctcsr(node, node+kz, ja, ia)
!
               coeff = two*cfun(x,y,z)
               ctr = (h*ffun(x,y,z) - coeff)*h/al(5)
               rhs(node) = rhs(node) + ctr * gamfun(side,x,y,z)
               ctr = cfun(x,y,z-hhalf) + ctr * betfun(side,x,y,z)
               coeff = cfun(x,y,z+hhalf)
               a(iau(node)) = a(iau(node)) - coeff + ctr
               a(nbr) = two*coeff
            end if
 29      continue
 28   continue
!
!     Similiarly for the top side of the boundary suface
!
      z = (nz - 1) * h
      side = 'z2'
      do 30 j = ly, uy
         y = (j-1)*h
         do 31 i = lx, ux
            x = (i-1)*h
            node = (nz-1)*kz + (j-1)*ky + i
!
            if (al(6) .eq. zero) then
               call clrow(node, a, ja, ia)
               a(iau(node)) = betfun(side,x,y,z)
               rhs(node) = gamfun(side,x,y,z)
            else
               nbr = lctcsr(node, node-kz, ja, ia)
!
               coeff = two*cfun(x,y,z)
               ctr = (coeff + h*ffun(x,y,z))*h/al(6)
               rhs(node) = rhs(node) - ctr * gamfun(side,x,y,z)
               ctr = cfun(x,y,z+hhalf) - ctr * betfun(side,x,y,z)
               coeff = cfun(x,y,z-hhalf)
               a(iau(node)) = a(iau(node)) - coeff + ctr
               a(nbr) = two*coeff
            end if
 31      continue
 30   continue
!
!     all set
!
      return
!-----------------------------------------------------------------------
      end
!-----end of fdaddbc----------------------------------------------------
!-----------------------------------------------------------------------
      subroutine clrow(i, a, ja, ia)
      integer i, ja(*), ia(*), k
      real *8 a(*)
!-----------------------------------------------------------------------
!     clear the row i to all zero, but still keep the structure of the
!     CSR matrix
!-----------------------------------------------------------------------
      do 10 k = ia(i), ia(i+1)-1
         a(k) = 0.0D0
 10   continue
!
      return
!-----end of clrow------------------------------------------------------
      end
!-----------------------------------------------------------------------
      function lctcsr(i,j,ja,ia)
      integer lctcsr, i, j, ja(*), ia(*), k
!-----------------------------------------------------------------------
!     locate the position of a matrix element in a CSR format
!     returns -1 if the desired element is zero
!-----------------------------------------------------------------------
      lctcsr = -1
      k = ia(i)
 10   if (k .lt. ia(i+1) .and. (lctcsr .eq. -1)) then
         if (ja(k) .eq. j) lctcsr = k
         k = k + 1
         goto 10
      end if
!
      return
!-----------------------------------------------------------------------
      end
!-----end of lctcsr-----------------------------------------------------

      subroutine csrcsc (n,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n+1),ja(*),jao(*)
      real*8  a(*),ao(*)
!-----------------------------------------------------------------------
! Compressed Sparse Row     to      Compressed Sparse Column
!
! (transposition operation)   Not in place. 
!----------------------------------------------------------------------- 
! -- not in place --
! this subroutine transposes a matrix stored in a, ja, ia format.
! ---------------
! on entry:
!----------
! n     = dimension of A.
! job   = integer to indicate whether to fill the values (job.eq.1) of the
!         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
!
! ipos  = starting position in ao, jao of the transposed matrix.
!         the iao array takes this into account (thus iao(1) is set to ipos.)
!         Note: this may be useful if one needs to append the data structure
!         of the transpose to that of A. In this case use for example
!                call csrcsc (n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
!         for any other normal usage, enter ipos=1.
! a     = real array of length nnz (nnz=number of nonzero elements in input 
!         matrix) containing the nonzero elements.
! ja    = integer array of length nnz containing the column positions
!         of the corresponding elements in a.
! ia    = integer of size n+1. ia(k) contains the position in a, ja of
!         the beginning of the k-th row.
!
! on return:
! ---------- 
! output arguments:
! ao    = real array of size nzz containing the "a" part of the transpose
! jao   = integer array of size nnz containing the column indices.
! iao   = integer array of size n+1 containing the "ia" index array of
!         the transpose. 
!
!----------------------------------------------------------------------- 
      call csrcsc2 (n,n,job,ipos,a,ja,ia,ao,jao,iao)
      end
!-----------------------------------------------------------------------
      subroutine csrcsc2 (n,n2,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n2+1),ja(*),jao(*)
      real*8  a(*),ao(*)
!-----------------------------------------------------------------------
! Compressed Sparse Row     to      Compressed Sparse Column
!
! (transposition operation)   Not in place. 
!----------------------------------------------------------------------- 
! Rectangular version.  n is number of rows of CSR matrix,
!                       n2 (input) is number of columns of CSC matrix.
!----------------------------------------------------------------------- 
! -- not in place --
! this subroutine transposes a matrix stored in a, ja, ia format.
! ---------------
! on entry:
!----------
! n     = number of rows of CSR matrix.
! n2    = number of columns of CSC matrix.
! job   = integer to indicate whether to fill the values (job.eq.1) of the
!         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
!
! ipos  = starting position in ao, jao of the transposed matrix.
!         the iao array takes this into account (thus iao(1) is set to ipos.)
!         Note: this may be useful if one needs to append the data structure
!         of the transpose to that of A. In this case use for example
!                call csrcsc2 (n,n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
!         for any other normal usage, enter ipos=1.
! a     = real array of length nnz (nnz=number of nonzero elements in input 
!         matrix) containing the nonzero elements.
! ja    = integer array of length nnz containing the column positions
!         of the corresponding elements in a.
! ia    = integer of size n+1. ia(k) contains the position in a, ja of
!         the beginning of the k-th row.
!
! on return:
! ---------- 
! output arguments:
! ao    = real array of size nzz containing the "a" part of the transpose
! jao   = integer array of size nnz containing the column indices.
! iao   = integer array of size n+1 containing the "ia" index array of
!         the transpose. 
!
!----------------------------------------------------------------------- 
!----------------- compute lengths of rows of transp(A) ----------------
      do 1 i=1,n2+1
         iao(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i), ia(i+1)-1 
            j = ja(k)+1
            iao(j) = iao(j)+1
 2       continue 
 3    continue
!---------- compute pointers from lengths ------------------------------
      iao(1) = ipos 
      do 4 i=1,n2
         iao(i+1) = iao(i) + iao(i+1)
 4    continue
!--------------- now do the actual copying ----------------------------- 
      do 6 i=1,n
         do 62 k=ia(i),ia(i+1)-1 
            j = ja(k) 
            next = iao(j)
            if (job .eq. 1)  ao(next) = a(k)
            jao(next) = i
            iao(j) = next+1
 62      continue
 6    continue
!-------------------------- reshift iao and leave ---------------------- 
      do 7 i=n2,1,-1
         iao(i+1) = iao(i)
 7    continue
      iao(1) = ipos
!--------------- end of csrcsc2 ---------------------------------------- 
!-----------------------------------------------------------------------
      end

