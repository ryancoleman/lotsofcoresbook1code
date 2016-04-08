        module scalapack_mod

        public :: DLEN_, CTXT_, M_, N_, MB_, NB_, RSRC_, CSRC_, LLD_
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,      &
     &                   LLD_, MB_, M_, NB_, N_, RSRC_                      &
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,      &
     &                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,     &
     &                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )


        interface

        real*8 function dlamch( cmach )
        implicit none
        character cmach
        end function dlamch

        integer function iceil(inum,idenom)
        implicit none
        integer inum,idenom
        end function iceil

        integer function numroc( n, nb, iproc, isrcproc, nprocs )
        implicit none
        integer n,nb,iproc,isrcproc,nprocs
        end function numroc

        real*8 function pzlange(norm,m,n,A,ia,ja,descA,work)
        implicit none
        integer m,n,ia,ja
        character norm
        integer descA(*)
        real*8 work(*)
        complex*16 A(*)
        end function pzlange

        real*4 function pclange(norm,m,n,A,ia,ja,descA,work)
        implicit none
        integer m,n,ia,ja
        character norm
        integer descA(*)
        real*4 work(*)
        complex*8 A(*)
        end function pclange


        subroutine infog2l(grindx,gcindx,desc,nprow,npcol,                 &
     &                     myrow,mycol,lrindx,lcindx,rsrc,csrc)
        implicit none
        integer grindx,gcindx,nprow,npcol,myrow,mycol
        integer lrindx,lcindx, rsrc,csrc
        integer desc(*)
        end subroutine infog2l

        subroutine pcgeadd( trans, m,n,alpha,A,ia,ja,descA,                &
     &                      beta, B,ib,jb,descB )
        implicit none
        character trans
        integer m,n,ia,ja,ib,jb
        integer descA(*), descB(*)
        complex*8 A(*), B(*)
        complex*8 alpha, beta
        end subroutine pcgeadd


        subroutine pzgeadd( trans, m,n,alpha,A,ia,ja,descA,                &
     &                      beta, B,ib,jb,descB )
        implicit none
        character trans
        integer m,n,ia,ja,ib,jb
        integer descA(*), descB(*)
        complex*16 A(*), B(*)
        complex*16 alpha, beta
        end subroutine pzgeadd

        subroutine pzgemv( trans, m,n, alpha, A,ia,ja,descA,              &
     &        X,ix,jx,descX,incx,  beta, Y,iy,jy,descY,incy )
        implicit none
        character trans
        integer m,n,  ia,ja, ix,jx,incx, iy,jy,incy
        integer descA(*), descX(*), descY(*)
        complex*16 A(*), X(*), Y(*)
        complex*16 alpha, beta
        end subroutine pzgemv

        subroutine pzgemm( transA, transB, m,n,k,                         &
     &     alpha, A,ia,ja,descA,  B,ib,jb,descB,                          &
     &     beta, C,ic,jc,descC )
        implicit none
        character transA, transB
        integer m,n,k, ia,ja, ib,jb, ic,jc
        integer descA(*),descB(*),descC(*)
        complex*16 alpha,beta
        complex*16 A(*),B(*),C(*)
        end subroutine pzgemm

        subroutine pcgetrs( trans, n, nrhs, A, ia, ja, descA, ipiv,       &
     &                      B, ib, jb, descB, info )
        implicit none
        character trans
        integer n,nrhs, ia,ja, ib,jb, info
        integer descA(*), descB(*), ipiv(*)
        complex*8 A(*), B(*)
        end subroutine pcgetrs

        subroutine pcgetrf( m, n, A, ia, ja, descA, ipiv, info)
        implicit none
        integer m, n, ia, ja, info
        integer descA(*), ipiv(*)
        complex*8 A(*)
        end subroutine pcgetrf

        subroutine blacs_gridinfo( icontxt, nprow,npcol, myrow,mycol )
        implicit none
        integer icontxt, nprow,npcol, myrow,mycol
        end subroutine blacs_gridinfo

        integer function indxg2p(indxglob,nb,iproc,isrcproc,nprocs)
        implicit none
        integer indxglob,nb,iproc,isrcproc,nprocs
        end function indxg2p


        subroutine psgetrs( trans, n, nrhs, A, ia, ja, descA, ipiv,       &
     &                      B, ib, jb, descB, info )
        implicit none
        character trans
        integer n,nrhs, ia,ja, ib,jb, info
        integer descA(*), descB(*), ipiv(*)
        real A(*), B(*)
        end subroutine psgetrs

        subroutine psgetrf( m, n, A, ia, ja, descA, ipiv, info)
        implicit none
        integer m, n, ia, ja, info
        integer descA(*), ipiv(*)
        real A(*)
        end subroutine psgetrf



        double precision function pdlange(norm,m,n,A,ia,ja,descA,work)
        implicit none
        integer m,n,ia,ja
        character norm
        integer descA(*)
        double precision work(*)
        double precision A(*)
        end function pdlange



        subroutine pdgeadd( trans, m,n,alpha,A,ia,ja,descA,                &
     &                      beta, B,ib,jb,descB )
        implicit none
        character trans
        integer m,n,ia,ja,ib,jb
        integer descA(*), descB(*)
        double precision A(*), B(*)
        double precision alpha, beta
        end subroutine pdgeadd

        subroutine pdgemv( trans, m,n, alpha, A,ia,ja,descA,              &
     &        X,ix,jx,descX,incx,  beta, Y,iy,jy,descY,incy )
        implicit none
        character trans
        integer m,n,  ia,ja, ix,jx,incx, iy,jy,incy
        integer descA(*), descX(*), descY(*)
        double precision A(*), X(*), Y(*)
        double precision alpha, beta
        end subroutine pdgemv

        subroutine pdgemm( transA, transB, m,n,k,                         &
     &     alpha, A,ia,ja,descA,  B,ib,jb,descB,                          &
     &     beta, C,ic,jc,descC )
        implicit none
        character transA, transB
        integer m,n,k, ia,ja, ib,jb, ic,jc
        integer descA(*),descB(*),descC(*)
        double precision alpha,beta
        double precision A(*),B(*),C(*)
        end subroutine pdgemm



        real function pslange(norm,m,n,A,ia,ja,descA,work)
        implicit none
        integer m,n,ia,ja
        character norm
        integer descA(*)
        real work(*)
        real A(*)
        end function pslange



        subroutine psgeadd( trans, m,n,alpha,A,ia,ja,descA,                &
     &                      beta, B,ib,jb,descB )
        implicit none
        character trans
        integer m,n,ia,ja,ib,jb
        integer descA(*), descB(*)
        real A(*), B(*)
        real alpha, beta
        end subroutine psgeadd

        subroutine psgemv( trans, m,n, alpha, A,ia,ja,descA,              &
     &        X,ix,jx,descX,incx,  beta, Y,iy,jy,descY,incy )
        implicit none
        character trans
        integer m,n,  ia,ja, ix,jx,incx, iy,jy,incy
        integer descA(*), descX(*), descY(*)
        real A(*), X(*), Y(*)
        real alpha, beta
        end subroutine psgemv

        subroutine psgemm( transA, transB, m,n,k,                         &
     &     alpha, A,ia,ja,descA,  B,ib,jb,descB,                          &
     &     beta, C,ic,jc,descC )
        implicit none
        character transA, transB
        integer m,n,k, ia,ja, ib,jb, ic,jc
        integer descA(*),descB(*),descC(*)
        real alpha,beta
        real A(*),B(*),C(*)
        end subroutine psgemm

        end interface

        end module scalapack_mod
