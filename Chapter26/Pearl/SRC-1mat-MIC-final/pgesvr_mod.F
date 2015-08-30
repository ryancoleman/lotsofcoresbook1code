        module pgesvr_mod
        use scalapack_mod
        use prof_mod

        implicit none
        private
        public :: pzgesvr, pdgesvr
        public :: pconvertz2c, pconvertc2z
        public :: pconvertd2s, pconverts2d
        public :: convertz2c, convertc2z
        public :: convertd2s, converts2d
*       ----------------------------------------------
*       use reduced precision and iterative refinement
*       ----------------------------------------------

        integer, parameter :: idebug = 1

        real*4, parameter :: zero4 = 0.0d0
        real*8, parameter :: zero8 = 0.0d0
        integer, parameter :: sp = selected_real_kind(6,30)
        integer, parameter :: dp = selected_real_kind(15,300)




        contains

        subroutine pzgecopy( m,n, A,ia,ja,descA, B,ib,jb,descB)
        implicit none
        integer, intent(in) :: m,n,ia,ja,ib,jb
        integer, intent(in) :: descA(*)
        integer, intent(in) :: descB(*)
        complex(kind=dp), intent(in) :: A(*)
        complex(kind=dp), intent(inout) :: B(*)

        complex(kind=dp) :: alpha,beta

        alpha = 1.0d0
        beta = 0.0d0
        call pzgeadd( 'N', m,n, alpha, A,ia,ja,descA, beta,                 &
     &                 B,ib,jb,descB)

        return
        end subroutine pzgecopy



        subroutine pcgecopy( m,n, A,ia,ja,descA, B,ib,jb,descB)
        implicit none
        integer, intent(in) :: m,n,ia,ja,ib,jb
        integer, intent(in) :: descA(*)
        integer, intent(in) :: descB(*)
        complex(kind=sp), intent(in) :: A(*)
        complex(kind=sp), intent(inout) :: B(*)

        complex(kind=sp) :: alpha,beta

        alpha = 1.0d0
        beta = 0.0d0
        call pcgeadd( 'N', m,n, alpha, A,ia,ja,descA, beta,                 &
     &                 B,ib,jb,descB)

        return
        end subroutine pcgecopy


        subroutine pdgecopy( m,n, A,ia,ja,descA, B,ib,jb,descB)
        implicit none
        integer, intent(in) :: m,n,ia,ja,ib,jb
        integer, intent(in) :: descA(*)
        integer, intent(in) :: descB(*)
        real(kind=dp), intent(in) :: A(*)
        real(kind=dp), intent(inout) :: B(*)

        real(kind=dp) :: alpha,beta

        alpha = 1.0d0
        beta = 0.0d0
        call pdgeadd( 'N', m,n, alpha, A,ia,ja,descA, beta,                 &
     &                 B,ib,jb,descB)

        return
        end subroutine pdgecopy



        subroutine psgecopy( m,n, A,ia,ja,descA, B,ib,jb,descB)
        implicit none
        integer, intent(in) :: m,n,ia,ja,ib,jb
        integer, intent(in) :: descA(*)
        integer, intent(in) :: descB(*)
        real(kind=sp), intent(in) :: A(*)
        real(kind=sp), intent(inout) :: B(*)

        real(kind=sp) :: alpha,beta

        alpha = 1.0d0
        beta = 0.0d0
        call psgeadd( 'N', m,n, alpha, A,ia,ja,descA, beta,                 &
     &                 B,ib,jb,descB)

        return
        end subroutine psgecopy

        subroutine  pconvertz2c( m,n,A,ia,ja,descA, B,ib,jb,descB )
        implicit none

        integer, intent(in) :: m,n, ia,ja, ib,jb
        integer, intent(in) :: descA(DLEN_), descB(DLEN_)
        complex(kind=dp) :: A(*)
        complex(kind=sp) :: B(*)
*       ---------------
*       local variables
*       ---------------
        integer :: i, nnb, mmb
        integer :: istart,iend,isize
        integer :: jstart,jend,jsize
        integer :: iia,jja,  iib,jjb
        integer :: lld, Locp, Locq, lrindx, lcindx, rsrc, csrc
        integer :: myrow,mycol, nprow,npcol

        integer :: descAtmp(DLEN_)

        complex(kind=dp), dimension(:), allocatable :: Atmp
        complex(kind=sp), dimension(:), allocatable :: Artmp

        logical :: isok 
        integer :: Asize, info, ierr
        integer, parameter :: rratio = 2
        integer, parameter :: cratio = 2

#ifdef USE_PROF
        call profstart('pconvertz2c')
#endif

        call blacs_gridinfo( descA(CTXT_), nprow,npcol,myrow,mycol)
        nnb = min(n, rratio*descA(NB_)*npcol)
        mmb = min(m, cratio*descA(MB_)*nprow)

*       --------------------------------------
*       approximately align Atmp with A(ia,ja) 
*       to discourage avoid communication
*       --------------------------------------
        call infog2l( ia,ja,descA, nprow,npcol, myrow,mycol,                &
     &         lrindx, lcindx, rsrc, csrc )
        Locp = numroc( mmb, descA(MB_), myrow, rsrc, nprow)
        Locq = numroc( nnb, descA(NB_), mycol, csrc, npcol)
        lld = max(1, Locp)

*       -----------------------------
*       setup descriptor and storage
*       -----------------------------
        Asize = max(1, Locp*Locq )
        allocate( Atmp( Asize ), Artmp(Asize), stat=ierr )
        if (ierr.ne.0) then
           write(*,*)'pconvertz2c:allocate(Atmp) returns ierr=',ierr
        endif

        do i=1,Asize
          Atmp(i) = 0.0d0
        enddo
        do i=1,Asize
          Artmp(i) =  0.0
        enddo

        call descinit( descAtmp, mmb, nnb, descA(MB_), descA(NB_),           &
     &         rsrc, csrc, descA(CTXT_), lld, info )
        isok = (info.eq.0)
        if (.not.isok) then
                write(*,*) 'pzconvertz2c: descinit return info = ', info
                stop 'error in pzconvertz2c'
        endif


        do jstart=1,n,nnb
        do istart=1,m,mmb

           jend = min(n,  jstart+nnb-1)
           jsize = jend - jstart + 1

           iend = min(m,  istart+mmb-1)
           isize = iend - istart + 1

           iia = ia + (istart-1)
           jja = ja + (jstart-1)
           call pzgecopy( isize,jsize,A,iia,jja,                            &
     &            descA, Atmp,1,1,descAtmp )

           call convertz2c( Locp,Locq, Atmp, descAtmp(LLD_),              &
     &              Artmp, descAtmp(LLD_) )

           iib = ib + (istart-1)
           jjb = jb + (jstart-1)
           call pcgecopy( isize,jsize, Artmp, 1,1, descAtmp,               &
     &                 B, iib,jjb, descB )

        enddo
        enddo


        deallocate( Atmp, Artmp, stat=ierr )
        if (ierr.ne.0) then
           write(*,*) 'pconvertz2c:deallocate(Atmp) return ierr=',ierr
        endif

#ifdef USE_PROF
        call profend('pconvertz2c')
#endif

        return

        end subroutine pconvertz2c


        subroutine  pconvertd2s( m,n,A,ia,ja,descA, B,ib,jb,descB )
        implicit none

        integer, intent(in) :: m,n, ia,ja, ib,jb
        integer, intent(in) :: descA(DLEN_), descB(DLEN_)
        real(kind=dp) :: A(*)
        real(kind=sp) :: B(*)
*       ---------------
*       local variables
*       ---------------
        integer :: i, nnb, mmb
        integer :: istart,iend,isize
        integer :: jstart,jend,jsize
        integer :: iia,jja,  iib,jjb
        integer :: lld, Locp, Locq, lrindx, lcindx, rsrc, csrc
        integer :: myrow,mycol, nprow,npcol

        integer :: descAtmp(DLEN_)

        real(kind=dp), dimension(:), allocatable :: Atmp
        real(kind=sp), dimension(:), allocatable :: Artmp

        logical :: isok 
        integer :: Asize, info, ierr
        integer, parameter :: rratio = 2
        integer, parameter :: cratio = 2

#ifdef USE_PROF
        call profstart('pconvertd2s')
#endif

        call blacs_gridinfo( descA(CTXT_), nprow,npcol,myrow,mycol)
        nnb = min(n, rratio*descA(NB_)*npcol)
        mmb = min(m, cratio*descA(MB_)*nprow)

*       --------------------------------------
*       approximately align Atmp with A(ia,ja) 
*       to discourage avoid communication
*       --------------------------------------
        call infog2l( ia,ja,descA, nprow,npcol, myrow,mycol,                &
     &         lrindx, lcindx, rsrc, csrc )
        Locp = numroc( mmb, descA(MB_), myrow, rsrc, nprow)
        Locq = numroc( nnb, descA(NB_), mycol, csrc, npcol)
        lld = max(1, Locp)

*       -----------------------------
*       setup descriptor and storage
*       -----------------------------
        Asize = max(1, Locp*Locq )
        allocate( Atmp( Asize ), Artmp(Asize), stat=ierr )
        if (ierr.ne.0) then
           write(*,*)'pconvertd2s:allocate(Atmp) returns ierr=',ierr
        endif

        do i=1,Asize
          Atmp(i) = 0.0d0
        enddo
        do i=1,Asize
          Artmp(i) =  0.0
        enddo

        call descinit( descAtmp, mmb, nnb, descA(MB_), descA(NB_),           &
     &         rsrc, csrc, descA(CTXT_), lld, info )
        isok = (info.eq.0)
        if (.not.isok) then
                write(*,*) 'pzconvertd2s: descinit return info = ', info
                stop 'error in pdconvertd2s'
        endif


        do jstart=1,n,nnb
        do istart=1,m,mmb

           jend = min(n,  jstart+nnb-1)
           jsize = jend - jstart + 1

           iend = min(m,  istart+mmb-1)
           isize = iend - istart + 1

           iia = ia + (istart-1)
           jja = ja + (jstart-1)
           call pdgecopy( isize,jsize,A,iia,jja,                            &
     &            descA, Atmp,1,1,descAtmp )

           call convertd2s( Locp,Locq, Atmp, descAtmp(LLD_),                &
     &              Artmp, descAtmp(LLD_) )

           iib = ib + (istart-1)
           jjb = jb + (jstart-1)
           call psgecopy( isize,jsize, Artmp, 1,1, descAtmp,               &
     &                 B, iib,jjb, descB )

        enddo
        enddo


        deallocate( Atmp, Artmp, stat=ierr )
        if (ierr.ne.0) then
           write(*,*) 'pconvertd2s:deallocate(Atmp) return ierr=',ierr
        endif

#ifdef USE_PROF
        call profend('pconvertd2s')
#endif

        return

        end subroutine pconvertd2s



        subroutine  pconvertc2z( m,n,A,ia,ja,descA, B,ib,jb,descB )
        implicit none

        integer, intent(in) :: m,n, ia,ja, ib,jb
        integer, intent(in) :: descA(DLEN_), descB(DLEN_)
        complex(kind=sp) :: A(*)
        complex(kind=dp) :: B(*)
*       ---------------
*       local variables
*       ---------------
        integer :: i, nnb, mmb
        integer :: istart,iend,isize
        integer :: jstart,jend,jsize
        integer :: iia,jja,  iib,jjb
        integer :: lld, Locp, Locq, lrindx, lcindx, rsrc, csrc
        integer :: myrow,mycol, nprow,npcol

        integer :: descAtmp(DLEN_)

        complex(kind=dp), dimension(:), allocatable :: Atmp
        complex(kind=sp), dimension(:), allocatable :: Artmp

        logical :: isok 
        integer :: Asize, info, ierr
        integer, parameter :: rratio = 2
        integer, parameter :: cratio = 2

#ifdef USE_PROF
        call profstart('pconvertc2z')
#endif

        call blacs_gridinfo( descA(CTXT_), nprow,npcol,myrow,mycol)
        nnb = min(n, rratio*descA(NB_)*npcol)
        mmb = min(m, cratio*descA(MB_)*nprow)

*       --------------------------------------
*       approximately align Atmp with A(ia,ja) 
*       to discourage avoid communication
*       --------------------------------------
        call infog2l( ia,ja,descA, nprow,npcol, myrow,mycol,                &
     &         lrindx, lcindx, rsrc, csrc )
        Locp = numroc( mmb, descA(MB_), myrow, rsrc, nprow)
        Locq = numroc( nnb, descA(NB_), mycol, csrc, npcol)
        lld = max(1, Locp)

*       -----------------------------
*       setup descriptor and storage
*       -----------------------------
        Asize = max(1, Locp*Locq )
        allocate( Atmp( Asize ), Artmp(Asize), stat=ierr )
        if (ierr.ne.0) then
          write(*,*) 'pconvertc2z: allocate(Atmp) return ierr=',ierr
        endif

        do i=1,Asize
          Atmp(i) = 0.0d0
        enddo
        do i=1,Asize
          Artmp(i) = 0.0d0
        enddo

        call descinit( descAtmp, mmb, nnb, descA(MB_), descA(NB_),           &
     &         rsrc, csrc, descA(CTXT_), lld, info )
        isok = (info.eq.0)
        if (.not.isok) then
                write(*,*) 'pzc2z: descinit return info = ', info
                stop 'error in pzc2z'
        endif


        do jstart=1,n,nnb
        do istart=1,m,mmb

           jend = min(n,  jstart+nnb-1)
           jsize = jend - jstart + 1

           iend = min(m,  istart+mmb-1)
           isize = iend - istart + 1

           iia = ia + (istart-1)
           jja = ja + (jstart-1)
           call pcgecopy( isize,jsize,A,iia,jja,                            &
     &            descA, Artmp,1,1,descAtmp )

           call convertc2z( Locp,Locq, Artmp, descAtmp(LLD_),              &
     &              Atmp, descAtmp(LLD_) )

           iib = ib + (istart-1)
           jjb = jb + (jstart-1)
           call pzgecopy( isize,jsize, Atmp, 1,1, descAtmp,               &
     &                 B, iib,jjb, descB )

        enddo
        enddo


        deallocate( Atmp, Artmp, stat=ierr )
        if (ierr.ne.0) then
          write(*,*) 'pconvertc2z:deallocate(Atmp) return ierr=',ierr
        endif

#ifdef USE_PROF
        call profend('pconvertc2z')
#endif

        return

        end subroutine pconvertc2z



        subroutine  pconverts2d( m,n,A,ia,ja,descA, B,ib,jb,descB )
        implicit none

        integer, intent(in) :: m,n, ia,ja, ib,jb
        integer, intent(in) :: descA(DLEN_), descB(DLEN_)
        real(kind=sp) :: A(*)
        real(kind=dp) :: B(*)
*       ---------------
*       local variables
*       ---------------
        integer :: i, nnb, mmb
        integer :: istart,iend,isize
        integer :: jstart,jend,jsize
        integer :: iia,jja,  iib,jjb
        integer :: lld, Locp, Locq, lrindx, lcindx, rsrc, csrc
        integer :: myrow,mycol, nprow,npcol

        integer :: descAtmp(DLEN_)

        real(kind=dp), dimension(:), allocatable :: Atmp
        real(kind=sp), dimension(:), allocatable :: Artmp

        logical :: isok 
        integer :: Asize, info, ierr
        integer, parameter :: rratio = 2
        integer, parameter :: cratio = 2

#ifdef USE_PROF
        call profstart('pconverts2d')
#endif

        call blacs_gridinfo( descA(CTXT_), nprow,npcol,myrow,mycol)
        nnb = min(n, rratio*descA(NB_)*npcol)
        mmb = min(m, cratio*descA(MB_)*nprow)

*       --------------------------------------
*       approximately align Atmp with A(ia,ja) 
*       to discourage avoid communication
*       --------------------------------------
        call infog2l( ia,ja,descA, nprow,npcol, myrow,mycol,                &
     &         lrindx, lcindx, rsrc, csrc )
        Locp = numroc( mmb, descA(MB_), myrow, rsrc, nprow)
        Locq = numroc( nnb, descA(NB_), mycol, csrc, npcol)
        lld = max(1, Locp)

*       -----------------------------
*       setup descriptor and storage
*       -----------------------------
        Asize = max(1, Locp*Locq )
        allocate( Atmp( Asize ), Artmp(Asize), stat=ierr )
        if (ierr.ne.0) then
          write(*,*) 'pconverts2d: allocate(Atmp) return ierr=',ierr
        endif

        do i=1,Asize
          Atmp(i) = 0.0d0
        enddo
        do i=1,Asize
          Artmp(i) = 0.0d0
        enddo

        call descinit( descAtmp, mmb, nnb, descA(MB_), descA(NB_),           &
     &         rsrc, csrc, descA(CTXT_), lld, info )
        isok = (info.eq.0)
        if (.not.isok) then
                write(*,*) 'pzc2z: descinit return info = ', info
                stop 'error in pzc2z'
        endif


        do jstart=1,n,nnb
        do istart=1,m,mmb

           jend = min(n,  jstart+nnb-1)
           jsize = jend - jstart + 1

           iend = min(m,  istart+mmb-1)
           isize = iend - istart + 1

           iia = ia + (istart-1)
           jja = ja + (jstart-1)
           call psgecopy( isize,jsize,A,iia,jja,                            &
     &            descA, Artmp,1,1,descAtmp )

           call converts2d( Locp,Locq, Artmp, descAtmp(LLD_),              &
     &              Atmp, descAtmp(LLD_) )

           iib = ib + (istart-1)
           jjb = jb + (jstart-1)
           call pdgecopy( isize,jsize, Atmp, 1,1, descAtmp,               &
     &                 B, iib,jjb, descB )

        enddo
        enddo


        deallocate( Atmp, Artmp, stat=ierr )
        if (ierr.ne.0) then
          write(*,*) 'pconverts2d:deallocate(Atmp) return ierr=',ierr
        endif

#ifdef USE_PROF
        call profend('pconverts2d')
#endif

        return

        end subroutine pconverts2d


        subroutine convertz2c( m,n, A, ldA, B, ldB )
        integer, intent(in) :: m,n,ldA,ldB
        complex(kind=dp),intent(in) :: A(ldA,*)
        complex(kind=sp),intent(inout) :: B(ldB,*)

        integer :: i,j

        do j=1,n
        do i=1,m
          B(i,j) = A(i,j)
        enddo
        enddo

        return
        end subroutine convertz2c




        subroutine convertc2z( m,n, A, ldA, B, ldB )
        integer, intent(in) :: m,n,ldA,ldB
        complex(kind=sp),intent(in) :: A(ldA,*)
        complex(kind=dp),intent(inout) :: B(ldB,*)

        integer :: i,j

        do j=1,n
        do i=1,m
          B(i,j) = A(i,j)
        enddo
        enddo

        return
        end subroutine convertc2z



        subroutine convertd2s( m,n, A, ldA, B, ldB )
        integer, intent(in) :: m,n,ldA,ldB
        real(kind=dp),intent(in) :: A(ldA,*)
        real(kind=sp),intent(inout) :: B(ldB,*)

        integer :: i,j

        do j=1,n
        do i=1,m
          B(i,j) = A(i,j)
        enddo
        enddo

        return
        end subroutine convertd2s




        subroutine converts2d( m,n, A, ldA, B, ldB )
        integer, intent(in) :: m,n,ldA,ldB
        real(kind=sp),intent(in) :: A(ldA,*)
        real(kind=dp),intent(inout) :: B(ldB,*)

        integer :: i,j

        do j=1,n
        do i=1,m
          B(i,j) = A(i,j)
        enddo
        enddo

        return
        end subroutine converts2d


        subroutine pzgesvr( n, nrhs, A, ia, ja, descA,                     &
     &          ipiv, B, ib, jb, descB, info )
        integer, intent(in) :: n,nrhs, ia,ja, ib,jb
        integer, intent(in) :: descA(DLEN_), descB(DLEN_)
        integer, intent(inout) :: info

        integer, intent(inout) :: ipiv(*)
        complex(kind=dp), intent(inout) :: A(*)
        complex(kind=dp), intent(inout) :: B(*)
*       ---------------
*       local variables
*       ---------------
        logical :: isok, isdone, isroot

        integer :: descAr(DLEN_), descX(DLEN_)
        integer :: ictxt, myrow,mycol, nprow,npcol
        integer :: Locp, Locq, lrindx, lcindx, lld, rsrc,csrc
        integer :: Locp_B, Locq_B


        integer :: iter
        integer, parameter :: maxiter = 30

        integer :: Xsize, Asize, ierr
        complex(kind=sp), dimension(:), allocatable :: Ar, Xr
        complex(kind=dp), dimension(:), allocatable :: X, R
        real(kind=dp)  :: anorm, eps, bwdcte, xnorm
        real(kind=dp)  :: rnorm,rnorm_old,bwdmax

        real(kind=sp), dimension(:), allocatable :: rwork
        real(kind=dp), dimension(:), allocatable :: dwork
        integer :: ldwork 
        integer :: m,Mp0,Np0

        complex(kind=dp), parameter :: one = 1.0d0
        complex(kind=dp) :: alpha,beta

        integer :: i
        complex(kind=sp) :: work(1024*1024)
#ifdef USE_HPL
        integer :: hpl_lld, hpl_ineed
#endif

#ifdef USE_GPU
        logical, parameter :: use_scalapack = .false.
#else
        logical, parameter :: use_scalapack = .true.
#endif

        character(len=255) :: memsize_str
        integer  istatus, ivalue
        integer idebug, memsize, memsize_default
        parameter(idebug = 1,memsize_default=128*1024*1024)

*       ----------------------------------------------
*       norm = 'F','f','E', or 'e' for Frobenius norm
*       norm = 'M','m' for max norm
*       norm = '1','O','o' for 1-norm
*       norm = 'I','i', for infinity norm
*       ----------------------------------------------
        character :: normtype = 'I'   

        ictxt = descA(CTXT_)
        call blacs_gridinfo( ictxt, nprow,npcol,myrow,mycol )
        isroot = (myrow.eq.0).and.(mycol.eq.0)

        m = n
        Mp0 = iceil( m, descA(MB_) )*descA(MB_)
        Np0 = iceil( n, descA(NB_) )*descA(NB_)
        ldwork = max(Mp0,Np0) + 1
        allocate( dwork(ldwork),rwork(ldwork),stat=ierr)
        if (ierr.ne.0) then
           write(*,*) 'allocate(dwork) return ierr=',ierr
        endif

         anorm = pzlange(normtype, n,n, A,ia,ja,descA, dwork )
         eps = dlamch('P')
         bwdmax = 1.0d0
*        bwdcte = anorm * eps * min(4.0d0,sqrt(dble(n))/6.0d0)
*        bwdcte = anorm * eps * max(4.0d0,sqrt(dble(n))/6.0d0)
         bwdcte = anorm * eps * sqrt( dble(n) ) * bwdmax

      

*       ----------------------------------
*       make a copy with reduced precision
*       ----------------------------------


*       ------------------------------- 
*       slightly over-estimate storage
*       ------------------------------- 
        rsrc = myrow
        csrc = mycol
        Locp = numroc( n, descA(MB_), myrow, rsrc, nprow)
        Locq = numroc( n+1, descA(NB_), mycol, csrc, npcol)
        lld = max(1,Locp)
        Asize = max(1,Locp*Locq)

*       ----------------------------------------------------
*       set rsrc,csrc so data is exactly aligned
*       on the same processor, no communication is required
*       ----------------------------------------------------
        call infog2l( ia,ja, descA, nprow,npcol, myrow,mycol,               &
     &         lrindx, lcindx, rsrc, csrc )

*       ---------------------------------------
*       for simplicity, use same (rsrc,csrc) as 
*       processor owning (ia,ja)
*       ---------------------------------------


#ifdef USE_HPL
        call hpl_cmatinit( n, descA(MB_), hpl_lld, hpl_ineed )

        if (idebug .ge. 2) then
           write(*,*) 'myrow,mycol,hpl_lld,hpl_ineed ',                     &
     &              myrow,mycol,hpl_lld,hpl_ineed
           write(*,*) 'n,mb,nb ', n, descA(MB_),descA(NB_)
        endif

        call assert(hpl_lld.ge.lld, 'hpl_lld < lld ', hpl_lld )
        call assert(hpl_ineed.ge.Asize,'hpl_ineed < Asize ',hpl_ineed)

        lld = hpl_lld
        size = hpl_ineed
#endif

        allocate( Ar( Asize ), stat=ierr)
        if (ierr.ne.0) then
          write(*,*) 'allocate(Ar) return ierr=',ierr
        endif
        Ar(:) = 0
     

        call descinit( descAr, n,n+1, descA(MB_), descA(NB_),                 &
     &       rsrc, csrc, descA(CTXT_), lld, info )
        isok = (info.eq.0)
        if (.not.isok) then
                write(*,*) 'decinit for Ar return info=',info
                stop '** error in pzgesvr '
        endif

#ifdef USE_PROF
        call blacs_barrier( descAr(CTXT_),'A')
        call profstart('pzgesvr: pconvertz2c A')
#endif

        
        call pconvertz2c( n,n, A,ia,ja,descA, Ar,1,1,descAr )

#ifdef USE_PROF
        call blacs_barrier( descAr(CTXT_),'A')
        call profend('pzgesvr: pconvertz2c A')
#endif


*       ------------------
*       make a copy of rhs
*       ------------------
        rsrc = myrow
        csrc = mycol
        Locp_B = numroc( n, descB(MB_), myrow, rsrc, nprow)
        Locq_B = numroc( nrhs, descB(NB_), mycol, csrc, npcol)
        Xsize = max(1, Locp_B*Locq_B)
        allocate( X( Xsize ), R(Xsize), Xr(Xsize), stat=ierr )
        if (ierr.ne.0) then
          write(*,*) 'allocate(X) return ierr=',ierr
        endif

*       ----------------
*       setup descriptor
*       ----------------
        call infog2l( ib,jb, descB, nprow,npcol, myrow,mycol,              &
     &          lrindx, lcindx, rsrc, csrc )

        lld = max(1, Locp_B )
        call descinit( descX, n, nrhs,  descB(MB_), descB(NB_),            &
     &            rsrc, csrc, descB(CTXT_), lld, info )


        isok = (info.eq.0)
        if (.not.isok) then
                write(*,*) 'descinit for X return info=',info
                stop '** error in pzgesvr '
        endif

        call pzgecopy( n,nrhs, B,ib,jb,descB, X,1,1,descX )

        call pconvertz2c( n,nrhs, X,1,1,descX, Xr,1,1,descX )

        if ((idebug.ge.2) .and. isroot) then
*       ------------
*       double check
*       ------------
        rnorm = pclange(normtype,n,nrhs,Xr,1,1,descX,rwork)
        if ((myrow.eq.0).and.(mycol.eq.0)) then
          write(*,*) 'norm(Xr) after c2z ', rnorm
        endif


        rnorm = pzlange(normtype, n,nrhs, B,ib,jb,descB, dwork )
        if ((myrow.eq.0).and.(mycol.eq.0)) then
          write(*,*) 'norm(B) ', rnorm
        endif
        call pconvertc2z( n,nrhs, Xr,1,1,descX, R,1,1,descX )


        rnorm = pzlange(normtype,n,nrhs,R,1,1,descX,dwork)
        if ((myrow.eq.0).and.(mycol.eq.0)) then
          write(*,*) 'norm(R) after c2z ', rnorm
        endif

        alpha = -1.0d0
        beta = 1.0d0
        call pzgeadd( 'N', n,nrhs, alpha, X,1,1,descX,                     &
     &                 beta,R,1,1,descX)
        rnorm = pzlange( normtype, n,nrhs, R,1,1,descX,dwork)
        if ((myrow.eq.0).and.(mycol.eq.0)) then
          write(*,*) 'rnorm from conversion ',rnorm
        endif

        endif
*       ---------------------------------------------
*       perform LU factorization in reduced precision
*       ---------------------------------------------

#ifdef USE_HPL

#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profstart('hpl_pcgesv')
#endif

        call hpl_pcgesv(n, Ar, descAr,ipiv,info)
        call assert(info.eq.0,'hpl_pcgesv return info',info)
        
#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profend('hpl_pcgesv')
#endif

#else




#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profstart('pcgetrf')
#endif
        if (use_scalapack) then
          call pcgetrf( n,n, Ar, 1,1, descAr, ipiv, info )
        else
          memsize = memsize_default
          call getenv("MEMSIZE",memsize_str)
          if (len(trim(memsize_str)).ge.1) then
            ivalue = 128
            read(memsize_str,*,iostat=istatus) ivalue
            if ((istatus.eq.0) .and. (ivalue .ge.1)) then
              memsize = ivalue*1024*1024
            endif
          endif

          if ((myrow.eq.0).and.(mycol.eq.0)) then
                    write(*,*) 'memsize = ', memsize
          endif
!debug           call pclaprnt(n,n,Ar,1,1,descAr,0,0,'Ar',17,work)
          if ((myrow.eq.0).and.(mycol.eq.0)) then
             write(*,*) 'before calling pcgetrf_ooc2'
          endif
          call pcgetrf_ooc2(n,n,Ar,1,1,descAr,ipiv, memsize,info)
          !call pcgetrf(n,n,Ar,1,1,descAr,ipiv, info)

!debug
!debug           if ((myrow.eq.0).and.(mycol.eq.0)) then
!debug             write(*,*) 'ipiv'
!debug             do i=1,n
!debug               write(*,*) 'ipiv(',i,') = ',ipiv(i)
!debug             enddo
!debug           endif
!debug           call pclaprnt(n,n,Ar,1,1,descAr,0,0,'LU',16,work)



        endif
#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profend('pcgetrf')
#endif


#endif

        isok = (info.eq.0)
        if (.not.isok) then
*          ---------------------------------
*          system may be too ill conditioned
*          ---------------------------------
           if ((myrow.eq.0).and.(mycol.eq.0)) then
             write(*,*) 'pcgetrf return info=',info
           endif
           info  = -(n+1)
           goto 999
        endif
*       -----------------------
*       generate  initial guess
*       -----------------------
        call pconvertz2c( n,nrhs, B,ib,jb,descB, Xr,1,1,descX )

#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profstart('pcgetrs')
#endif
        call pcgetrs( 'N', n, nrhs, Ar, 1,1, descAr, ipiv,                 &
     &        Xr, 1,1, descX, info )
#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profend('pcgetrs')
#endif


        call pconvertc2z( n,nrhs, Xr,1,1,descX,X,1,1,descX)

        if ((idebug.ge.2).and.isroot) then
                write(*,*) 'initial B(1:5) ', B(1:5)
                write(*,*) 'initial Xr(1:5) ', Xr(1:5)
                write(*,*) 'initial X(1:5) ', X(1:5)
        endif

*       ----------------        
*       compute residual
*       ----------------        
        call pzgecopy( n,nrhs, B, ib,jb, descB, R, 1,1,descX )

        alpha = -one
        beta = one

#ifdef USE_PROF
        call profstart('pzgesvr:pzgemm')
#endif
        call pzgemm('N','N',n,nrhs,n, alpha, A,ia,ja,descA,                 &
     &        X,1,1,descX,  beta, R, 1,1, descX )
#ifdef USE_PROF
        call profend('pzgesvr:pzgemm')
#endif

        xnorm = pzlange(normtype,n,nrhs,B,ib,jb,descB,dwork)
        rnorm_old = xnorm
        do iter=1,maxiter
*          ---------------------
*          check for convergence
*          ---------------------
           rnorm =  pzlange( normtype, n,nrhs, R,1,1,descX, dwork )
           isdone = (rnorm.le.xnorm*bwdcte) .and. (iter.ge.6)

           if ((idebug.ge.1).and.isroot) then
             write(*,*) 'iter, rnorm ', iter,rnorm
           endif

           if (isdone) then
                   exit
           endif

           if ((iter.ge.7).and.(rnorm.gt.rnorm_old)) then
              exit
           endif
*          ------------------
*          compute correction
*          ------------------
           call pconvertz2c( n,nrhs, R,1,1, descX, Xr,1,1,descX )
#ifdef USE_PROF
           call profstart('pcgetrs')
#endif
           call pcgetrs( 'N', n,nrhs, Ar, 1,1, descAr, ipiv,               &
     &              Xr,1,1,descX,info )
#ifdef USE_PROF
           call profend('pcgetrs')
#endif
*          ----------------
*          update solution
*          ----------------
           call pconvertc2z( n,nrhs, Xr, 1,1,descX, R,1,1,descX )
           if ((idebug.ge.2) .and.(isroot)) then
             rnorm = pclange(normtype,n,nrhs, Xr,1,1,descX,rwork)
             write(*,*) 'norm(Xr correction) ', rnorm

             rnorm = pzlange(normtype,n,nrhs, R,1,1,descX,dwork)
             write(*,*) 'norm( correction ) ',  rnorm
             write(*,*) 'R(1:5) ', R(1:5)
             write(*,*) 'Xr(1:5) ', Xr(1:5)
             write(*,*) 'X(1:5) ', X(1:5)
           endif

           alpha =   one
           beta = one
           call pzgeadd( 'N', n,nrhs, alpha, R,1,1,descX,                  &
     &                    beta, X,1,1,descX )

*          ----------------
*          compute residual
*          ----------------
           call pzgecopy( n,nrhs, B,ib,jb,descB, R, 1,1, descX )
           alpha = -one
           beta = one
#ifdef USE_PROF
           call profstart('pzgesvr:pzgemm')
#endif
           call pzgemm('N','N',n,nrhs,n,                                    &
     &          alpha, A,ia,ja,descA, X,1,1,descX,                          &
     &          beta, R,1,1,descX )
#ifdef USE_PROF
           call profend('pzgesvr:pzgemm')
#endif
         rnorm_old = rnorm
        enddo

*       ---------------------
*       finally copy solution
*       ---------------------
        call pzgecopy( n,nrhs, X,1,1,descX, B,ib,jb,descB )
        info = 0
        if (.not.isdone) then
                info = -(n+1)
        endif

 999    continue

        deallocate( Ar, X, R, Xr, dwork,rwork,stat=ierr)
        if (ierr.ne.0) then
         write(*,*) 'deallocate(Ar) return ierr=',ierr
        endif

        return

        end subroutine pzgesvr



        subroutine pdgesvr( n, nrhs, A, ia, ja, descA,                     &
     &          ipiv, B, ib, jb, descB, info )
        integer, intent(in) :: n,nrhs, ia,ja, ib,jb
        integer, intent(in) :: descA(DLEN_), descB(DLEN_)
        integer, intent(inout) :: info

        integer, intent(inout) :: ipiv(*)
        real(kind=dp), intent(inout) :: A(*)
        real(kind=dp), intent(inout) :: B(*)
*       ---------------
*       local variables
*       ---------------
        logical :: isok, isdone, isroot

        integer :: descAr(DLEN_), descX(DLEN_)
        integer :: ictxt, myrow,mycol, nprow,npcol
        integer :: Locp, Locq, lrindx, lcindx, lld, rsrc,csrc
        integer :: Locp_B, Locq_B


        integer :: iter
        integer, parameter :: maxiter = 30

        integer :: Xsize, Asize, ierr
        real(kind=sp), dimension(:), allocatable :: Ar, Xr
        real(kind=dp), dimension(:), allocatable :: X, R
        real(kind=dp)  :: anorm, eps, bwdcte, xnorm
        real(kind=dp)  :: rnorm,rnorm_old,bwdmax

        real(kind=sp), dimension(:), allocatable :: rwork
        real(kind=dp), dimension(:), allocatable :: dwork
        integer :: ldwork 
        integer :: m,Mp0,Np0

        real(kind=dp), parameter :: one = 1.0d0
        real(kind=dp) :: alpha,beta

#ifdef USE_HPL
        integer :: hpl_lld, hpl_ineed
#endif

*       ----------------------------------------------
*       norm = 'F','f','E', or 'e' for Frobenius norm
*       norm = 'M','m' for max norm
*       norm = '1','O','o' for 1-norm
*       norm = 'I','i', for infinity norm
*       ----------------------------------------------
        character :: normtype = 'I'   

        ictxt = descA(CTXT_)
        call blacs_gridinfo( ictxt, nprow,npcol,myrow,mycol )
        isroot = (myrow.eq.0).and.(mycol.eq.0)

        m = n
        Mp0 = iceil( m, descA(MB_) )*descA(MB_)
        Np0 = iceil( n, descA(NB_) )*descA(NB_)
        ldwork = max(Mp0,Np0) + 1
        allocate( dwork(ldwork),rwork(ldwork),stat=ierr)
        if (ierr.ne.0) then
           write(*,*) 'allocate(dwork) return ierr=',ierr
        endif

         anorm = pdlange(normtype, n,n, A,ia,ja,descA, dwork )
         eps = dlamch('P')
         bwdmax = 1.0d0
*        bwdcte = anorm * eps * min(4.0d0,sqrt(dble(n))/6.0d0)
*        bwdcte = anorm * eps * max(4.0d0,sqrt(dble(n))/6.0d0)
         bwdcte = anorm * eps * sqrt( dble(n) ) * bwdmax

      

*       ----------------------------------
*       make a copy with reduced precision
*       ----------------------------------

*       ----------------------------------------------------
*       set rsrc,csrc so data is exactly aligned
*       on the same processor, no communication is required
*       ----------------------------------------------------
        call infog2l( ia,ja, descA, nprow,npcol, myrow,mycol,               &
     &         lrindx, lcindx, rsrc, csrc )

*       ---------------------------------------
*       for simplicity, use (rsrc,csrc) = (0,0)
*       ---------------------------------------
        rsrc = 0
        csrc = 0

        Locp = numroc( n, descA(MB_), myrow, rsrc, nprow)
        Locq = numroc( n+1, descA(NB_), mycol, csrc, npcol)
        lld = max(1,Locp)
        Asize = max(1,Locp*Locq)


#ifdef USE_HPL
        call hpl_smatinit( n, descA(MB_), hpl_lld, hpl_ineed )

        if (idebug .ge. 2) then
           write(*,*) 'myrow,mycol,hpl_lld,hpl_ineed ',                     &
     &              myrow,mycol,hpl_lld,hpl_ineed
           write(*,*) 'n,mb,nb ', n, descA(MB_),descA(NB_)
        endif

        call assert(hpl_lld.ge.lld, 'hpl_lld < lld ', hpl_lld )
        call assert(hpl_ineed.ge.Asize,'hpl_ineed < Asize ',hpl_ineed)

        lld = hpl_lld
        Asize = hpl_ineed
#endif

        allocate( Ar( Asize ), stat=ierr)
        if (ierr.ne.0) then
          write(*,*) 'allocate(Ar) return ierr=',ierr
        endif
        Ar(:) = 0
     

        call descinit( descAr, n,n+1, descA(MB_), descA(NB_),                 &
     &       rsrc, csrc, descA(CTXT_), lld, info )
        isok = (info.eq.0)
        if (.not.isok) then
                write(*,*) 'decinit for Ar return info=',info
                stop '** error in pdgesvr '
        endif

#ifdef USE_PROF
        call blacs_barrier( descAr(CTXT_),'A')
        call profstart('pzgesvr: pconvertd2s A')
#endif

        
        call pconvertd2s( n,n, A,ia,ja,descA, Ar,1,1,descAr )

#ifdef USE_PROF
        call blacs_barrier( descAr(CTXT_),'A')
        call profend('pzgesvr: pconvertd2s A')
#endif


*       ------------------
*       make a copy of rhs
*       ------------------
        call infog2l( ib,jb, descB, nprow,npcol, myrow,mycol,              &
     &          lrindx, lcindx, rsrc, csrc )
        Locp_B = numroc( n, descB(MB_), myrow, rsrc, nprow)
        Locq_B = numroc( nrhs, descB(NB_), mycol, csrc, npcol)
        Xsize = max(1, Locp_B*Locq_B)
        allocate( X( Xsize ), R(Xsize), Xr(Xsize), stat=ierr )
        if (ierr.ne.0) then
          write(*,*) 'allocate(X) return ierr=',ierr
        endif

*       ----------------
*       setup descriptor
*       ----------------
        lld = max(1, Locp_B )
        call descinit( descX, n, nrhs,  descB(MB_), descB(NB_),            &
     &            rsrc, csrc, descB(CTXT_), lld, info )


        isok = (info.eq.0)
        if (.not.isok) then
                write(*,*) 'descinit for X return info=',info
                stop '** error in pdgesvr '
        endif

        call pdgecopy( n,nrhs, B,ib,jb,descB, X,1,1,descX )

        call pconvertd2s( n,nrhs, X,1,1,descX, Xr,1,1,descX )

        if ((idebug.ge.2) .and. isroot) then
*       ------------
*       double check
*       ------------
        rnorm = pslange(normtype,n,nrhs,Xr,1,1,descX,rwork)
        write(*,*) 'norm(Xr) after c2z ', rnorm


        rnorm = pdlange(normtype, n,nrhs, B,ib,jb,descB, dwork )
        write(*,*) 'norm(B) ', rnorm
        call pconverts2d( n,nrhs, Xr,1,1,descX, R,1,1,descX )


        rnorm = pdlange(normtype,n,nrhs,R,1,1,descX,dwork)
        write(*,*) 'norm(R) after c2z ', rnorm

        alpha = -1.0d0
        beta = 1.0d0
        call pdgeadd( 'N', n,nrhs, alpha, X,1,1,descX,                     &
     &                 beta,R,1,1,descX)
        rnorm = pdlange( normtype, n,nrhs, R,1,1,descX,dwork)
        write(*,*) 'rnorm from conversion ',rnorm

        endif
*       ---------------------------------------------
*       perform LU factorization in reduced precision
*       ---------------------------------------------

#ifdef USE_HPL

#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profstart('hpl_psgesv')
#endif

        call hpl_psgesv(n, Ar, descAr,ipiv,info)
        call assert(info.eq.0,'hpl_psgesv return info',info)
        
#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profend('hpl_psgesv')
#endif

#else




#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profstart('psgetrf')
#endif
        call psgetrf( n,n, Ar, 1,1, descAr, ipiv, info )
#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profend('psgetrf')
#endif


#endif

        isok = (info.eq.0)
        if (.not.isok) then
*          ---------------------------------
*          system may be too ill conditioned
*          ---------------------------------
           info  = -(n+1)
           goto 999
        endif
*       -----------------------
*       generate  initial guess
*       -----------------------
        call pconvertd2s( n,nrhs, B,ib,jb,descB, Xr,1,1,descX )

#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profstart('psgetrs')
#endif
        call psgetrs( 'N', n, nrhs, Ar, 1,1, descAr, ipiv,                 &
     &        Xr, 1,1, descX, info )
#ifdef USE_PROF
        call blacs_barrier(ictxt,'A')
        call profend('psgetrs')
#endif


        call pconverts2d( n,nrhs, Xr,1,1,descX,X,1,1,descX)

        if ((idebug.ge.2).and.isroot) then
                write(*,*) 'initial B(1:5) ', B(1:5)
                write(*,*) 'initial Xr(1:5) ', Xr(1:5)
                write(*,*) 'initial X(1:5) ', X(1:5)
        endif

*       ----------------        
*       compute residual
*       ----------------        
        call pdgecopy( n,nrhs, B, ib,jb, descB, R, 1,1,descX )

        alpha = -one
        beta = one

#ifdef USE_PROF
        call profstart('pzgesvr:pzgemm')
#endif
        call pdgemm('N','N',n,nrhs,n, alpha, A,ia,ja,descA,                 &
     &        X,1,1,descX,  beta, R, 1,1, descX )
#ifdef USE_PROF
        call profend('pzgesvr:pzgemm')
#endif

        xnorm = pdlange(normtype,n,nrhs,B,ib,jb,descB,dwork)
        rnorm_old = xnorm
        do iter=1,maxiter
*          ---------------------
*          check for convergence
*          ---------------------
           rnorm =  pdlange( normtype, n,nrhs, R,1,1,descX, dwork )
           isdone = (rnorm.le.xnorm*bwdcte) .and. (iter.ge.6)

           if ((idebug.ge.1).and.isroot) then
             write(*,*) 'iter,rnorm,rnorm_old ',iter,rnorm,rnorm_old
           endif

           if (isdone) then
                   exit
           endif

           if ((iter.ge.7).and.(rnorm.gt.rnorm_old)) then
              exit
           endif
*          ------------------
*          compute correction
*          ------------------
           call pconvertd2s( n,nrhs, R,1,1, descX, Xr,1,1,descX )
#ifdef USE_PROF
           call profstart('psgetrs')
#endif
           call psgetrs( 'N', n,nrhs, Ar, 1,1, descAr, ipiv,               &
     &              Xr,1,1,descX,info )
#ifdef USE_PROF
           call profend('psgetrs')
#endif
*          ----------------
*          update solution
*          ----------------
           call pconverts2d( n,nrhs, Xr, 1,1,descX, R,1,1,descX )
           if ((idebug.ge.2) .and.(isroot)) then
             rnorm = pslange(normtype,n,nrhs, Xr,1,1,descX,rwork)
             write(*,*) 'norm(Xr correction) ', rnorm

             rnorm = pdlange(normtype,n,nrhs, R,1,1,descX,dwork)
             write(*,*) 'norm( correction ) ',  rnorm
             write(*,*) 'R(1:5) ', R(1:5)
             write(*,*) 'Xr(1:5) ', Xr(1:5)
             write(*,*) 'X(1:5) ', X(1:5)
           endif

           alpha =   one
           beta = one
           call pdgeadd( 'N', n,nrhs, alpha, R,1,1,descX,                  &
     &                    beta, X,1,1,descX )

*          ----------------
*          compute residual
*          ----------------
           call pdgecopy( n,nrhs, B,ib,jb,descB, R, 1,1, descX )
           alpha = -one
           beta = one
#ifdef USE_PROF
           call profstart('pzgesvr:pzgemm')
#endif
           call pdgemm('N','N',n,nrhs,n,                                    &
     &          alpha, A,ia,ja,descA, X,1,1,descX,                          &
     &          beta, R,1,1,descX )
#ifdef USE_PROF
           call profend('pzgesvr:pzgemm')
#endif
         rnorm_old = rnorm
        enddo

*       ---------------------
*       finally copy solution
*       ---------------------
        call pdgecopy( n,nrhs, X,1,1,descX, B,ib,jb,descB )
        info = 0
        if (.not.isdone) then
                info = -(n+1)
        endif

 999    continue

        deallocate( Ar, X, R, Xr, dwork,rwork,stat=ierr)
        if (ierr.ne.0) then
         write(*,*) 'deallocate(Ar) return ierr=',ierr
        endif

        return

        end subroutine pdgesvr

        subroutine pdgesvr_old( n, nrhs, A, ia, ja, descA,                     &
     &          ipiv, B, ib, jb, descB, info )
        integer, intent(in) :: n,nrhs, ia,ja, ib,jb
        integer, intent(in) :: descA(DLEN_), descB(DLEN_)
        integer, intent(inout) :: info

        integer, intent(inout) :: ipiv(*)
        real(kind=dp), intent(inout) :: A(*)
        real(kind=dp), intent(inout) :: B(*)
*       ---------------
*       local variables
*       ---------------
        logical :: isok, isdone, isroot

        integer :: descAr(DLEN_), descX(DLEN_)
        integer :: ictxt, myrow,mycol, nprow,npcol
        integer :: Locp, Locq, lrindx, lcindx, lld, rsrc,csrc
        integer :: Locp_B, Locq_B


        integer :: iter
        integer, parameter :: maxiter = 30

        integer :: Xsize, Asize, ierr
        real(kind=sp), dimension(:), allocatable :: Ar, Xr
        real(kind=dp), dimension(:), allocatable :: X, R
        real(kind=dp)  :: anorm, eps, bwdcte, xnorm
        real(kind=dp)  :: rnorm,rnorm_old,bwdmax

        real(kind=sp), dimension(:), allocatable :: rwork
        real(kind=dp), dimension(:), allocatable :: dwork
        integer :: ldwork 
        integer :: m,Mp0,Np0

        real(kind=dp), parameter :: one = 1.0d0
        real(kind=dp) :: alpha,beta

*       ----------------------------------------------
*       norm = 'F','f','E', or 'e' for Frobenius norm
*       norm = 'M','m' for max norm
*       norm = '1','O','o' for 1-norm
*       norm = 'I','i', for infinity norm
*       ----------------------------------------------
        character :: normtype = 'I'   

        ictxt = descA(CTXT_)
        call blacs_gridinfo( ictxt, nprow,npcol,myrow,mycol )
        isroot = (myrow.eq.0).and.(mycol.eq.0)

        m = n
        Mp0 = iceil( m, descA(MB_) )*descA(MB_)
        Np0 = iceil( n, descA(NB_) )*descA(NB_)
        ldwork = max(Mp0,Np0) + 1
        allocate( dwork(ldwork),rwork(ldwork),stat=ierr)
        if (ierr.ne.0) then
           write(*,*) 'allocate(dwork) return ierr=',ierr
        endif

         anorm = pdlange(normtype, n,n, A,ia,ja,descA, dwork )
         eps = dlamch('P')
         bwdmax = 1.0d0
*        bwdcte = anorm * eps * min(4.0d0,sqrt(dble(n))/6.0d0)
*        bwdcte = anorm * eps * max(4.0d0,sqrt(dble(n))/6.0d0)
         bwdcte = anorm * eps * sqrt( dble(n) ) * bwdmax

      

*       ----------------------------------
*       make a copy with reduced precision
*       ----------------------------------

*       ----------------------------------------------------
*       set rsrc,csrc so data is exactly aligned
*       on the same processor, no communication is required
*       ----------------------------------------------------
        call infog2l( ia,ja, descA, nprow,npcol, myrow,mycol,               &
     &         lrindx, lcindx, rsrc, csrc )

*       ---------------------------------------
*       for simplicity, use (rsrc,csrc) = (0,0)
*       ---------------------------------------
        rsrc = 0
        csrc = 0

        Locp = numroc( n, descA(MB_), myrow, rsrc, nprow)
        Locq = numroc( n, descA(NB_), mycol, csrc, npcol)
        Asize = max(1,Locp*Locq)

        allocate( Ar( Asize ), stat=ierr)
        if (ierr.ne.0) then
          write(*,*) 'allocate(Ar) return ierr=',ierr
        endif
     

        lld = max(1,Locp)
        call descinit( descAr, n,n, descA(MB_), descA(NB_),                 &
     &       rsrc, csrc, descA(CTXT_), lld, info )
        isok = (info.eq.0)
        if (.not.isok) then
                write(*,*) 'decinit for Ar return info=',info
                stop '** error in pdgesvr '
        endif

        call pconvertd2s( n,n, A,ia,ja,descA, Ar,1,1,descAr )



*       ------------------
*       make a copy of rhs
*       ------------------
        call infog2l( ib,jb, descB, nprow,npcol, myrow,mycol,              &
     &          lrindx, lcindx, rsrc, csrc )
        Locp_B = numroc( n, descB(MB_), myrow, rsrc, nprow)
        Locq_B = numroc( nrhs, descB(NB_), mycol, csrc, npcol)
        Xsize = max(1, Locp_B*Locq_B)
        allocate( X( Xsize ), R(Xsize), Xr(Xsize), stat=ierr )
        if (ierr.ne.0) then
          write(*,*) 'allocate(X) return ierr=',ierr
        endif

*       ----------------
*       setup descriptor
*       ----------------
        lld = max(1, Locp_B )
        call descinit( descX, n, nrhs,  descB(MB_), descB(NB_),            &
     &            rsrc, csrc, descB(CTXT_), lld, info )


        isok = (info.eq.0)
        if (.not.isok) then
                write(*,*) 'descinit for X return info=',info
                stop '** error in pdgesvr '
        endif

        call pdgecopy( n,nrhs, B,ib,jb,descB, X,1,1,descX )

        call pconvertd2s( n,nrhs, X,1,1,descX, Xr,1,1,descX )

        if ((idebug.ge.2) .and. isroot) then
*       ------------
*       double check
*       ------------
        rnorm = pslange(normtype,n,nrhs,Xr,1,1,descX,rwork)
        write(*,*) 'norm(Xr) after c2z ', rnorm


        rnorm = pdlange(normtype, n,nrhs, B,ib,jb,descB, dwork )
        write(*,*) 'norm(B) ', rnorm
        call pconverts2d( n,nrhs, Xr,1,1,descX, R,1,1,descX )


        rnorm = pdlange(normtype,n,nrhs,R,1,1,descX,dwork)
        write(*,*) 'norm(R) after c2z ', rnorm

        alpha = -1.0d0
        beta = 1.0d0
        call pdgeadd( 'N', n,nrhs, alpha, X,1,1,descX,                     &
     &                 beta,R,1,1,descX)
        rnorm = pdlange( normtype, n,nrhs, R,1,1,descX,dwork)
        write(*,*) 'rnorm from conversion ',rnorm

        endif
*       ---------------------------------------------
*       perform LU factorization in reduced precision
*       ---------------------------------------------
#ifdef USE_PROF
        call profstart('pcgetrf')
#endif
        call psgetrf( n,n, Ar, 1,1, descAr, ipiv, info )
#ifdef USE_PROF
        call profend('pcgetrf')
#endif

        isok = (info.eq.0)
        if (.not.isok) then
*          ---------------------------------
*          system may be too ill conditioned
*          ---------------------------------
           goto 999
        endif
*       -----------------------
*       generate  initial guess
*       -----------------------
        call pconvertd2s( n,nrhs, B,ib,jb,descB, Xr,1,1,descX )

#ifdef USE_PROF
        call profstart('psgetrs')
#endif
        call psgetrs( 'N', n, nrhs, Ar, 1,1, descAr, ipiv,                 &
     &        Xr, 1,1, descX, info )
#ifdef USE_PROF
        call profend('psgetrs')
#endif


        call pconverts2d( n,nrhs, Xr,1,1,descX,X,1,1,descX)

        if ((idebug.ge.2).and.isroot) then
                write(*,*) 'initial B(1:5) ', B(1:5)
                write(*,*) 'initial Xr(1:5) ', Xr(1:5)
                write(*,*) 'initial X(1:5) ', X(1:5)
        endif

*       ----------------        
*       compute residual
*       ----------------        
        call pdgecopy( n,nrhs, B, ib,jb, descB, R, 1,1,descX )

        alpha = -one
        beta = one

#ifdef USE_PROF
        call profstart('pzgesvr:pdgemm')
#endif
        call pdgemm('N','N',n,nrhs,n, alpha, A,ia,ja,descA,                 &
     &        X,1,1,descX,  beta, R, 1,1, descX )
#ifdef USE_PROF
        call profend('pzgesvr:pdgemm')
#endif

        xnorm = pdlange(normtype,n,nrhs,B,ib,jb,descB,dwork)
        rnorm_old = xnorm
        do iter=1,maxiter
*          ---------------------
*          check for convergence
*          ---------------------
           rnorm =  pdlange( normtype, n,nrhs, R,1,1,descX, dwork )
           isdone = (rnorm.le.xnorm*bwdcte) .and. (iter.ge.6)

           if ((idebug.ge.1).and.isroot) then
             write(*,*) 'iter,rnorm,rnorm_old ',iter,rnorm,rnorm_old
           endif

           if (isdone) then
                   exit
           endif

           if ((iter.ge.maxiter).and.(rnorm.gt.rnorm_old)) then
              exit
           endif
*          ------------------
*          compute correction
*          ------------------
           call pconvertd2s( n,nrhs, R,1,1, descX, Xr,1,1,descX )
#ifdef USE_PROF
           call profstart('psgetrs')
#endif
           call psgetrs( 'N', n,nrhs, Ar, 1,1, descAr, ipiv,               &
     &              Xr,1,1,descX,info )
#ifdef USE_PROF
           call profend('psgetrs')
#endif
*          ----------------
*          update solution
*          ----------------
           call pconverts2d( n,nrhs, Xr, 1,1,descX, R,1,1,descX )
           if ((idebug.ge.2) .and.(isroot)) then
             rnorm = pslange(normtype,n,nrhs, Xr,1,1,descX,rwork)
             write(*,*) 'norm(Xr correction) ', rnorm

             rnorm = pdlange(normtype,n,nrhs, R,1,1,descX,dwork)
             write(*,*) 'norm( correction ) ',  rnorm
             write(*,*) 'R(1:5) ', R(1:5)
             write(*,*) 'Xr(1:5) ', Xr(1:5)
             write(*,*) 'X(1:5) ', X(1:5)
           endif

           alpha =   one
           beta = one
           call pdgeadd( 'N', n,nrhs, alpha, R,1,1,descX,                  &
     &                    beta, X,1,1,descX )

*          ----------------
*          compute residual
*          ----------------
           call pdgecopy( n,nrhs, B,ib,jb,descB, R, 1,1, descX )
           alpha = -one
           beta = one
#ifdef USE_PROF
           call profstart('pzgesvr:pdgemm')
#endif
           call pdgemm('N','N',n,nrhs,n,                                    &
     &          alpha, A,ia,ja,descA, X,1,1,descX,                          &
     &          beta, R,1,1,descX )
#ifdef USE_PROF
           call profend('pzgesvr:pdgemm')
#endif
         rnorm_old = rnorm
        enddo

*       ---------------------
*       finally copy solution
*       ---------------------
        call pdgecopy( n,nrhs, X,1,1,descX, B,ib,jb,descB )
        info = 0
        if (.not.isdone) then
                info = -(n+1)
        endif

 999    continue

        deallocate( Ar, X, R, Xr, dwork,rwork,stat=ierr)
        if (ierr.ne.0) then
         write(*,*) 'deallocate(Ar) return ierr=',ierr
        endif

        return

        end subroutine pdgesvr_old

        end module pgesvr_mod
