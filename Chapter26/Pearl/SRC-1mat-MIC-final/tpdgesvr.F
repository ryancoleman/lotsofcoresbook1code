        program tpdgesvr
        use scalapack_mod
        use crandom_mod
        ! use prof_mod
        use pgesvr_mod

        implicit none
*        include 'desc.h'
        integer idebug
        parameter(idebug=1)
*       ------------------------------------------------
*       simple program to test pdgetrf versus variants of
*       the algorithm
*       ------------------------------------------------
        integer descA(DLEN_), descB(DLEN_)

        character*80 ctitle

        logical isroot,isok 
        logical randomize_seed, well_conditioned
        integer mb,nb,m,n,nrhs,rsrc,csrc
        integer Locp,Locq,lld, Asize,Bsize,ipivsize

        integer, parameter :: dp = kind(1.0d0)
        integer, parameter :: sp = kind(1.0)

        real(kind=dp), allocatable, dimension(:) :: A
        real(kind=dp), allocatable, dimension(:) :: X,B,B2,Borg
        integer isize,ntimes,i
        integer, allocatable, dimension(:) :: seed_array
        real(kind=dp), allocatable, dimension(:) :: rseed_array

        integer, allocatable, dimension(:) :: ipiv

        integer nprow,npcol,myprow,mypcol,info,icontext
        integer ia,ja,ib,jb
        real(kind=dp) alpha, beta
        real(kind=dp) errB
        real(kind=dp) t1, t2

        integer  irank, nproc
        integer ierr
        integer indev
        parameter(indev=15)
        logical, parameter :: use_input = .true.

#ifdef USE_HPL
        real(kind=sp), allocatable, dimension(:) :: Ar
#endif

       

#ifdef USE_GPU
       call cublasInit()
#endif

        ctitle = ' '

        well_conditioned = .false.
        randomize_seed = .false.

#ifdef USE_PROF
        call profinit()
#endif

*       ----------------
*       setup blacs grid
*       ----------------
        call blacs_pinfo(irank,nproc)
        call blacs_get(-1,0,icontext)
        do nprow=max(1,int(sqrt(dble(nproc)))),1,-1
          npcol = nproc/nprow
          if (nprow*npcol.eq.nproc) exit
        enddo
        call blacs_gridinit( icontext, 'Row-major',nprow,npcol)
        call blacs_gridinfo( icontext, nprow,npcol, myprow,mypcol)
        isroot = (myprow.eq.0).and.(mypcol.eq.0)


         
        if (use_input) then
         if (isroot) then
            open(indev,file='tpgesvr.dat',form='formatted',                 &
     &        status='old',iostat=info)
 
           rewind(indev)
           read(indev,*) m,nrhs 
           close(indev)
         endif

*          -----------------------------------
*          broadcast input to other processors
*          -----------------------------------
         if (isroot) then
             call igebs2d(icontext,'A',' ',1,1,m,1)
             call igebs2d(icontext,'A',' ',1,1,nrhs,1)
         else
             call igebr2d(icontext,'A',' ',1,1,m,1, 0,0)
             call igebr2d(icontext,'A',' ',1,1,nrhs,1, 0,0)
         endif

         else
          m=500
          nrhs=1
         endif






        if (isroot) then
        write(*,*) 'testing ',trim(ctitle)
        write(*,*) 'm,nrhs',m,nrhs
        write(*,*) 'nprow,npcol ',nprow,npcol
        endif

*       ---------------------
*       create matrices, A
*       ---------------------
        n = m

        mb = 99
        nb = 99
        rsrc = 0
        csrc = 0

        Locp = numroc( m, mb, myprow, rsrc, nprow )
        Locq = numroc( n, nb, mypcol, csrc, npcol )
        lld = max(1,Locp)
        Asize = max(1,lld * Locq)

#ifdef USE_HPL
        call blacs_barrier(icontext,'A')
        call hpl_sblacsinit( icontext )
#endif


        call descinit( descA, m,n, mb,nb,rsrc,csrc,icontext,lld,info)
        call assert(info.eq.0,                                                &
     &             'main:descinit(descA) return info = ',info)


        ipivsize = m + mb

        allocate( A(Asize), ipiv(ipivsize), stat=ierr)
        if (ierr.ne.0) then
          write(*,*) 'allocate(A(Asize)) return ierr=',ierr
          write(*,*) 'mb,nb ', mb,nb
          write(*,*) 'myprow,mypcol ', myprow,mypcol
          write(*,*) 'Locp, Locq ', Locp,Locq
        endif

        ipiv(1:ipivsize) = -1
        A(1:Asize) = -9.9e9

*       ----------
*       rhs vector
*       ----------
        Locp = numroc( m, mb, myprow, rsrc, nprow )
        Locq = numroc( nrhs, nb, mypcol, csrc, npcol )
        lld = max(1,Locp)

        call descinit( descB, m,nrhs, mb,nb,rsrc,csrc,icontext,lld,info)
        call assert(info.eq.0,                                              &
     &              'main:descinit(descB) return info = ',info)

        Bsize = max(1, lld*Locq )
        allocate( B2(Bsize),Borg(Bsize), B(Bsize), X(Bsize), stat=ierr )
        if (ierr.ne.0) then
          write(*,*) 'allocate( B(Bsize) ), return ierr=',ierr
        endif

        call crandom_number(Bsize, X)
        B(1:Bsize) = -9.9e9
        B2(1:Bsize) = -9.9e9
        Borg(1:Bsize) = -9.9e9




*       --------------------
*       generate random data
*       --------------------
        call blacs_barrier(icontext,'All')
#ifdef USE_PROF
        call blacs_barrier(icontext,'All')
        call profstart('random')
#endif
        ntimes = 1+(myprow + mypcol*nprow)

        if (randomize_seed) then
        call random_seed(size=isize)
        allocate( seed_array(isize), rseed_array(isize) )
        do i=1,ntimes
          call random_number(rseed_array)
        enddo
        seed_array = int(rseed_array*(2**28))
        call random_seed(put=seed_array)
        else
         do i=1,ntimes
           call crandom_number(Asize,A)
         enddo
        endif

        call crandom_number(Asize,A)

*         call crandom_number(Bsize,B)

#ifdef USE_PROF
        call blacs_barrier(icontext,'All')
        call profend('random')
#endif


        if (well_conditioned) then
#ifdef USE_PROF
                call profstart('well_conditioned')
#endif
                A(1:Asize) = -1.0d0
                do ia=1,min(m,n)
                   alpha = n+1
                   ja = ia
                   call pdelset( A,ia,ja,descA, alpha )
                enddo
#ifdef USE_PROF
                call profend('well_conditioned')
#endif
        endif

        beta = 0.0d0
        alpha = 1.0d0
        call pdgemm( 'N', 'N', m,nrhs,m, alpha, A,1,1,descA,                &
     &            X,1,1,descB, beta, B,1,1,descB )
        Borg(1:Bsize) = B(1:Bsize)
        B2(1:Bsize) = B(1:Bsize)


#ifdef USE_HPL
*       ----------------------------------
*       time solution in reduced precision
*       ----------------------------------
        allocate( Ar(Asize), stat=info)
        Ar = 0.0
        call assert(info.eq.0,'allocate Ar failed',info)
        call pconvertd2s( n,n,A,1,1,descA, Ar,1,1,descA )

        call blacs_barrier(icontext,'All')
        call cpu_time(t1)
#ifdef USE_PROF
        call profstart('tpdgesvr:psgetrf')
#endif

        call psgetrf(n,n,Ar,1,1,descA,ipiv,info)

        call blacs_barrier(icontext,'All')
        call cpu_time(t2)
#ifdef USE_PROF
        call profend('tpdgesvr:psgetrf')
#endif
        if (isroot) then
           write(*,*) 'tpdgesvr, psgetrf took ', t2-t1,' sec '
        endif

        deallocate( Ar, stat=info)
        call assert(info.eq.0,'deallocate Ar failed',info)

#endif


        
*       ----------
*       use pgesvr
*       ----------
        ia = 1
        ja = 1
        ib = 1
        jb = 1
        info = 0
        call blacs_barrier(icontext,'All')
        call cpu_time(t1)
#ifdef USE_PROF
        call profstart('pdgesvr')
#endif

        call pdgesvr( n,nrhs,A,ia,ja,descA,ipiv,B2,ib,jb,descB,info)

        call blacs_barrier(icontext,'All')
        call cpu_time(t2)
#ifdef USE_PROF
        call profend('pdgesvr')
#endif
        if (isroot) then 
          write(*,*) 'pdgesvr took  ', t2-t1,' sec'
        endif


        isok = (info.eq.0)
        call assert(isok,'pgesvr return info=',info)

*       ---------------------
*       use scalapack pdgetrf
*       ---------------------
        ia = 1
        ja = 1
        info = 0
        B(1:Bsize) = Borg(1:Bsize)

        call blacs_barrier(icontext,'All')
        call cpu_time(t1)
#ifdef USE_PROF
        call profstart('pdgetrf')
#endif

        call pdgetrf(m,n,A,ia,ja,descA,ipiv,info)

        call blacs_barrier(icontext,'All')
        call cpu_time(t2)
#ifdef USE_PROF
        call profend('pdgetrf')
#endif
        if (isroot) then
          write(*,*) 'pdgetrf took ', t2-t1,' sec'
        endif

        isok = (info.eq.0)
        call assert(isok,'pdgetrf return info=',info)

        ib = 1
        jb = 1

#ifdef USE_PROF
        call blacs_barrier(icontext,'A')
        call profstart('pdgetrs')
#endif
        call pdgetrs('NoTrans',m,nrhs,A,ia,ja,descA,ipiv,
     &              B,ib,jb,descB, info )
#ifdef USE_PROF
        call blacs_barrier(icontext,'A')
        call profend('pdgetrs')
#endif


        isok = (info.eq.0)
        call assert(isok,'pdgetrs return info=',info)



*       ------------
*       check result
*       ------------
        errB = maxval( abs(B(1:Bsize)-B2(1:Bsize)) )

        call dgamx2d(icontext,'All',' ',1,1,errB,1,-1,-1,-1,-1,-1)
        if (isroot) then
                write(*,*) 'difference in B is ',errB
        endif



        if (isroot) then
#ifdef USE_PROF
                call profstat()
                call profinit()
#endif
        endif


!       ---------------
!       prepare to exit
!       ---------------
        call blacs_barrier(icontext, 'All')
        call blacs_gridexit( icontext )
        call blacs_exit(0)
        stop
        end
