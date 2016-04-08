      subroutine  iafirst( ia,mb,myprow,rsrc,nprow,iafirst_div)
      integer ia,mb,myprow,rsrc,nprow
      integer iafirst_div
      implicit none

      integer iia,iproc0,iproc,iproc_shift,myprow_shift
      logical is_found

      iia = ia
      iproc0 = indxg2p(iia,mb,myprow,rsrc,nprow)
      iproc = iproc0
      is_found = (iproc0 .eq. myprow)
      if (is_found) then
          iafirst_div = iia
          return
      endif


      iproc_shift = iproc - rsrc
      do while (iproc_shift .lt. 0) 
         iproc_shift = iproc_shift + nprow
      enddo

      myprow_shift = myprow - rsrc
      do while (myprow_shift .lt. iproc_shift)
        myprow_shift = myprow_shift + nprow
      enddo

      iia  = ia + (mb - MOD( (ia-1),mb)) +                               &
     &            ((myprow_shift-1)-iproc_shift)*mb
      iafirst_div = iia
      return
      end subroutine iafirst

      subroutine ialast( ia,mb,myprow,rsrc,nprow, iia_ans )
      integer ia,mb,myprow,rsrc,nprow
      integer iia_ans
      implicit none

      iia = -1
      do i=0,(nprow*mb-1)
         iia = ia - i
         iproc = indxg2p(iia,mb,myprow,rsrc,nprow)
         is_found = (iproc .eq. myprow)
         if (is_found) then
             exit
          endif
       enddo
       iia_ans = iia
       return
       end subroutine ialast



