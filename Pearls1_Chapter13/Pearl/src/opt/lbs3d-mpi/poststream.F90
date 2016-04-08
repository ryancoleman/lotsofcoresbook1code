!-------------------------------------------------------------------------------
! Subroutine : PostStream
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Update f in nodes at processor boundaries.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Copyright 2013 Carlos Rosales Fernandez and The University of Texas at Austin.
! Copyright 2008 Carlos Rosales Fernandez, David S. Whyte and IHPC (A*STAR).
!
! This file is part of MP-LABS.
!
! MP-LABS is free software: you can redistribute it and/or modify it under the
! terms of the GNU GPL version 3 or (at your option) any later version.
!
! MP-LABS is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! MP-LABS, in the file COPYING.txt. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------

 SUBROUTINE PostStream

!  Common Variables
 USE Domain
 USE FluidParams, ONLY : phi
 USE LBMParams
 USE MPIParams
 USE MPI
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k
 INTEGER :: MPI_ERR, MPI_REQ_X(4), MPI_REQ_Y(4), MPI_REQ_Z(4)
 INTEGER :: MPI_REQ_EX(8), MPI_REQ_EY(8), MPI_REQ_EZ(8)
 INTEGER :: MPI_STAT(MPI_STATUS_SIZE,4), MPI_ESTAT(MPI_STATUS_SIZE,8)

!$OMP PARALLEL

!----------Exchange inward pointing f values ---------------------------------
!$OMP SINGLE
 CALL MPI_IRECV( f_east_rcv,  xsize, MPI_DOUBLE_PRECISION, east,  TAG1, &
                 MPI_COMM_VGRID, MPI_REQ_X(1), MPI_ERR )
 CALL MPI_IRECV( f_west_rcv,  xsize, MPI_DOUBLE_PRECISION, west,  TAG2, &
                 MPI_COMM_VGRID, MPI_REQ_X(2), MPI_ERR )
 CALL MPI_IRECV( f_north_rcv, ysize, MPI_DOUBLE_PRECISION, north, TAG3, &
                 MPI_COMM_VGRID, MPI_REQ_Y(1), MPI_ERR )
 CALL MPI_IRECV( f_south_rcv, ysize, MPI_DOUBLE_PRECISION, south, TAG4, &
                 MPI_COMM_VGRID, MPI_REQ_Y(2), MPI_ERR )
 CALL MPI_IRECV( f_top_rcv,   zsize, MPI_DOUBLE_PRECISION, top,   TAG5, &
                 MPI_COMM_VGRID, MPI_REQ_Z(1), MPI_ERR )
 CALL MPI_IRECV( f_bot_rcv,   zsize, MPI_DOUBLE_PRECISION, bot,   TAG6, &
                 MPI_COMM_VGRID, MPI_REQ_Z(2), MPI_ERR )
!$OMP END SINGLE

! X direction
!$OMP DO PRIVATE(j)
 DO k = 1, NZ
   DO j = 1, NY
     f_west_snd(j+NY*(k-1)) = f2( NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     f_east_snd(j+NY*(k-1)) = f1( NX+1 + NYG*( j + NZG*k ) + nxt )
   END DO
 END DO
!$OMP SINGLE
 CALL MPI_ISEND( f_west_snd, xsize , MPI_DOUBLE_PRECISION, west, TAG1, &
                 MPI_COMM_VGRID, MPI_REQ_X(3), MPI_ERR )
 CALL MPI_ISEND( f_east_snd, xsize , MPI_DOUBLE_PRECISION, east, TAG2, &
                 MPI_COMM_VGRID, MPI_REQ_X(4), MPI_ERR )
!$OMP END SINGLE

! Y direction
!$OMP DO PRIVATE(i)
 DO k = 1, NZ
   DO i = 1, NX
     f_south_snd(i+NX*(k-1)) = f4( i + NXG*( NYG*k ) + nxt )
   END DO
   DO i = 1, NX
     f_north_snd(i+NX*(k-1)) = f3( i + NXG*( NY+1 + NYG*k ) + nxt )
   END DO
 END DO
!$OMP SINGLE
 CALL MPI_ISEND( f_south_snd, ysize, MPI_DOUBLE_PRECISION, south, TAG3, &
                 MPI_COMM_VGRID, MPI_REQ_Y(3), MPI_ERR )
 CALL MPI_ISEND( f_north_snd, ysize, MPI_DOUBLE_PRECISION, north, TAG4, &
                 MPI_COMM_VGRID, MPI_REQ_Y(4), MPI_ERR )
!$OMP END SINGLE

! Z direction
!$OMP DO PRIVATE(i)
 DO j = 1, NY
   DO i = 1, NX
     f_bot_snd(i+NX*(j-1)) = f6( i + NXG*( j ) + nxt )
   END DO
   DO i = 1, NX
     f_top_snd(i+NX*(j-1)) = f5( i + NXG*( j + NYG*(NZ+1) ) + nxt )
   END DO
 END DO
!$OMP SINGLE
 CALL MPI_ISEND( f_bot_snd, zsize, MPI_DOUBLE_PRECISION, bot, TAG5, &
                 MPI_COMM_VGRID, MPI_REQ_Z(3), MPI_ERR )
 CALL MPI_ISEND( f_top_snd, zsize, MPI_DOUBLE_PRECISION, top, TAG6, &
                 MPI_COMM_VGRID, MPI_REQ_Z(4), MPI_ERR )

!-------- Update inward values of f in real nodes ----------------------------
! X direction
 CALL MPI_WAITALL(4, MPI_REQ_X, MPI_STAT, MPI_ERR)
!$OMP END SINGLE
!$OMP BARRIER
!$OMP DO PRIVATE(j)
 DO k = 1, NZ
   DO j = 1, NY
     f2( NX + NXG*( j + NYG*k ) + nxt ) = f_east_rcv(j+NY*(k-1))
   END DO
   DO j = 1, NY
     f1( 1 + NXG*( j + NYG*k ) + nxt ) = f_west_rcv(j+NY*(k-1))
   END DO
 END DO

! Y direction
!$OMP SINGLE
 CALL MPI_WAITALL(4, MPI_REQ_Y, MPI_STAT, MPI_ERR)
!$OMP END SINGLE
!$OMP BARRIER
!$OMP DO PRIVATE(i)
 DO k = 1, NZ
   DO i = 1, NX
     f4( i + NXG*( NY + NYG*k ) + nxt ) = f_north_rcv(i+NX*(k-1))
   END DO
   DO i = 1, NX
     f3( i + NXG*( 1 + NYG*k ) + nxt ) = f_south_rcv(i+NX*(k-1))
   END DO
 END DO

! Z direction
!$OMP SINGLE
 CALL MPI_WAITALL(4, MPI_REQ_Z, MPI_STAT, MPI_ERR)
!$OMP END SINGLE
!$OMP BARRIER
!$OMP DO PRIVATE(i)
 DO j = 1, NY
   DO i = 1, NX
     f6( i + NXG*( j + NYG*NZ ) + nxt ) = f_top_rcv(i+NX*(j-1))
   END DO
   DO i = 1, NX
     f5( i + NXG*( j + NYG ) + nxt ) = f_bot_rcv(i+NX*(j-1))
   END DO
 END DO

!
! At this stage all values of f needed for the next step have been updated
!


!-------- Exchange data for phi and inward g components -----------------------
!$OMP SINGLE
 CALL MPI_IRECV( east_rcv,  xsize6, MPI_DOUBLE_PRECISION, east,  TAG1, &
                 MPI_COMM_VGRID, MPI_REQ_X(1), MPI_ERR )
 CALL MPI_IRECV( west_rcv,  xsize6, MPI_DOUBLE_PRECISION, west,  TAG2, &
                 MPI_COMM_VGRID, MPI_REQ_X(2), MPI_ERR )
 CALL MPI_IRECV( north_rcv, ysize6, MPI_DOUBLE_PRECISION, north, TAG3, &
                 MPI_COMM_VGRID, MPI_REQ_Y(1), MPI_ERR )
 CALL MPI_IRECV( south_rcv, ysize6, MPI_DOUBLE_PRECISION, south, TAG4, &
                 MPI_COMM_VGRID, MPI_REQ_Y(2), MPI_ERR )
 CALL MPI_IRECV( top_rcv,   zsize6, MPI_DOUBLE_PRECISION, top,   TAG5, &
                 MPI_COMM_VGRID, MPI_REQ_Z(1), MPI_ERR )
 CALL MPI_IRECV( bot_rcv,   zsize6, MPI_DOUBLE_PRECISION, bot,   TAG6, &
                 MPI_COMM_VGRID, MPI_REQ_Z(2), MPI_ERR )
!$OMP END SINGLE


! X direction
!$OMP DO PRIVATE(j)
 DO k = 1, NZ
   DO j = 1, NY
     west_snd(j+NY*(k-1)) = f0( 1 + NXG*( j + NYG*k )       ) &
                          + f1( 1 + NXG*( j + NYG*k ) + nxt ) &
                          + f2( 1 + NXG*( j + NYG*k ) + nxt ) &
                          + f3( 1 + NXG*( j + NYG*k ) + nxt ) &
                          + f4( 1 + NXG*( j + NYG*k ) + nxt ) &
                          + f5( 1 + NXG*( j + NYG*k ) + nxt ) &
                          + f6( 1 + NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     west_snd(j+NY*(k-1+NZ)) = g2( NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     west_snd(j+NY*(k-1+2*NZ)) = g8( NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     west_snd(j+NY*(k-1+3*NZ)) = g10( NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     west_snd(j+NY*(k-1+4*NZ)) = g12( NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     west_snd(j+NY*(k-1+5*NZ)) = g14( NXG*( j + NYG*k ) + nxt )
   END DO

   DO j = 1, NY
     east_snd(j+NY*(k-1)) = f0( NX + NXG*( j + NYG*k )       ) &
                          + f1( NX + NXG*( j + NYG*k ) + nxt ) &
                          + f2( NX + NXG*( j + NYG*k ) + nxt ) &
                          + f3( NX + NXG*( j + NYG*k ) + nxt ) &
                          + f4( NX + NXG*( j + NYG*k ) + nxt ) &
                          + f5( NX + NXG*( j + NYG*k ) + nxt ) &
                          + f6( NX + NXG*( j + NYG*k ) + nxt )  
   END DO
   DO j = 1, NY
     east_snd(j+NY*(k-1+NZ)) = g1( NX+1 + NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     east_snd(j+NY*(k-1+2*NZ)) = g7( NX+1 + NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     east_snd(j+NY*(k-1+3*NZ)) = g9( NX+1 + NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     east_snd(j+NY*(k-1+4*NZ)) = g11( NX+1 + NXG*( j + NYG*k ) + nxt )
   END DO
   DO j = 1, NY
     east_snd(j+NY*(k-1+5*NZ)) = g13( NX+1 + NXG*( j + NYG*k ) + nxt )
   END DO
 END DO
!$OMP SINGLE
 CALL MPI_ISEND( west_snd, xsize6, MPI_DOUBLE_PRECISION, west, TAG1, &
                 MPI_COMM_VGRID, MPI_REQ_X(3), MPI_ERR )
 CALL MPI_ISEND( east_snd, xsize6, MPI_DOUBLE_PRECISION, east, TAG2, &
                 MPI_COMM_VGRID, MPI_REQ_X(4), MPI_ERR )
!$OMP END SINGLE

! Y direction
!$OMP DO PRIVATE(i)
 DO k = 1, NZ
   DO i = 1, NX
     south_snd(i+NX*(k-1)) = f0( i + NXG*( 1 + NYG*k )       ) &
                           + f1( i + NXG*( 1 + NYG*k ) + nxt ) &
                           + f2( i + NXG*( 1 + NYG*k ) + nxt ) &
                           + f3( i + NXG*( 1 + NYG*k ) + nxt ) &
                           + f4( i + NXG*( 1 + NYG*k ) + nxt ) &
                           + f5( i + NXG*( 1 + NYG*k ) + nxt ) &
                           + f6( i + NXG*( 1 + NYG*k ) + nxt )
   END DO
   DO i = 1, NX
     south_snd(i+NX*(k-1+NZ)) = g4( i + NXG*( NYG*k) + nxt )
   END DO
   DO i = 1, NX
     south_snd(i+NX*(k-1+2*NZ)) = g8( i + NXG*( NYG*k) + nxt )
   END DO
   DO i = 1, NX
     south_snd(i+NX*(k-1+3*NZ)) = g9( i + NXG*( NYG*k) + nxt )
   END DO
   DO i = 1, NX
     south_snd(i+NX*(k-1+4*NZ)) = g16( i + NXG*( NYG*k) + nxt )
   END DO
   DO i = 1, NX
     south_snd(i+NX*(k-1+5*NZ)) = g18( i + NXG*( NYG*k) + nxt )
   END DO

   DO i = 1, NX
     north_snd(i+NX*(k-1)) = f0( i + NXG*( NY + NYG*k )       ) &
                           + f1( i + NXG*( NY + NYG*k ) + nxt ) &
                           + f2( i + NXG*( NY + NYG*k ) + nxt ) &
                           + f3( i + NXG*( NY + NYG*k ) + nxt ) &
                           + f4( i + NXG*( NY + NYG*k ) + nxt ) &
                           + f5( i + NXG*( NY + NYG*k ) + nxt ) &
                           + f6( i + NXG*( NY + NYG*k ) + nxt )
   END DO
   DO i = 1, NX
     north_snd(i+NX*(k-1+NZ)) = g3( i + NXG*( NY+1 + NYG*k ) + nxt )
   END DO
   DO i = 1, NX
     north_snd(i+NX*(k-1+2*NZ)) = g7( i + NXG*( NY+1 + NYG*k ) + nxt )
   END DO
   DO i = 1, NX
     north_snd(i+NX*(k-1+3*NZ)) = g10( i + NXG*( NY+1 + NYG*k ) + nxt )
   END DO
   DO i = 1, NX
     north_snd(i+NX*(k-1+4*NZ)) = g15( i + NXG*( NY+1 + NYG*k ) + nxt )
   END DO
   DO i = 1, NX
     north_snd(i+NX*(k-1+5*NZ)) = g17( i + NXG*( NY+1 + NYG*k ) + nxt )
   END DO
 END DO
!$OMP SINGLE
 CALL MPI_ISEND( south_snd, ysize6, MPI_DOUBLE_PRECISION, south, TAG3, &
                 MPI_COMM_VGRID, MPI_REQ_Y(3), MPI_ERR )
 CALL MPI_ISEND( north_snd, ysize6, MPI_DOUBLE_PRECISION, north, TAG4, &
                 MPI_COMM_VGRID, MPI_REQ_Y(4), MPI_ERR )
!$OMP END SINGLE

! Z direction
!$OMP DO PRIVATE(i)
 DO j = 1, NY
   DO i = 1, NX
     bot_snd(i+NX*(j-1)) = f0( i + NXG*( j + NYG )       ) &
                         + f1( i + NXG*( j + NYG ) + nxt ) &
                         + f2( i + NXG*( j + NYG ) + nxt ) &
                         + f3( i + NXG*( j + NYG ) + nxt ) &
                         + f4( i + NXG*( j + NYG ) + nxt ) &
                         + f5( i + NXG*( j + NYG ) + nxt ) &
                         + f6( i + NXG*( j + NYG ) + nxt )
   END DO
   DO i = 1, NX
     bot_snd(i+NX*(j-1+NY)) = g6( i + NXG*j + nxt )
   END DO
   DO i = 1, NX
     bot_snd(i+NX*(j-1+2*NY)) = g12( i + NXG*j + nxt )
   END DO
   DO i = 1, NX
     bot_snd(i+NX*(j-1+3*NY)) = g13( i + NXG*j + nxt )
   END DO
   DO i = 1, NX
     bot_snd(i+NX*(j-1+4*NY)) = g16( i + NXG*j + nxt )
   END DO
   DO i = 1, NX
     bot_snd(i+NX*(j-1+5*NY)) = g17( i + NXG*j + nxt )
   END DO

   DO i = 1, NX
     top_snd(i+NX*(j-1)) = f0( i + NXG*( j + NYG*NZ )       ) &
                         + f1( i + NXG*( j + NYG*NZ ) + nxt ) &
                         + f2( i + NXG*( j + NYG*NZ ) + nxt ) &
                         + f3( i + NXG*( j + NYG*NZ ) + nxt ) &
                         + f4( i + NXG*( j + NYG*NZ ) + nxt ) &
                         + f5( i + NXG*( j + NYG*NZ ) + nxt ) &
                         + f6( i + NXG*( j + NYG*NZ ) + nxt )
   END DO
   DO i = 1, NX
     top_snd(i+NX*(j-1+NY)) = g5( i + NXG*( j + NYG*(NZ+1) ) + nxt )
   END DO
   DO i = 1, NX
     top_snd(i+NX*(j-1+2*NY)) = g11( i + NXG*( j + NYG*(NZ+1) ) + nxt )
   END DO
   DO i = 1, NX
     top_snd(i+NX*(j-1+3*NY)) = g14( i + NXG*( j + NYG*(NZ+1) ) + nxt )
   END DO
   DO i = 1, NX
     top_snd(i+NX*(j-1+4*NY)) = g15( i + NXG*( j + NYG*(NZ+1) ) + nxt )
   END DO
   DO i = 1, NX
     top_snd(i+NX*(j-1+5*NY)) = g18( i + NXG*( j + NYG*(NZ+1) ) + nxt )
   END DO
 END DO
!$OMP SINGLE
 CALL MPI_ISEND( bot_snd, zsize6, MPI_DOUBLE_PRECISION, bot, TAG5, &
                 MPI_COMM_VGRID, MPI_REQ_Z(3), MPI_ERR )
 CALL MPI_ISEND( top_snd, zsize6, MPI_DOUBLE_PRECISION, top, TAG6, &
                 MPI_COMM_VGRID, MPI_REQ_Z(4), MPI_ERR )

!----------------- EDGES -------------------------------------------------------
 CALL MPI_IRECV( tn_rcv, xedge2, MPI_DOUBLE_PRECISION, tn, TAG7,  &
                 MPI_COMM_VGRID, MPI_REQ_EX(1), MPI_ERR )
 CALL MPI_IRECV( ts_rcv, xedge2, MPI_DOUBLE_PRECISION, ts, TAG8,  &
                 MPI_COMM_VGRID, MPI_REQ_EX(2), MPI_ERR )
 CALL MPI_IRECV( bn_rcv, xedge2, MPI_DOUBLE_PRECISION, bn, TAG9,  &
                 MPI_COMM_VGRID, MPI_REQ_EX(3), MPI_ERR )
 CALL MPI_IRECV( bs_rcv, xedge2, MPI_DOUBLE_PRECISION, bs, TAG10, &
                 MPI_COMM_VGRID, MPI_REQ_EX(4), MPI_ERR )
 CALL MPI_IRECV( te_rcv, yedge2, MPI_DOUBLE_PRECISION, te, TAG11, &
                 MPI_COMM_VGRID, MPI_REQ_EY(1), MPI_ERR )
 CALL MPI_IRECV( tw_rcv, yedge2, MPI_DOUBLE_PRECISION, tw, TAG12, &
                 MPI_COMM_VGRID, MPI_REQ_EY(2), MPI_ERR )
 CALL MPI_IRECV( be_rcv, yedge2, MPI_DOUBLE_PRECISION, be, TAG13, &
                 MPI_COMM_VGRID, MPI_REQ_EY(3), MPI_ERR )
 CALL MPI_IRECV( bw_rcv, yedge2, MPI_DOUBLE_PRECISION, bw, TAG14, &
                 MPI_COMM_VGRID, MPI_REQ_EY(4), MPI_ERR )
 CALL MPI_IRECV( ne_rcv, zedge2, MPI_DOUBLE_PRECISION, ne, TAG15, &
                 MPI_COMM_VGRID, MPI_REQ_EZ(1), MPI_ERR )
 CALL MPI_IRECV( nw_rcv, zedge2, MPI_DOUBLE_PRECISION, nw, TAG16, &
                 MPI_COMM_VGRID, MPI_REQ_EZ(2), MPI_ERR )
 CALL MPI_IRECV( se_rcv, zedge2, MPI_DOUBLE_PRECISION, se, TAG17, &
                 MPI_COMM_VGRID, MPI_REQ_EZ(3), MPI_ERR )
 CALL MPI_IRECV( sw_rcv, zedge2, MPI_DOUBLE_PRECISION, sw, TAG18, &
                 MPI_COMM_VGRID, MPI_REQ_EZ(4), MPI_ERR )
!$OMP END SINGLE

! Edges along X
!$OMP DO
 DO i = 1, NX
   tn_snd(i)    = f0( i + NXG*( NY + NYG*NZ )       ) &
                + f1( i + NXG*( NY + NYG*NZ ) + nxt ) &
                + f2( i + NXG*( NY + NYG*NZ ) + nxt ) &
                + f3( i + NXG*( NY + NYG*NZ ) + nxt ) &
                + f4( i + NXG*( NY + NYG*NZ ) + nxt ) &
                + f5( i + NXG*( NY + NYG*NZ ) + nxt ) &
                + f6( i + NXG*( NY + NYG*NZ ) + nxt )
   tn_snd(i+NX) = g15( i + NXG*( NY+1 + NYG*(NZ+1) ) + nxt )

   ts_snd(i)    = f0( i + NXG*( 1 + NYG*NZ )       ) &
                + f1( i + NXG*( 1 + NYG*NZ ) + nxt ) &
                + f2( i + NXG*( 1 + NYG*NZ ) + nxt ) &
                + f3( i + NXG*( 1 + NYG*NZ ) + nxt ) &
                + f4( i + NXG*( 1 + NYG*NZ ) + nxt ) &
                + f5( i + NXG*( 1 + NYG*NZ ) + nxt ) &
                + f6( i + NXG*( 1 + NYG*NZ ) + nxt )
   ts_snd(i+NX) = g18( i + NXG*NYG*(NZ+1) + nxt )

   bs_snd(i)    = f0( i + NXG*( 1 + NYG )       ) &
                + f1( i + NXG*( 1 + NYG ) + nxt ) &
                + f2( i + NXG*( 1 + NYG ) + nxt ) &
                + f3( i + NXG*( 1 + NYG ) + nxt ) &
                + f4( i + NXG*( 1 + NYG ) + nxt ) &
                + f5( i + NXG*( 1 + NYG ) + nxt ) &
                + f6( i + NXG*( 1 + NYG ) + nxt )
   bs_snd(i+NX) = g16( i + nxt )

   bn_snd(i)    = f0( i + NXG*( NY + NYG )       ) &
                + f1( i + NXG*( NY + NYG ) + nxt ) &
                + f2( i + NXG*( NY + NYG ) + nxt ) &
                + f3( i + NXG*( NY + NYG ) + nxt ) &
                + f4( i + NXG*( NY + NYG ) + nxt ) &
                + f5( i + NXG*( NY + NYG ) + nxt ) &
                + f6( i + NXG*( NY + NYG ) + nxt )
   bn_snd(i+NX) = g17( i + NXG*( NY+1 ) + nxt )
 END DO
!$OMP SINGLE
 CALL MPI_ISEND( bs_snd, xedge2, MPI_DOUBLE_PRECISION, bs, TAG7,  &
                 MPI_COMM_VGRID, MPI_REQ_EX(5), MPI_ERR )
 CALL MPI_ISEND( bn_snd, xedge2, MPI_DOUBLE_PRECISION, bn, TAG8,  &
                 MPI_COMM_VGRID, MPI_REQ_EX(6), MPI_ERR )
 CALL MPI_ISEND( ts_snd, xedge2, MPI_DOUBLE_PRECISION, ts, TAG9,  &
                 MPI_COMM_VGRID, MPI_REQ_EX(7), MPI_ERR )
 CALL MPI_ISEND( tn_snd, xedge2, MPI_DOUBLE_PRECISION, tn, TAG10, &
                 MPI_COMM_VGRID, MPI_REQ_EX(8), MPI_ERR )
!$OMP END SINGLE

! Edges along Y
!$OMP DO
 DO j = 1, NY
   te_snd(j)    = f0( NX + NXG*( j + NYG*NZ )       ) &
                + f1( NX + NXG*( j + NYG*NZ ) + nxt ) &
                + f2( NX + NXG*( j + NYG*NZ ) + nxt ) &
                + f3( NX + NXG*( j + NYG*NZ ) + nxt ) &
                + f4( NX + NXG*( j + NYG*NZ ) + nxt ) &
                + f5( NX + NXG*( j + NYG*NZ ) + nxt ) &
                + f6( NX + NXG*( j + NYG*NZ ) + nxt )
   te_snd(j+NY) = g11( NX+1 + NXG*( j + NYG*(NZ+1) ) + nxt )

   tw_snd(j)    = f0( 1 + NXG*( j + NYG*NZ )       ) &
                + f1( 1 + NXG*( j + NYG*NZ ) + nxt ) &
                + f2( 1 + NXG*( j + NYG*NZ ) + nxt ) &
                + f3( 1 + NXG*( j + NYG*NZ ) + nxt ) &
                + f4( 1 + NXG*( j + NYG*NZ ) + nxt ) &
                + f5( 1 + NXG*( j + NYG*NZ ) + nxt ) &
                + f6( 1 + NXG*( j + NYG*NZ ) + nxt )
   tw_snd(j+NY) = g14( NXG*( j + NYG*(NZ+1) ) + nxt )

   bw_snd(j)    = f0( 1 + NXG*( j + NYG )       ) &
                + f1( 1 + NXG*( j + NYG ) + nxt ) &
                + f2( 1 + NXG*( j + NYG ) + nxt ) &
                + f3( 1 + NXG*( j + NYG ) + nxt ) &
                + f4( 1 + NXG*( j + NYG ) + nxt ) &
                + f5( 1 + NXG*( j + NYG ) + nxt ) &
                + f6( 1 + NXG*( j + NYG ) + nxt )
   bw_snd(j+NY) = g12( NXG*j + nxt )

   be_snd(j)    = f0( NX + NXG*( j + NYG )       ) &
                + f1( NX + NXG*( j + NYG ) + nxt ) &
                + f2( NX + NXG*( j + NYG ) + nxt ) &
                + f3( NX + NXG*( j + NYG ) + nxt ) &
                + f4( NX + NXG*( j + NYG ) + nxt ) &
                + f5( NX + NXG*( j + NYG ) + nxt ) &
                + f6( NX + NXG*( j + NYG ) + nxt )
   be_snd(j+NY) = g13( NX + 1 + NXG*j + nxt )
 END DO
!$OMP SINGLE
 CALL MPI_ISEND( bw_snd, yedge2, MPI_DOUBLE_PRECISION, bw, TAG11, &
                 MPI_COMM_VGRID, MPI_REQ_EY(5), MPI_ERR )
 CALL MPI_ISEND( be_snd, yedge2, MPI_DOUBLE_PRECISION, be, TAG12, &
                 MPI_COMM_VGRID, MPI_REQ_EY(6), MPI_ERR )
 CALL MPI_ISEND( tw_snd, yedge2, MPI_DOUBLE_PRECISION, tw, TAG13, &
                 MPI_COMM_VGRID, MPI_REQ_EY(7), MPI_ERR )
 CALL MPI_ISEND( te_snd, yedge2, MPI_DOUBLE_PRECISION, te, TAG14, &
                 MPI_COMM_VGRID, MPI_REQ_EY(8), MPI_ERR )
!$OMP END SINGLE

! Edges along Z
!$OMP DO
 DO k = 1, NZ
   ne_snd(k)    = f0( NX + NXG*( NY + NYG*k )       ) &
                + f1( NX + NXG*( NY + NYG*k ) + nxt ) &
                + f2( NX + NXG*( NY + NYG*k ) + nxt ) &
                + f3( NX + NXG*( NY + NYG*k ) + nxt ) &
                + f4( NX + NXG*( NY + NYG*k ) + nxt ) &
                + f5( NX + NXG*( NY + NYG*k ) + nxt ) &
                + f6( NX + NXG*( NY + NYG*k ) + nxt )
   ne_snd(k+NZ) = g7( NX+1 + NXG*( NY+1 + NYG*k ) + nxt )

   nw_snd(k)    = f0( 1 + NXG*( NY + NYG*k )       ) &
                + f1( 1 + NXG*( NY + NYG*k ) + nxt ) &
                + f2( 1 + NXG*( NY + NYG*k ) + nxt ) &
                + f3( 1 + NXG*( NY + NYG*k ) + nxt ) &
                + f4( 1 + NXG*( NY + NYG*k ) + nxt ) &
                + f5( 1 + NXG*( NY + NYG*k ) + nxt ) &
                + f6( 1 + NXG*( NY + NYG*k ) + nxt )
   nw_snd(k+NZ) = g10( NXG*( NY+1 + NYG*k ) + nxt )

   sw_snd(k)    = f0( 1 + NXG*( 1 +NYG*k )       ) &
                + f1( 1 + NXG*( 1 +NYG*k ) + nxt ) &
                + f2( 1 + NXG*( 1 +NYG*k ) + nxt ) &
                + f3( 1 + NXG*( 1 +NYG*k ) + nxt ) &
                + f4( 1 + NXG*( 1 +NYG*k ) + nxt ) &
                + f5( 1 + NXG*( 1 +NYG*k ) + nxt ) &
                + f6( 1 + NXG*( 1 +NYG*k ) + nxt )
   sw_snd(k+NZ) = g8( NXG*NYG*k + nxt )

   se_snd(k)    = f0( NX + NXG*( 1 + NYG*k )       ) &
                + f1( NX + NXG*( 1 + NYG*k ) + nxt ) &
                + f2( NX + NXG*( 1 + NYG*k ) + nxt ) &
                + f3( NX + NXG*( 1 + NYG*k ) + nxt ) &
                + f4( NX + NXG*( 1 + NYG*k ) + nxt ) &
                + f5( NX + NXG*( 1 + NYG*k ) + nxt ) &
                + f6( NX + NXG*( 1 + NYG*k ) + nxt )
   se_snd(k+NZ) = g9( NX+1 + NXG*NYG*k + nxt )
 END DO
!$OMP SINGLE
 CALL MPI_ISEND( sw_snd, zedge2, MPI_DOUBLE_PRECISION, sw, TAG15, &
                 MPI_COMM_VGRID, MPI_REQ_EZ(5), MPI_ERR )
 CALL MPI_ISEND( se_snd, zedge2, MPI_DOUBLE_PRECISION, se, TAG16, &
                 MPI_COMM_VGRID, MPI_REQ_EZ(6), MPI_ERR )
 CALL MPI_ISEND( nw_snd, zedge2, MPI_DOUBLE_PRECISION, nw, TAG17, &
                 MPI_COMM_VGRID, MPI_REQ_EZ(7), MPI_ERR )
 CALL MPI_ISEND( ne_snd, zedge2, MPI_DOUBLE_PRECISION, ne, TAG18, &
                 MPI_COMM_VGRID, MPI_REQ_EZ(8), MPI_ERR )

!-------- Update order parameter and inward values of g in real nodes ----------
! X direction
 CALL MPI_WAITALL(4, MPI_REQ_X, MPI_STAT, MPI_ERR)
!$OMP END SINGLE
!$OMP BARRIER
!$OMP DO PRIVATE(j)
 DO k = 1, NZ
   DO j = 1, NY
     phi( NX+1 + NXG*( j + NYG*k ) )     = east_rcv(j+NY*(k-1))
   END DO
   DO j = 1, NY
     g2( NX + NXG*( j + NYG*k ) + nxt )  = east_rcv(j+NY*(k-1+NZ))
   END DO
   DO j = 1, NY
     g8( NX + NXG*( j + NYG*k ) + nxt )  = east_rcv(j+NY*(k-1+2*NZ))
   END DO
   DO j = 1, NY
     g10( NX + NXG*( j + NYG*k ) + nxt ) = east_rcv(j+NY*(k-1+3*NZ))
   END DO
   DO j = 1, NY
     g12( NX + NXG*( j + NYG*k ) + nxt ) = east_rcv(j+NY*(k-1+4*NZ))
   END DO
   DO j = 1, NY
     g14( NX + NXG*( j + NYG*k ) + nxt ) = east_rcv(j+NY*(k-1+5*NZ))
   END DO

   DO j = 1, NY
     phi( NXG*( j + NYG*k ) ) = west_rcv(j+NY*(k-1))
   END DO
   DO j = 1, NY
     g1( 1 + NXG*( j + NYG*k ) + nxt )  = west_rcv(j+NY*(k-1+NZ))
   END DO
   DO j = 1, NY
     g7( 1 + NXG*( j + NYG*k ) + nxt )  = west_rcv(j+NY*(k-1+2*NZ))
   END DO
   DO j = 1, NY
     g9( 1 + NXG*( j + NYG*k ) + nxt )  = west_rcv(j+NY*(k-1+3*NZ))
   END DO
   DO j = 1, NY
     g11( 1 + NXG*( j + NYG*k ) + nxt ) = west_rcv(j+NY*(k-1+4*NZ))
   END DO
   DO j = 1, NY
     g13( 1 + NXG*( j + NYG*k ) + nxt ) = west_rcv(j+NY*(k-1+5*NZ))
   END DO
 END DO

! Y direction
!$OMP SINGLE
 CALL MPI_WAITALL(4, MPI_REQ_Y, MPI_STAT, MPI_ERR)
!$OMP END SINGLE
!$OMP BARRIER
!$OMP DO PRIVATE(i)
 DO k = 1, NZ
   DO i = 1, NX
     phi( i + NXG*( NY+1 + NYG*k ) )     = north_rcv(i+NX*(k-1))
   END DO
   DO i = 1, NX
     g4( i + NXG*( NY + NYG*k ) + nxt ) = north_rcv(i+NX*(k-1+NZ))
   END DO
   DO i = 1, NX
     g8( i + NXG*( NY + NYG*k ) + nxt ) = north_rcv(i+NX*(k-1+2*NZ))
   END DO
   DO i = 1, NX
     g9( i + NXG*( NY + NYG*k ) + nxt ) = north_rcv(i+NX*(k-1+3*NZ))
   END DO
   DO i = 1, NX
     g16( i + NXG*( NY + NYG*k ) + nxt ) = north_rcv(i+NX*(k-1+4*NZ))
   END DO
   DO i = 1, NX
     g18( i + NXG*( NY + NYG*k ) + nxt ) = north_rcv(i+NX*(k-1+5*NZ))
   END DO

   DO i = 1, NX
     phi( i + NXG*( NYG*k ) )     = south_rcv(i+NX*(k-1))
   END DO
   DO i = 1, NX
     g3( i + NXG*( 1 + NYG*k ) + nxt ) = south_rcv(i+NX*(k-1+NZ))
   END DO
   DO i = 1, NX
     g7( i + NXG*( 1 + NYG*k ) + nxt ) = south_rcv(i+NX*(k-1+2*NZ))
   END DO
   DO i = 1, NX
     g10( i + NXG*( 1 + NYG*k ) + nxt ) = south_rcv(i+NX*(k-1+3*NZ))
   END DO
   DO i = 1, NX
     g15( i + NXG*( 1 + NYG*k ) + nxt ) = south_rcv(i+NX*(k-1+4*NZ))
   END DO
   DO i = 1, NX
     g17( i + NXG*( 1 + NYG*k ) + nxt ) = south_rcv(i+NX*(k-1+5*NZ))
   END DO
 END DO

! Z direction
!$OMP SINGLE
 CALL MPI_WAITALL(4, MPI_REQ_Z, MPI_STAT, MPI_ERR)
!$OMP END SINGLE
!$OMP BARRIER
!$OMP DO PRIVATE(i)
 DO j = 1, NY
   DO i = 1, NX
     phi( i + NXG*( j + NYG*(NZ+1) ) )     = top_rcv(i+NX*(j-1))
   END DO
   DO i = 1, NX
     g6( i + NXG*( j + NYG*NZ ) + nxt ) = top_rcv(i+NX*(j-1+NY))
   END DO
   DO i = 1, NX
     g12( i + NXG*( j + NYG*NZ ) + nxt ) = top_rcv(i+NX*(j-1+2*NY))
   END DO
   DO i = 1, NX
     g13( i + NXG*( j + NYG*NZ ) + nxt ) = top_rcv(i+NX*(j-1+3*NY))
   END DO
   DO i = 1, NX
     g16( i + NXG*( j + NYG*NZ ) + nxt ) = top_rcv(i+NX*(j-1+4*NY))
   END DO
   DO i = 1, NX
     g17( i + NXG*( j + NYG*NZ ) + nxt ) = top_rcv(i+NX*(j-1+5*NY))
   END DO

   DO i = 1, NX
     phi( i + NXG*j )     = bot_rcv(i+NX*(j-1))
   END DO
   DO i = 1, NX
     g5( i + NXG*( j + NYG ) + nxt ) = bot_rcv(i+NX*(j-1+NY))
   END DO
   DO i = 1, NX
     g11( i + NXG*( j + NYG ) + nxt ) = bot_rcv(i+NX*(j-1+2*NY))
   END DO
   DO i = 1, NX
     g14( i + NXG*( j + NYG ) + nxt ) = bot_rcv(i+NX*(j-1+3*NY))
   END DO
   DO i = 1, NX
     g15( i + NXG*( j + NYG ) + nxt ) = bot_rcv(i+NX*(j-1+4*NY))
   END DO
   DO i = 1, NX
     g18( i + NXG*( j + NYG ) + nxt ) = bot_rcv(i+NX*(j-1+5*NY))
   END DO
 END DO

! Edges along x
!$OMP SINGLE
 CALL MPI_WAITALL(8, MPI_REQ_EX, MPI_ESTAT, MPI_ERR)
!$OMP END SINGLE
!$OMP BARRIER
!$OMP DO
 DO i = 1, NX
   phi( i + NXG*( NY+1 + NYG*(NZ+1) ) ) = tn_rcv(i)
   g16( i + NXG*( NY + NYG*NZ ) + nxt ) = tn_rcv(i+NX)

   phi( i + NXG*NYG*(NZ+1) )           = ts_rcv(i)
   g17( i + NXG*( 1 + NYG*NZ ) + nxt ) = ts_rcv(i+NX)

   phi( i )                         = bs_rcv(i)
   g15( i + NXG*( 1 + NYG ) + nxt ) = bs_rcv(i+NX)

   phi( i + NXG*(NY+1) )             = bn_rcv(i)
   g18( i + NXG*( NY + NYG ) + nxt ) = bn_rcv(i+NX)
 END DO

! Edges along y
!$OMP SINGLE
 CALL MPI_WAITALL(8, MPI_REQ_EY, MPI_ESTAT, MPI_ERR)
!$OMP END SINGLE
!$OMP BARRIER
!$OMP DO
 DO j = 1, NY
   phi( NX + 1 + NXG*(j + NYG*(NZ+1) ) ) = te_rcv(j)
   g12( NX + NXG*( j + NYG*NZ ) + nxt )  = te_rcv(j+NY)

   phi( NXG*( j + NYG*(NZ+1) ) )       = tw_rcv(j)
   g13( 1 + NXG*( j + NYG*NZ ) + nxt ) = tw_rcv(j+NY)

   phi( NXG*j )                     = bw_rcv(j)
   g11( 1 + NXG*( j + NYG ) + nxt ) = bw_rcv(j+NY)

   phi( NX+1 + NXG*j )               = be_rcv(j)
   g14( NX + NXG*( j + NYG ) + nxt ) = be_rcv(j+NY)
 END DO

! Edges along z
!$OMP SINGLE
 CALL MPI_WAITALL(8, MPI_REQ_EZ, MPI_ESTAT, MPI_ERR)
!$OMP END SINGLE
!$OMP BARRIER
!$OMP DO
 DO k = 1, NZ
   phi( NX+1 + NXG*( NY+1 + NYG*k ) )  = ne_rcv(k)
   g8( NX + NXG*( NY + NYG*k ) + nxt ) = ne_rcv(k+NZ)

   phi( NXG*( NY+1 + NYG*k ) )        = nw_rcv(k)
   g9( 1 + NXG*( NY + NYG*k ) + nxt ) = nw_rcv(k+NZ)

   phi( NXG*NYG*k )                  = sw_rcv(k)
   g7( 1 + NXG*( 1 + NYG*k ) + nxt ) = sw_rcv(k+NZ)

   phi( NX+1 + NXG*NYG*k )             = se_rcv(k)
   g10( NX + NXG*( 1 + NYG*k ) + nxt ) = se_rcv(k+NZ)
 END DO

!$OMP END PARALLEL

 RETURN
 END SUBROUTINE PostStream

