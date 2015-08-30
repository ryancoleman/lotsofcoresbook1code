!-------------------------------------------------------------------------------
! Subroutine : PostCollision
! Revision   : 1.5 (2014/09/02)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Update g in nodes at processor boundaries and outward f in ghosts.
! This is faster than using a neighbor array in collision.
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

 SUBROUTINE PostCollision

!  Common Variables
 USE Domain
 USE LBMParams
 USE MPIParams
 USE MPI
 IMPLICIT NONE

!  Local variables
 INTEGER :: i, j, k
 INTEGER :: MPI_ERR
 INTEGER :: MPI_REQ_X(4), MPI_REQ_Y(4), MPI_REQ_Z(4)
 INTEGER :: MPI_STAT(MPI_STATUS_SIZE,4)

!$OMP PARALLEL

!-------- Exchange data for outward f components ------------------
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
     f_west_snd(j+NY*(k-1)) = f1( 1 + NXG*( j + NYG*k ) + now )
   END DO
   DO j = 1, NY
     f_east_snd(j+NY*(k-1)) = f2( NX + NXG*( j + NYG*k ) + now )
   END DO
 END DO
!$OMP SINGLE
 CALL MPI_ISEND( f_west_snd, xsize, MPI_DOUBLE_PRECISION, west, TAG1, &
                 MPI_COMM_VGRID, MPI_REQ_X(3), MPI_ERR )
 CALL MPI_ISEND( f_east_snd, xsize, MPI_DOUBLE_PRECISION, east, TAG2, &
                 MPI_COMM_VGRID, MPI_REQ_X(4), MPI_ERR )
!$OMP END SINGLE

! Y direction
!$OMP DO PRIVATE(i)
 DO k = 1, NZ
   DO i = 1, NX
     f_south_snd(i+NX*(k-1)) = f3( i + NXG*( 1 + NYG*k ) + now )
   END DO
   DO i = 1, NX
     f_north_snd(i+NX*(k-1)) = f4(i + NXG*( NY + NYG*k) + now)
   END DO
 END DO
!$OMP SINGLE
 CALL MPI_ISEND( f_south_snd, ysize , MPI_DOUBLE_PRECISION, south, TAG3, &
                 MPI_COMM_VGRID, MPI_REQ_Y(3), MPI_ERR )
 CALL MPI_ISEND( f_north_snd, ysize , MPI_DOUBLE_PRECISION, north, TAG4, &
                 MPI_COMM_VGRID, MPI_REQ_Y(4), MPI_ERR )
!$OMP END SINGLE

! Z direction
!$OMP DO PRIVATE(i)
 DO j = 1, NY
   DO i = 1, NX
     f_bot_snd(i+NX*(j-1)) = f5( i + NXG*( j + NYG ) + now )
   END DO
   DO i = 1, NX
     f_top_snd(i+NX*(j-1)) = f6( i + NXG*( j + NYG*NZ ) + now )
   END DO
 END DO
!$OMP SINGLE
 CALL MPI_ISEND( f_bot_snd, zsize , MPI_DOUBLE_PRECISION, bot, TAG5, &
                 MPI_COMM_VGRID, MPI_REQ_Z(3), MPI_ERR )
 CALL MPI_ISEND( f_top_snd, zsize , MPI_DOUBLE_PRECISION, top, TAG6, &
                 MPI_COMM_VGRID, MPI_REQ_Z(4), MPI_ERR )

!-------- Update outward f components for the ghosts -------------------------
! X direction
 CALL MPI_WAITALL(4, MPI_REQ_X, MPI_STAT, MPI_ERR)
!$OMP END SINGLE
!$OMP BARRIER
!$OMP DO PRIVATE(j)
 DO k = 1, NZ
   DO j = 1, NY
     f1( (NX+1) + NXG*( j + NYG*k ) + now ) = f_east_rcv(j+NY*(k-1))
   END DO
   DO j = 1, NY
     f2( NXG*( j + NYG*k) + now ) = f_west_rcv(j+NY*(k-1))
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
     f3( i + NXG*( (NY+1) + NYG*k ) + now ) = f_north_rcv(i+NX*(k-1))
   END DO
   DO i = 1, NX
     f4( i + NXG*( NYG*k ) + now ) = f_south_rcv(i+NX*(k-1))
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
     f5( i + NXG*( j + NYG*(NZ+1) ) + now ) = f_top_rcv(i+NX*(j-1))
   END DO
   DO i = 1, NX
     f6( i + NXG*j + now ) = f_bot_rcv(i+NX*(j-1))
   END DO
 END DO

!
! At this stage we are ready to do the relaxation step in stream.f90
!

!$OMP END PARALLEL

 RETURN
 END SUBROUTINE PostCollision

