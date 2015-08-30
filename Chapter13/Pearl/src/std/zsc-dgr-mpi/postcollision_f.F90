!-------------------------------------------------------------------------------
! Subroutine : PostCollisionF
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Update outward f in ghost layers using MPI.
!> @details
!! Copy inward f values from the real nodes into outward f values in the ghosts.
!! These f values are required for the relaxation step in "stream.f".
!!
!! MPI exchange is done using non-blocking receive and blocking send pairs
!! (completed by a wait statement) in alternate directions to avoid deadlock.
!! The components of f that need to be exchanged are packed into buffer arrays
!! before the exchange.

!-------------------------------------------------------------------------------
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

 SUBROUTINE PostCollisionF

! Common Variables
 USE Domain,    ONLY : xl_f, xlg_f, xsize_f, xu_f, xug_f, yl_f, ylg_f, ysize_f, yu_f, yug_f
 USE LBMParams, ONLY : fcol
 USE MPIParams
 USE MPI
 IMPLICIT NONE

! Local Variables
 INTEGER :: count, i, j
 INTEGER :: MPI_ERR, MPI_REQ
 INTEGER :: status(MPI_STATUS_SIZE)


!-------- Save data for inward f components ------------------------------------
! X DIRECTION
 count = 1
 DO j = yl_f, yu_f
   f_west_snd(count) = fcol(xl_f,j,1)
   f_east_snd(count) = fcol(xu_f,j,3)
   count = count + 1
 END DO

! Y DIRECTION
 count = 1
 DO i = xl_f, xu_f
   f_south_snd(count) = fcol(i,yl_f,2)
   f_north_snd(count) = fcol(i,yu_f,4)
   count = count + 1
 END DO


!-------- Exchange information using MPI calls ---------------------------------
! X DIRECTION
 CALL MPI_IRECV(f_east_rcv, xsize_f, MPI_DOUBLE_PRECISION, east, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_west_snd, xsize_f, MPI_DOUBLE_PRECISION, west, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_west_rcv, xsize_f, MPI_DOUBLE_PRECISION, west, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_east_snd, xsize_f, MPI_DOUBLE_PRECISION, east, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! Y DIRECTION
 CALL MPI_IRECV(f_north_rcv, ysize_f, MPI_DOUBLE_PRECISION, north, TAG3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_south_snd, ysize_f, MPI_DOUBLE_PRECISION, south, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_south_rcv, ysize_f, MPI_DOUBLE_PRECISION, south, TAG4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_north_snd, ysize_f, MPI_DOUBLE_PRECISION, north, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)


!-------- Update outward pointing f components for the ghosts ------------------
! x direction
 count = 1
 DO j = yl_f, yu_f
   fcol(xug_f,j,1) = f_east_rcv(count)
   fcol(xlg_f,j,3) = f_west_rcv(count)
   count = count + 1
 END DO

! y direction
 count = 1
 DO i = xl_f, xu_f
   fcol(i,yug_f,2) = f_north_rcv(count)
   fcol(i,ylg_f,4) = f_south_rcv(count)
   count = count + 1
 END DO

 RETURN
 END SUBROUTINE PostCollisionF

