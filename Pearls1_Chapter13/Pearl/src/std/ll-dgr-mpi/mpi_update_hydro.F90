!-------------------------------------------------------------------------------
! Subroutine : MPI_UpdateHydro
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Update values of the velocity and the pressure in the ghost layers.
!> @details
!! Copy values of the velocity and the pressure from the real nodes at the
!! processor boundaries to the ghost layers. Only the arrays corresponding to
!! the coarse grid (pressure distribution function) are updated.

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

 SUBROUTINE MPI_UpdateHydro

! Common Variables
 USE Domain
 USE FluidParams
 USE LBMParams
 USE MPIParams
 USE MPI
 IMPLICIT NONE

! Local Variables
 INTEGER :: count, i, j
 INTEGER :: MPI_ERR, MPI_REQ
 INTEGER :: status(MPI_STATUS_SIZE)


!-------- Save velocity and pressure from the real nodes at the boundaries -----
! X DIRECTION
 count = 1
 DO j = yl, yu
   hydro_west_snd(count    ) = u(xl,j,1)
   hydro_west_snd(count + 1) = u(xl,j,2)
   hydro_west_snd(count + 2) = p(xl,j)

   hydro_east_snd(count    ) = u(xu,j,1)
   hydro_east_snd(count + 1) = u(xu,j,2)
   hydro_east_snd(count + 2) = p(xu,j)
   count = count + 3
 END DO

! Y DIRECTION
 count = 1
 DO i = xl, xu
   hydro_south_snd(count    ) = u(i,yl,1)
   hydro_south_snd(count + 1) = u(i,yl,2)
   hydro_south_snd(count + 2) = p(i,yl)

   hydro_north_snd(count    ) = u(i,yu,1)
   hydro_north_snd(count + 1) = u(i,yu,2)
   hydro_north_snd(count + 2) = p(i,yu)
   count = count + 3
 END DO

! NORTH-EAST DIAGONAL
 hydro_ne_snd(1) = u(xu,yu,1)
 hydro_ne_snd(2) = u(xu,yu,2)
 hydro_ne_snd(3) = p(xu,yu)

! NORTH-WEST DIAGONAL
 hydro_nw_snd(1) = u(xl,yu,1)
 hydro_nw_snd(2) = u(xl,yu,2)
 hydro_nw_snd(3) = p(xl,yu)

! SOUTH-WEST DIAGONAL
 hydro_sw_snd(1) = u(xl,yl,1)
 hydro_sw_snd(2) = u(xl,yl,2)
 hydro_sw_snd(3) = p(xl,yl)

! SOUTH-EAST DIAGONAL
 hydro_se_snd(1) = u(xu,yl,1)
 hydro_se_snd(2) = u(xu,yl,2)
 hydro_se_snd(3) = p(xu,yl)

!-------- Exchange information using MPI calls ---------------------------------
! X DIRECTION
 CALL MPI_IRECV(hydro_east_rcv, xsize3, MPI_DOUBLE_PRECISION, east, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( hydro_west_snd, xsize3, MPI_DOUBLE_PRECISION, west, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(hydro_west_rcv, xsize3, MPI_DOUBLE_PRECISION, west, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( hydro_east_snd, xsize3, MPI_DOUBLE_PRECISION, east, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! Y DIRECTION
 CALL MPI_IRECV(hydro_north_rcv, ysize3, MPI_DOUBLE_PRECISION, north, TAG3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( hydro_south_snd, ysize3, MPI_DOUBLE_PRECISION, south, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(hydro_south_rcv, ysize3, MPI_DOUBLE_PRECISION, south, TAG4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( hydro_north_snd, ysize3, MPI_DOUBLE_PRECISION, north, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! DIAGONALS
 CALL MPI_IRECV(hydro_ne_rcv, 3, MPI_DOUBLE_PRECISION, ne, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( hydro_sw_snd, 3, MPI_DOUBLE_PRECISION, sw, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(hydro_nw_rcv, 3, MPI_DOUBLE_PRECISION, nw, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( hydro_se_snd, 3, MPI_DOUBLE_PRECISION, se, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(hydro_se_rcv, 3, MPI_DOUBLE_PRECISION, se, TAG3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( hydro_nw_snd, 3, MPI_DOUBLE_PRECISION, nw, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(hydro_sw_rcv, 3, MPI_DOUBLE_PRECISION, sw, TAG4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( hydro_ne_snd, 3, MPI_DOUBLE_PRECISION, ne, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)


!-------- Update velocity and pressure at the ghost nodes ----------------------
! X DIRECTION
 count = 1
 DO j = yl, yu
   u(xug,j,1) = hydro_east_rcv(count)
   u(xug,j,2) = hydro_east_rcv(count + 1)
   p(xug,j)   = hydro_east_rcv(count + 2)

   u(xlg,j,1) = hydro_west_rcv(count)
   u(xlg,j,2) = hydro_west_rcv(count + 1)
   p(xlg,j)   = hydro_west_rcv(count + 2)
   count = count + 3
 END DO

! Y DIRECTION
 count = 1
 DO i = xl, xu
   u(i,yug,1) = hydro_north_rcv(count)
   u(i,yug,2) = hydro_north_rcv(count + 1)
   p(i,yug)   = hydro_north_rcv(count + 2)

   u(i,ylg,1) = hydro_south_rcv(count)
   u(i,ylg,2) = hydro_south_rcv(count + 1)
   p(i,ylg)   = hydro_south_rcv(count + 2)
   count = count + 3
 END DO

! NORTH-EAST DIAGONAL
 u(xug,yug,1) = hydro_ne_rcv(1)
 u(xug,yug,2) = hydro_ne_rcv(2)
 p(xug,yug)   = hydro_ne_rcv(3)

! NORTH-WEST DIAGONAL
 u(xlg,yug,1) = hydro_nw_rcv(1)
 u(xlg,yug,2) = hydro_nw_rcv(2)
 p(xlg,yug)   = hydro_nw_rcv(3)

! SOUTH-WEST DIAGONAL
 u(xlg,ylg,1) = hydro_sw_rcv(1)
 u(xlg,ylg,2) = hydro_sw_rcv(2)
 p(xlg,ylg)   = hydro_sw_rcv(3)

! SOUTH-EAST DIAGONAL
 u(xug,ylg,1) = hydro_se_rcv(1)
 u(xug,ylg,2) = hydro_se_rcv(2)
 p(xug,ylg)   = hydro_se_rcv(3)


 RETURN
 END SUBROUTINE MPI_UpdateHydro

