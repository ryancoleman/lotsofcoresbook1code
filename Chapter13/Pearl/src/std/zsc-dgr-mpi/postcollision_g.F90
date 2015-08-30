!-------------------------------------------------------------------------------
! Subroutine : PostCollisionG
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Update g in nodes at processor boundaries and (u,p) at ghosts using MPI.
!> @details
!! Copy outward g values from the ghosts into inward g values in the real nodes.
!! Copy values of pressure and velocity from the real nodes into the ghost
!! layers. These (u,p) values are required for the interpolation step in
!! "update.f".
!!
!! MPI exchange is done using non-blocking receive and blocking send pairs
!! (completed by a wait statement) in alternate directions to avoid deadlock.
!! The components of g, u, and p that need to be exchanged are packed into
!!  buffer arrays before the exchange.

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

 SUBROUTINE PostCollisionG

! Common Variables
 USE Domain,      ONLY : nxt, xl, xlg, xsize6, xu, xug, yl, ylg, ysize6, yu, yug
 USE FluidParams, ONLY : p, u
 USE LBMParams,   ONLY : g
 USE MPIParams
 USE MPI
 IMPLICIT NONE

! Local Variables
 INTEGER :: count, i, j
 INTEGER :: MPI_ERR, MPI_REQ
 INTEGER :: status(MPI_STATUS_SIZE)


!-------- Save data for outward g components from the ghost nodes and ----------
!-------- velocity and pressure from the real nodes at the boundaries ----------
! X DIRECTION
 count = 1
 DO j = yl, yu
   hydro_west_snd(count)     = g(xlg,j,3,nxt)
   hydro_west_snd(count + 1) = g(xlg,j,6,nxt)
   hydro_west_snd(count + 2) = g(xlg,j,7,nxt)
   hydro_west_snd(count + 3) = u(xl,j,1)
   hydro_west_snd(count + 4) = u(xl,j,2)
   hydro_west_snd(count + 5) = p(xl,j)

   hydro_east_snd(count)     = g(xug,j,1,nxt)
   hydro_east_snd(count + 1) = g(xug,j,5,nxt)
   hydro_east_snd(count + 2) = g(xug,j,8,nxt)
   hydro_east_snd(count + 3) = u(xu,j,1)
   hydro_east_snd(count + 4) = u(xu,j,2)
   hydro_east_snd(count + 5) = p(xu,j)
   count = count + 6
 END DO

! Y DIRECTION
 count = 1
 DO i = xl, xu
   hydro_south_snd(count)     = g(i,ylg,4,nxt)
   hydro_south_snd(count + 1) = g(i,ylg,7,nxt)
   hydro_south_snd(count + 2) = g(i,ylg,8,nxt)
   hydro_south_snd(count + 3) = u(i,yl,1)
   hydro_south_snd(count + 4) = u(i,yl,2)
   hydro_south_snd(count + 5) = p(i,yl)

   hydro_north_snd(count)     = g(i,yug,2,nxt)
   hydro_north_snd(count + 1) = g(i,yug,5,nxt)
   hydro_north_snd(count + 2) = g(i,yug,6,nxt)
   hydro_north_snd(count + 3) = u(i,yu,1)
   hydro_north_snd(count + 4) = u(i,yu,2)
   hydro_north_snd(count + 5) = p(i,yu)
   count = count + 6
 END DO

! NORTH-EAST DIAGONAL
 hydro_ne_snd(1) = g(xug,yug,5,nxt)
 hydro_ne_snd(2) = u(xu,yu,1)
 hydro_ne_snd(3) = u(xu,yu,2)
 hydro_ne_snd(4) = p(xu,yu)

! NORTH-WEST DIAGONAL
 hydro_nw_snd(1) = g(xlg,yug,6,nxt)
 hydro_nw_snd(2) = u(xl,yu,1)
 hydro_nw_snd(3) = u(xl,yu,2)
 hydro_nw_snd(4) = p(xl,yu)

! SOUTH-WEST DIAGONAL
 hydro_sw_snd(1) = g(xlg,ylg,7,nxt)
 hydro_sw_snd(2) = u(xl,yl,1)
 hydro_sw_snd(3) = u(xl,yl,2)
 hydro_sw_snd(4) = p(xl,yl)

! SOUTH-EAST DIAGONAL
 hydro_se_snd(1) = g(xug,ylg,8,nxt)
 hydro_se_snd(2) = u(xu,yl,1)
 hydro_se_snd(3) = u(xu,yl,2)
 hydro_se_snd(4) = p(xu,yl)

!-------- Exchange information using MPI calls ---------------------------------
! X DIRECTION
 CALL MPI_IRECV(hydro_east_rcv, xsize6, MPI_DOUBLE_PRECISION, east, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( hydro_west_snd, xsize6, MPI_DOUBLE_PRECISION, west, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(hydro_west_rcv, xsize6, MPI_DOUBLE_PRECISION, west, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( hydro_east_snd, xsize6, MPI_DOUBLE_PRECISION, east, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! Y DIRECTION
 CALL MPI_IRECV(hydro_north_rcv, ysize6, MPI_DOUBLE_PRECISION, north, TAG3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( hydro_south_snd, ysize6, MPI_DOUBLE_PRECISION, south, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(hydro_south_rcv, ysize6, MPI_DOUBLE_PRECISION, south, TAG4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( hydro_north_snd, ysize6, MPI_DOUBLE_PRECISION, north, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! DIAGONALS
 CALL MPI_IRECV(hydro_ne_rcv, 4, MPI_DOUBLE_PRECISION, ne, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( hydro_sw_snd, 4, MPI_DOUBLE_PRECISION, sw, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(hydro_nw_rcv, 4, MPI_DOUBLE_PRECISION, nw, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( hydro_se_snd, 4, MPI_DOUBLE_PRECISION, se, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(hydro_se_rcv, 4, MPI_DOUBLE_PRECISION, se, TAG3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( hydro_nw_snd, 4, MPI_DOUBLE_PRECISION, nw, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(hydro_sw_rcv, 4, MPI_DOUBLE_PRECISION, sw, TAG4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( hydro_ne_snd, 4, MPI_DOUBLE_PRECISION, ne, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

!-------- Update inward components of g in the real nodes at the boundaries ----
!-------- and velocity and pressure at the ghost nodes                      ----
! X DIRECTION
 count = 1
 DO j = yl, yu
   g(xu,j,3,nxt) = hydro_east_rcv(count)
   g(xu,j,6,nxt) = hydro_east_rcv(count + 1)
   g(xu,j,7,nxt) = hydro_east_rcv(count + 2)
   u(xug,j,1)    = hydro_east_rcv(count + 3)
   u(xug,j,2)    = hydro_east_rcv(count + 4)
   p(xug,j)      = hydro_east_rcv(count + 5)

   g(xl,j,1,nxt) = hydro_west_rcv(count)
   g(xl,j,5,nxt) = hydro_west_rcv(count + 1)
   g(xl,j,8,nxt) = hydro_west_rcv(count + 2)
   u(xlg,j,1)    = hydro_west_rcv(count + 3)
   u(xlg,j,2)    = hydro_west_rcv(count + 4)
   p(xlg,j)      = hydro_west_rcv(count + 5)
   count = count + 6
 END DO

! Y DIRECTION
 count = 1
 DO i = xl, xu
   g(i,yu,4,nxt) = hydro_north_rcv(count)
   g(i,yu,7,nxt) = hydro_north_rcv(count + 1)
   g(i,yu,8,nxt) = hydro_north_rcv(count + 2)
   u(i,yug,1)    = hydro_north_rcv(count + 3)
   u(i,yug,2)    = hydro_north_rcv(count + 4)
   p(i,yug)      = hydro_north_rcv(count + 5)

   g(i,yl,2,nxt) = hydro_south_rcv(count)
   g(i,yl,5,nxt) = hydro_south_rcv(count + 1)
   g(i,yl,6,nxt) = hydro_south_rcv(count + 2)
   u(i,ylg,1)    = hydro_south_rcv(count + 3)
   u(i,ylg,2)    = hydro_south_rcv(count + 4)
   p(i,ylg)      = hydro_south_rcv(count + 5)
   count = count + 6
 END DO

! NORTH-EAST DIAGONAL
 g(xu,yu,7,nxt) = hydro_ne_rcv(1)
 u(xug,yug,1)   = hydro_ne_rcv(2)
 u(xug,yug,2)   = hydro_ne_rcv(3)
 p(xug,yug)     = hydro_ne_rcv(4)

! NORTH-WEST DIAGONAL
 g(xl,yu,8,nxt) = hydro_nw_rcv(1)
 u(xlg,yug,1)   = hydro_nw_rcv(2)
 u(xlg,yug,2)   = hydro_nw_rcv(3)
 p(xlg,yug)     = hydro_nw_rcv(4)

! SOUTH-WEST DIAGONAL
 g(xl,yl,5,nxt) = hydro_sw_rcv(1)
 u(xlg,ylg,1)   = hydro_sw_rcv(2)
 u(xlg,ylg,2)   = hydro_sw_rcv(3)
 p(xlg,ylg)     = hydro_sw_rcv(4)

! SOUTH-EAST DIAGONAL
 g(xu,yl,6,nxt) = hydro_se_rcv(1)
 u(xug,ylg,1)   = hydro_se_rcv(2)
 u(xug,ylg,2)   = hydro_se_rcv(3)
 p(xug,ylg)     = hydro_se_rcv(4)


 RETURN
 END SUBROUTINE PostCollisionG

