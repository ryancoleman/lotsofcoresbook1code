!-------------------------------------------------------------------------------
! Subroutine : PostCollision
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Update g in nodes at processor boundaries and outward f in ghosts using MPI.
!> @details
!! Copy outward g values from the ghosts into inward g values in the real nodes.
!! Copy inward f values from the real nodes into outward f values in the ghosts.
!! These f values are required for the relaxation step in "stream.f".
!!
!! MPI exchange is done using non-blocking receive and blocking send pairs
!! (completed by a wait statement) in alternate directions to avoid deadlock.
!! The components of g and f that need to be exchanged are packed into buffer
!! arrays before the exchange.

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

 SUBROUTINE PostCollision

!  Common Variables
 USE Domain
 USE LBMParams, ONLY : f, g
 USE MPIParams
 USE MPI
 IMPLICIT NONE

!  Local variables
 INTEGER :: count, i, j, k
 INTEGER :: MPI_ERR
 INTEGER :: MPI_REQ(4)
 INTEGER :: status(MPI_STATUS_SIZE)

!-------- Exchange data for outward f and inward g components ------------------

! X direction
 CALL MPI_IRECV(f_east_rcv, xsize , MPI_DOUBLE_PRECISION, east, TAG1, MPI_COMM_VGRID, MPI_REQ(1), MPI_ERR)
 CALL MPI_IRECV(g_east_rcv, xsize5, MPI_DOUBLE_PRECISION, east, TAG2, MPI_COMM_VGRID, MPI_REQ(2), MPI_ERR)
 CALL MPI_IRECV(f_west_rcv, xsize , MPI_DOUBLE_PRECISION, west, TAG3, MPI_COMM_VGRID, MPI_REQ(3), MPI_ERR)
 CALL MPI_IRECV(g_west_rcv, xsize5, MPI_DOUBLE_PRECISION, west, TAG4, MPI_COMM_VGRID, MPI_REQ(4), MPI_ERR)
 count = 1
 DO k = zl, zu
   DO j = yl, yu
     f_west_snd(count)          = f(xl ,j,k, 1,now)
     g_west_snd(count)          = g(xlg,j,k, 2,nxt)
     g_west_snd(count + xsize)  = g(xlg,j,k, 8,nxt)
     g_west_snd(count + xsize2) = g(xlg,j,k,10,nxt)
     g_west_snd(count + xsize3) = g(xlg,j,k,12,nxt)
     g_west_snd(count + xsize4) = g(xlg,j,k,14,nxt)

     f_east_snd(count)          = f(xu ,j,k, 2,now)
     g_east_snd(count)          = g(xug,j,k, 1,nxt)
     g_east_snd(count + xsize)  = g(xug,j,k, 7,nxt)
     g_east_snd(count + xsize2) = g(xug,j,k, 9,nxt)
     g_east_snd(count + xsize3) = g(xug,j,k,11,nxt)
     g_east_snd(count + xsize4) = g(xug,j,k,13,nxt)
     count = count + 1
   END DO
 END DO
 CALL MPI_SEND(f_west_snd, xsize , MPI_DOUBLE_PRECISION, west, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(g_west_snd, xsize5, MPI_DOUBLE_PRECISION, west, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(f_east_snd, xsize , MPI_DOUBLE_PRECISION, east, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(g_east_snd, xsize5, MPI_DOUBLE_PRECISION, east, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(1), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(2), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(3), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(4), status, MPI_ERR)

! Y direction
 CALL MPI_IRECV(f_north_rcv, ysize , MPI_DOUBLE_PRECISION, north, TAG1, MPI_COMM_VGRID, MPI_REQ(1), MPI_ERR)
 CALL MPI_IRECV(g_north_rcv, ysize5, MPI_DOUBLE_PRECISION, north, TAG2, MPI_COMM_VGRID, MPI_REQ(2), MPI_ERR)
 CALL MPI_IRECV(f_south_rcv, ysize , MPI_DOUBLE_PRECISION, south, TAG3, MPI_COMM_VGRID, MPI_REQ(3), MPI_ERR)
 CALL MPI_IRECV(g_south_rcv, ysize5, MPI_DOUBLE_PRECISION, south, TAG4, MPI_COMM_VGRID, MPI_REQ(4), MPI_ERR)
 count = 1
 DO k = zl, zu
   DO i = xl, xu
     f_south_snd(count)          = f(i,yl ,k, 3,now)
     g_south_snd(count)          = g(i,ylg,k, 4,nxt)
     g_south_snd(count + ysize)  = g(i,ylg,k, 8,nxt)
     g_south_snd(count + ysize2) = g(i,ylg,k, 9,nxt)
     g_south_snd(count + ysize3) = g(i,ylg,k,16,nxt)
     g_south_snd(count + ysize4) = g(i,ylg,k,18,nxt)

     f_north_snd(count)          = f(i,yu ,k, 4,now)
     g_north_snd(count)          = g(i,yug,k, 3,nxt)
     g_north_snd(count + ysize)  = g(i,yug,k, 7,nxt)
     g_north_snd(count + ysize2) = g(i,yug,k,10,nxt)
     g_north_snd(count + ysize3) = g(i,yug,k,15,nxt)
     g_north_snd(count + ysize4) = g(i,yug,k,17,nxt)
     count = count + 1
   END DO
 END DO
 CALL MPI_SEND(f_south_snd, ysize , MPI_DOUBLE_PRECISION, south, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(g_south_snd, ysize5, MPI_DOUBLE_PRECISION, south, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(f_north_snd, ysize , MPI_DOUBLE_PRECISION, north, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(g_north_snd, ysize5, MPI_DOUBLE_PRECISION, north, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(1), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(2), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(3), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(4), status, MPI_ERR)

! Z direction
 CALL MPI_IRECV(f_top_rcv, zsize , MPI_DOUBLE_PRECISION, top, TAG1, MPI_COMM_VGRID, MPI_REQ(1), MPI_ERR)
 CALL MPI_IRECV(g_top_rcv, zsize5, MPI_DOUBLE_PRECISION, top, TAG2, MPI_COMM_VGRID, MPI_REQ(2), MPI_ERR)
 CALL MPI_IRECV(f_bot_rcv, zsize , MPI_DOUBLE_PRECISION, bot, TAG3, MPI_COMM_VGRID, MPI_REQ(3), MPI_ERR)
 CALL MPI_IRECV(g_bot_rcv, zsize5, MPI_DOUBLE_PRECISION, bot, TAG4, MPI_COMM_VGRID, MPI_REQ(4), MPI_ERR)
 count = 1
 DO j = yl, yu
   DO i = xl, xu
     f_bot_snd(count)          = f(i,j,zl , 5,now)
     g_bot_snd(count)          = g(i,j,zlg, 6,nxt)
     g_bot_snd(count + zsize)  = g(i,j,zlg,12,nxt)
     g_bot_snd(count + zsize2) = g(i,j,zlg,13,nxt)
     g_bot_snd(count + zsize3) = g(i,j,zlg,16,nxt)
     g_bot_snd(count + zsize4) = g(i,j,zlg,17,nxt)

     f_top_snd(count)          = f(i,j,zu , 6,now)
     g_top_snd(count)          = g(i,j,zug, 5,nxt)
     g_top_snd(count + zsize)  = g(i,j,zug,11,nxt)
     g_top_snd(count + zsize2) = g(i,j,zug,14,nxt)
     g_top_snd(count + zsize3) = g(i,j,zug,15,nxt)
     g_top_snd(count + zsize4) = g(i,j,zug,18,nxt)
     count = count + 1
   END DO
 END DO
 CALL MPI_SEND(f_bot_snd, zsize , MPI_DOUBLE_PRECISION, bot, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(g_bot_snd, zsize5, MPI_DOUBLE_PRECISION, bot, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(f_top_snd, zsize , MPI_DOUBLE_PRECISION, top, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(g_top_snd, zsize5, MPI_DOUBLE_PRECISION, top, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(1), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(2), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(3), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(4), status, MPI_ERR)

! Edges along X
 CALL MPI_IRECV(g_tn_rcv, xedge, MPI_DOUBLE_PRECISION, tn, TAG1, MPI_COMM_VGRID, MPI_REQ(1), MPI_ERR)
 CALL MPI_IRECV(g_ts_rcv, xedge, MPI_DOUBLE_PRECISION, ts, TAG2, MPI_COMM_VGRID, MPI_REQ(2), MPI_ERR)
 CALL MPI_IRECV(g_bn_rcv, xedge, MPI_DOUBLE_PRECISION, bn, TAG3, MPI_COMM_VGRID, MPI_REQ(3), MPI_ERR)
 CALL MPI_IRECV(g_bs_rcv, xedge, MPI_DOUBLE_PRECISION, bs, TAG4, MPI_COMM_VGRID, MPI_REQ(4), MPI_ERR)
 count = 1
 DO i = xl, xu
   g_tn_snd(count) = g(i,yug,zug,15,nxt)
   g_ts_snd(count) = g(i,ylg,zug,18,nxt)
   g_bs_snd(count) = g(i,ylg,zlg,16,nxt)
   g_bn_snd(count) = g(i,yug,zlg,17,nxt)
   count = count + 1
 END DO
 CALL MPI_SEND(g_bs_snd, xedge, MPI_DOUBLE_PRECISION, bs, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(g_bn_snd, xedge, MPI_DOUBLE_PRECISION, bn, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(g_ts_snd, xedge, MPI_DOUBLE_PRECISION, ts, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(g_tn_snd, xedge, MPI_DOUBLE_PRECISION, tn, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(1), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(2), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(3), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(4), status, MPI_ERR)

! Edges along Y
 CALL MPI_IRECV(g_te_rcv, yedge, MPI_DOUBLE_PRECISION, te, TAG1, MPI_COMM_VGRID, MPI_REQ(1), MPI_ERR)
 CALL MPI_IRECV(g_tw_rcv, yedge, MPI_DOUBLE_PRECISION, tw, TAG2, MPI_COMM_VGRID, MPI_REQ(2), MPI_ERR)
 CALL MPI_IRECV(g_be_rcv, yedge, MPI_DOUBLE_PRECISION, be, TAG3, MPI_COMM_VGRID, MPI_REQ(3), MPI_ERR)
 CALL MPI_IRECV(g_bw_rcv, yedge, MPI_DOUBLE_PRECISION, bw, TAG4, MPI_COMM_VGRID, MPI_REQ(4), MPI_ERR)
 count = 1
 DO j = yl, yu
   g_te_snd(count) = g(xug,j,zug,11,nxt)
   g_tw_snd(count) = g(xlg,j,zug,14,nxt)
   g_bw_snd(count) = g(xlg,j,zlg,12,nxt)
   g_be_snd(count) = g(xug,j,zlg,13,nxt)
   count = count + 1
 END DO
 CALL MPI_SEND(g_bw_snd, yedge, MPI_DOUBLE_PRECISION, bw, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(g_be_snd, yedge, MPI_DOUBLE_PRECISION, be, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(g_tw_snd, yedge, MPI_DOUBLE_PRECISION, tw, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(g_te_snd, yedge, MPI_DOUBLE_PRECISION, te, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(1), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(2), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(3), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(4), status, MPI_ERR)

! Edges along Z
 CALL MPI_IRECV(g_ne_rcv, zedge, MPI_DOUBLE_PRECISION, ne, TAG1, MPI_COMM_VGRID, MPI_REQ(1), MPI_ERR)
 CALL MPI_IRECV(g_nw_rcv, zedge, MPI_DOUBLE_PRECISION, nw, TAG2, MPI_COMM_VGRID, MPI_REQ(2), MPI_ERR)
 CALL MPI_IRECV(g_se_rcv, zedge, MPI_DOUBLE_PRECISION, se, TAG3, MPI_COMM_VGRID, MPI_REQ(3), MPI_ERR)
 CALL MPI_IRECV(g_sw_rcv, zedge, MPI_DOUBLE_PRECISION, sw, TAG4, MPI_COMM_VGRID, MPI_REQ(4), MPI_ERR)
 count = 1
 DO k = zl, zu
   g_ne_snd(count) = g(xug,yug,k, 7,nxt)
   g_nw_snd(count) = g(xlg,yug,k,10,nxt)
   g_sw_snd(count) = g(xlg,ylg,k, 8,nxt)
   g_se_snd(count) = g(xug,ylg,k, 9,nxt)
   count = count + 1
 END DO
 CALL MPI_SEND(g_sw_snd, zedge, MPI_DOUBLE_PRECISION, sw, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(g_se_snd, zedge, MPI_DOUBLE_PRECISION, se, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(g_nw_snd, zedge, MPI_DOUBLE_PRECISION, nw, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(g_ne_snd, zedge, MPI_DOUBLE_PRECISION, ne, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(1), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(2), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(3), status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ(4), status, MPI_ERR)

!-------- Update outward f components for the ghosts ---------------------------
! X direction
 count = 1
 DO k = zl, zu
   DO j = yl, yu
     f(xug,j,k,1,now) = f_east_rcv(count)
     f(xlg,j,k,2,now) = f_west_rcv(count)
     count = count + 1
   END DO
 END DO

! Y direction
 count = 1
 DO k = zl, zu
   DO i = xl, xu
     f(i,yug,k,3,now) = f_north_rcv(count)
     f(i,ylg,k,4,now) = f_south_rcv(count)
     count = count + 1
   END DO
 END DO

! Z direction
 count = 1
 DO j = yl, yu
   DO i = xl, xu
     f(i,j,zug,5,now) = f_top_rcv(count)
     f(i,j,zlg,6,now) = f_bot_rcv(count)
     count = count + 1
   END DO
 END DO

 RETURN
 END SUBROUTINE PostCollision

