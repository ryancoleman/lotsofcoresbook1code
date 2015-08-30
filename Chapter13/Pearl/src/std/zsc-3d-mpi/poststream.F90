!-------------------------------------------------------------------------------
! Subroutine : PostStream
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Update f in nodes at processor boundaries using MPI.
!> @details
!! Copy outward f values from the ghosts into inward f values in the real nodes,
!! then update order parameter phi and the values values for g and f in the
!! ghost layer.
!!
!! MPI exchange is done using non-blocking receive and blocking send pairs
!! (completed by a wait statement) in alternate directions to avoid deadlock.
!! The componets of g and f that need to be exchanged are packed into buffer
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

 SUBROUTINE PostStream

!  Common Variables
 USE Domain
 USE FluidParams, ONLY : phi
 USE LBMParams,   ONLY : f, g
 USE MPIParams
 USE MPI
 IMPLICIT NONE

!  Local Variables
 INTEGER :: count, i, j, k
 INTEGER :: MPI_ERR, MPI_REQ1, MPI_REQ2, MPI_REQ3, MPI_REQ4
 INTEGER :: status(MPI_STATUS_SIZE)

!----------Exchange inward pointing f values -----------------------------------
! X direction
 CALL MPI_IRECV(f_east_rcv, xsize , MPI_DOUBLE_PRECISION, east, TAG1, MPI_COMM_VGRID, MPI_REQ1, MPI_ERR)
 CALL MPI_IRECV(f_west_rcv, xsize , MPI_DOUBLE_PRECISION, west, TAG2, MPI_COMM_VGRID, MPI_REQ2, MPI_ERR)
 count = 1
 DO k = zl, zu
   DO j = yl, yu
     f_west_snd(count) = f(xlg,j,k,2,nxt)
     f_east_snd(count) = f(xug,j,k,1,nxt)
     count = count + 1
   END DO
 END DO
 CALL MPI_SEND(f_west_snd, xsize , MPI_DOUBLE_PRECISION, west, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(f_east_snd, xsize , MPI_DOUBLE_PRECISION, east, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ1, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ2, status, MPI_ERR)

! Y direction
 CALL MPI_IRECV(f_north_rcv, ysize, MPI_DOUBLE_PRECISION, north, TAG1, MPI_COMM_VGRID, MPI_REQ1, MPI_ERR)
 CALL MPI_IRECV(f_south_rcv, ysize, MPI_DOUBLE_PRECISION, south, TAG2, MPI_COMM_VGRID, MPI_REQ2, MPI_ERR)
 count = 1
 DO k = zl, zu
   DO i = xl, xu
     f_south_snd(count) = f(i,ylg,k,4,nxt)
     f_north_snd(count) = f(i,yug,k,3,nxt)
     count = count + 1
   END DO
 END DO
 CALL MPI_SEND(f_south_snd, ysize, MPI_DOUBLE_PRECISION, south, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(f_north_snd, ysize, MPI_DOUBLE_PRECISION, north, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ1, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ2, status, MPI_ERR)

! Z direction
 CALL MPI_IRECV(f_top_rcv, zsize, MPI_DOUBLE_PRECISION, top, TAG1, MPI_COMM_VGRID, MPI_REQ1, MPI_ERR)
 CALL MPI_IRECV(f_bot_rcv, zsize, MPI_DOUBLE_PRECISION, bot, TAG2, MPI_COMM_VGRID, MPI_REQ2, MPI_ERR)
 count = 1
 DO j = yl, yu
   DO i = xl, xu
     f_bot_snd(count) = f(i,j,zlg,6,nxt)
     f_top_snd(count)    = f(i,j,zug,5,nxt)
     count = count + 1
   END DO
 END DO
 CALL MPI_SEND(f_bot_snd, zsize, MPI_DOUBLE_PRECISION, bot, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(f_top_snd, zsize, MPI_DOUBLE_PRECISION, top, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ1, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ2, status, MPI_ERR)


!-------- Update inward values of f and g in real nodes ------------------------
! X direction
 count = 1
 DO k = zl, zu
   DO j = yl, yu
     f(xu,j,k, 2,nxt) = f_east_rcv(count)
     g(xu,j,k, 2,nxt) = g_east_rcv(count)
     g(xu,j,k, 8,nxt) = g_east_rcv(count + xsize)
     g(xu,j,k,10,nxt) = g_east_rcv(count + xsize2)
     g(xu,j,k,12,nxt) = g_east_rcv(count + xsize3)
     g(xu,j,k,14,nxt) = g_east_rcv(count + xsize4)

     f(xl,j,k, 1,nxt)  = f_west_rcv(count)
     g(xl,j,k, 1,nxt)  = g_west_rcv(count)
     g(xl,j,k, 7,nxt)  = g_west_rcv(count + xsize)
     g(xl,j,k, 9,nxt)  = g_west_rcv(count + xsize2)
     g(xl,j,k,11,nxt)  = g_west_rcv(count + xsize3)
     g(xl,j,k,13,nxt)  = g_west_rcv(count + xsize4)
     count = count + 1
   END DO
 END DO

! Y direction
 count = 1
 DO k = zl, zu
   DO i = xl, xu
     f(i,yu,k, 4,nxt) = f_north_rcv(count)
     g(i,yu,k, 4,nxt) = g_north_rcv(count)
     g(i,yu,k, 8,nxt) = g_north_rcv(count + ysize)
     g(i,yu,k, 9,nxt) = g_north_rcv(count + ysize2)
     g(i,yu,k,16,nxt) = g_north_rcv(count + ysize3)
     g(i,yu,k,18,nxt) = g_north_rcv(count + ysize4)

     f(i,yl,k, 3,nxt)  = f_south_rcv(count)
     g(i,yl,k, 3,nxt)  = g_south_rcv(count)
     g(i,yl,k, 7,nxt)  = g_south_rcv(count + ysize)
     g(i,yl,k,10,nxt)  = g_south_rcv(count + ysize2)
     g(i,yl,k,15,nxt)  = g_south_rcv(count + ysize3)
     g(i,yl,k,17,nxt)  = g_south_rcv(count + ysize4)
     count = count + 1
   END DO
 END DO

! Z direction
 count = 1
 DO j = yl, yu
   DO i = xl, xu
     f(i,j,zu, 6,nxt) = f_top_rcv(count)
     g(i,j,zu, 6,nxt) = g_top_rcv(count)
     g(i,j,zu,12,nxt) = g_top_rcv(count + zsize)
     g(i,j,zu,13,nxt) = g_top_rcv(count + zsize2)
     g(i,j,zu,16,nxt) = g_top_rcv(count + zsize3)
     g(i,j,zu,17,nxt) = g_top_rcv(count + zsize4)

     f(i,j,zl, 5,nxt)  = f_bot_rcv(count)
     g(i,j,zl, 5,nxt)  = g_bot_rcv(count)
     g(i,j,zl,11,nxt)  = g_bot_rcv(count + zsize)
     g(i,j,zl,14,nxt)  = g_bot_rcv(count + zsize2)
     g(i,j,zl,15,nxt)  = g_bot_rcv(count + zsize3)
     g(i,j,zl,18,nxt)  = g_bot_rcv(count + zsize4)
     count = count + 1
   END DO
 END DO

!  Edges along x
 count = 1
 DO i = xl, xu
   g(i,yu,zu,16,nxt) = g_tn_rcv(count)
   g(i,yl,zu,17,nxt) = g_ts_rcv(count)
   g(i,yl,zl,15,nxt) = g_bs_rcv(count)
   g(i,yu,zl,18,nxt) = g_bn_rcv(count)
   count = count + 1
 END DO

!  Edges along y
 count = 1
 DO j = yl, yu
   g(xu,j,zu,12,nxt) = g_te_rcv(count)
   g(xl,j,zu,13,nxt) = g_tw_rcv(count)
   g(xl,j,zl,11,nxt) = g_bw_rcv(count)
   g(xu,j,zl,14,nxt) = g_be_rcv(count)
   count = count + 1
 END DO

!  Edges along z
 count = 1
 DO k = zl, zu
   g(xu,yu,k, 8,nxt) = g_ne_rcv(count)
   g(xl,yu,k, 9,nxt) = g_nw_rcv(count)
   g(xl,yl,k, 7,nxt) = g_sw_rcv(count)
   g(xu,yl,k,10,nxt) = g_se_rcv(count)
   count = count + 1
 END DO

!-------- Update order parameter for all ghost nodes ---------------------------
! X direction
 CALL MPI_IRECV(phi_east_rcv, xsize, MPI_DOUBLE_PRECISION, east, TAG1, MPI_COMM_VGRID, MPI_REQ1, MPI_ERR)
 CALL MPI_IRECV(phi_west_rcv, xsize, MPI_DOUBLE_PRECISION, west, TAG2, MPI_COMM_VGRID, MPI_REQ2, MPI_ERR)
 count = 1
 DO k = zl, zu
   DO j = yl, yu
     phi_west_snd(count) = SUM(f(xl,j,k,:,nxt))
     phi_east_snd(count) = SUM(f(xu,j,k,:,nxt))
     count = count + 1
   END DO
 END DO
 CALL MPI_SEND(phi_west_snd, xsize, MPI_DOUBLE_PRECISION, west, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(phi_east_snd, xsize, MPI_DOUBLE_PRECISION, east, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ1, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ2, status, MPI_ERR)

! Y direction
 CALL MPI_IRECV(phi_north_rcv, ysize, MPI_DOUBLE_PRECISION, north, TAG1, MPI_COMM_VGRID, MPI_REQ1, MPI_ERR)
 CALL MPI_IRECV(phi_south_rcv, ysize, MPI_DOUBLE_PRECISION, south, TAG2, MPI_COMM_VGRID, MPI_REQ2, MPI_ERR)
 count = 1
 DO k = zl, zu
   DO i = xl, xu
     phi_south_snd(count) = SUM(f(i,yl,k,:,nxt))
     phi_north_snd(count) = SUM(f(i,yu,k,:,nxt))
     count = count + 1
   END DO
 END DO
 CALL MPI_SEND(phi_south_snd, ysize, MPI_DOUBLE_PRECISION, south, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(phi_north_snd, ysize, MPI_DOUBLE_PRECISION, north, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ1, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ2, status, MPI_ERR)

! Z direction
 CALL MPI_IRECV(phi_top_rcv, zsize, MPI_DOUBLE_PRECISION, top, TAG1, MPI_COMM_VGRID, MPI_REQ1, MPI_ERR)
 CALL MPI_IRECV(phi_bot_rcv, zsize, MPI_DOUBLE_PRECISION, bot, TAG2, MPI_COMM_VGRID, MPI_REQ2, MPI_ERR)
 count = 1
 DO j = yl, yu
   DO i = xl, xu
     phi_bot_snd(count) = SUM(f(i,j,zl,:,nxt))
     phi_top_snd(count)    = SUM(f(i,j,zu,:,nxt))
     count = count + 1
   END DO
 END DO
 CALL MPI_SEND(phi_bot_snd, zsize, MPI_DOUBLE_PRECISION, bot, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(phi_top_snd, zsize, MPI_DOUBLE_PRECISION, top, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ1, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ2, status, MPI_ERR)

! Edges along x
 CALL MPI_IRECV(phi_tn_rcv, xedge, MPI_DOUBLE_PRECISION, tn, TAG1, MPI_COMM_VGRID, MPI_REQ1, MPI_ERR)
 CALL MPI_IRECV(phi_bs_rcv, xedge, MPI_DOUBLE_PRECISION, bs, TAG2, MPI_COMM_VGRID, MPI_REQ2, MPI_ERR)
 CALL MPI_IRECV(phi_bn_rcv, xedge, MPI_DOUBLE_PRECISION, bn, TAG3, MPI_COMM_VGRID, MPI_REQ3, MPI_ERR)
 CALL MPI_IRECV(phi_ts_rcv, xedge, MPI_DOUBLE_PRECISION, ts, TAG4, MPI_COMM_VGRID, MPI_REQ4, MPI_ERR)
 count = 1
 DO i = xl, xu
   phi_tn_snd(count) = SUM(f(i,yu,zu,:,nxt))
   phi_ts_snd(count) = SUM(f(i,yl,zu,:,nxt))
   phi_bs_snd(count) = SUM(f(i,yl,zl,:,nxt))
   phi_bn_snd(count) = SUM(f(i,yu,zl,:,nxt))
   count = count + 1
 END DO
 CALL MPI_SEND(phi_bs_snd, xedge, MPI_DOUBLE_PRECISION, bs, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(phi_tn_snd, xedge, MPI_DOUBLE_PRECISION, tn, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(phi_ts_snd, xedge, MPI_DOUBLE_PRECISION, ts, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(phi_bn_snd, xedge, MPI_DOUBLE_PRECISION, bn, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ1, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ2, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ3, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ4, status, MPI_ERR)
 
! Edges along y
 CALL MPI_IRECV(phi_te_rcv, yedge, MPI_DOUBLE_PRECISION, te, TAG1, MPI_COMM_VGRID, MPI_REQ1, MPI_ERR)
 CALL MPI_IRECV(phi_bw_rcv, yedge, MPI_DOUBLE_PRECISION, bw, TAG2, MPI_COMM_VGRID, MPI_REQ2, MPI_ERR)
 CALL MPI_IRECV(phi_be_rcv, yedge, MPI_DOUBLE_PRECISION, be, TAG3, MPI_COMM_VGRID, MPI_REQ3, MPI_ERR)
 CALL MPI_IRECV(phi_tw_rcv, yedge, MPI_DOUBLE_PRECISION, tw, TAG4, MPI_COMM_VGRID, MPI_REQ4, MPI_ERR)
 count = 1
 DO j = yl, yu
   phi_te_snd(count) = SUM(f(xu,j,zu,:,nxt))
   phi_tw_snd(count) = SUM(f(xl,j,zu,:,nxt))
   phi_bw_snd(count) = SUM(f(xl,j,zl,:,nxt))
   phi_be_snd(count) = SUM(f(xu,j,zl,:,nxt))
   count = count + 1
 END DO
 CALL MPI_SEND(phi_bw_snd, yedge, MPI_DOUBLE_PRECISION, bw, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(phi_te_snd, yedge, MPI_DOUBLE_PRECISION, te, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(phi_tw_snd, yedge, MPI_DOUBLE_PRECISION, tw, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(phi_be_snd, yedge, MPI_DOUBLE_PRECISION, be, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ1, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ2, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ3, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ4, status, MPI_ERR)
 
! Edges along z
 CALL MPI_IRECV(phi_ne_rcv, zedge, MPI_DOUBLE_PRECISION, ne, TAG1, MPI_COMM_VGRID, MPI_REQ1, MPI_ERR)
 CALL MPI_IRECV(phi_sw_rcv, zedge, MPI_DOUBLE_PRECISION, sw, TAG2, MPI_COMM_VGRID, MPI_REQ2, MPI_ERR)
 CALL MPI_IRECV(phi_se_rcv, zedge, MPI_DOUBLE_PRECISION, se, TAG3, MPI_COMM_VGRID, MPI_REQ3, MPI_ERR)
 CALL MPI_IRECV(phi_nw_rcv, zedge, MPI_DOUBLE_PRECISION, nw, TAG4, MPI_COMM_VGRID, MPI_REQ4, MPI_ERR)
 count = 1
 DO k = zl, zu
   phi_ne_snd(count) = SUM(f(xu,yu,k,:,nxt))
   phi_nw_snd(count) = SUM(f(xl,yu,k,:,nxt))
   phi_sw_snd(count) = SUM(f(xl,yl,k,:,nxt))
   phi_se_snd(count) = SUM(f(xu,yl,k,:,nxt))
   count = count + 1
 END DO
 CALL MPI_SEND(phi_sw_snd, zedge, MPI_DOUBLE_PRECISION, sw, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(phi_ne_snd, zedge, MPI_DOUBLE_PRECISION, ne, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(phi_nw_snd, zedge, MPI_DOUBLE_PRECISION, nw, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_SEND(phi_se_snd, zedge, MPI_DOUBLE_PRECISION, se, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ1, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ2, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ3, status, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ4, status, MPI_ERR)

!-------- Update order parameter values ----------------------------------------
! X direction
 count = 1
 DO k = zl, zu
   DO j = yl, yu
     phi(xug,j,k) = phi_east_rcv(count)
     phi(xlg,j,k) = phi_west_rcv(count)
     count = count + 1
   END DO
 END DO

! Y direction
 count = 1
 DO k = zl, zu
   DO i = xl, xu
     phi(i,yug,k) = phi_north_rcv(count)
     phi(i,ylg,k) = phi_south_rcv(count)
     count = count + 1
   END DO
 END DO

! Z direction
 count = 1
 DO j = yl, yu
   DO i = xl, xu
     phi(i,j,zug) = phi_top_rcv(count)
     phi(i,j,zlg) = phi_bot_rcv(count)
     count = count + 1
   END DO
 END DO

! Edges along x
 count = 1
 DO i = xl, xu
   phi(i,yug,zug) = phi_tn_rcv(count)
   phi(i,ylg,zug) = phi_ts_rcv(count)
   phi(i,ylg,zlg) = phi_bs_rcv(count)
   phi(i,yug,zlg) = phi_bn_rcv(count)
   count = count + 1
 END DO

! Edges along y
 count = 1
 DO j = yl, yu
   phi(xug,j,zug) = phi_te_rcv(count)
   phi(xlg,j,zug) = phi_tw_rcv(count)
   phi(xlg,j,zlg) = phi_bw_rcv(count)
   phi(xug,j,zlg) = phi_be_rcv(count)
   count = count + 1
 END DO

! Edges along z
 count = 1
 DO k = zl, zu
   phi(xug,yug,k) = phi_ne_rcv(count)
   phi(xlg,yug,k) = phi_nw_rcv(count)
   phi(xlg,ylg,k) = phi_sw_rcv(count)
   phi(xug,ylg,k) = phi_se_rcv(count)
   count = count + 1
 END DO

 RETURN
 END SUBROUTINE PostStream
