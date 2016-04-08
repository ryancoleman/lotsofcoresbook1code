!-------------------------------------------------------------------------------
! Subroutine : MemAlloc
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Allocate (FLAG == 1) or deallocate (FLAG != 1) memory for common arrays.
!> @details
!! If FLAG is 1 then memory is allocated for the following common arrays:
!!
!! - \b phi (order parameter)
!! - \b ni  (near neighbors)
!! - \b f   (order parameter distribution function)
!! - \b g   (momentum distribution function)
!!
!! and the MPI buffer arrays required by the parallel D3Q7/D3Q19 Zheng-Shu-Chew
!! multiphase LBM. Memory is deallocated when FLAG is different from one.

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

 SUBROUTINE MemAlloc(FLAG)

!  Common Variables
 USE Domain
 USE FluidParams, ONLY : phi
 USE LBMParams,   ONLY : f, fdim, g, gdim
 USE MPIParams
 IMPLICIT NONE

 !   Input/Output Parameters
 INTEGER, INTENT (IN) :: FLAG

 IF(FLAG == 1) THEN

! Near neighbors and order parameter arrays
 ALLOCATE(  ni(xlg:xug,ylg:yug,zlg:zug,1:6) )
 ALLOCATE( phi(xlg:xug,ylg:yug,zlg:zug) )

! Distribution function arrays
 ALLOCATE( f(xlg:xug,ylg:yug,zlg:zug,0:fdim,0:1) )
 ALLOCATE( g(xlg:xug,ylg:yug,zlg:zug,0:gdim,0:1) )

! Allocation for data exchange in the x direction
   xsize  = (yu - yl + 1)*(zu - zl + 1)
   xsize2 = xsize + xsize
   xsize3 = xsize2 + xsize
   xsize4 = xsize3 + xsize
   xsize5 = xsize4 + xsize

   ALLOCATE( phi_west_snd(1:xsize) )
   ALLOCATE( phi_east_snd(1:xsize) )
   ALLOCATE( phi_west_rcv(1:xsize) )
   ALLOCATE( phi_east_rcv(1:xsize) )
   
   ALLOCATE( f_west_snd(1:xsize) )
   ALLOCATE( f_east_snd(1:xsize) )
   ALLOCATE( f_west_rcv(1:xsize) )
   ALLOCATE( f_east_rcv(1:xsize) )

   ALLOCATE( g_west_snd(1:xsize5) )
   ALLOCATE( g_east_snd(1:xsize5) )
   ALLOCATE( g_west_rcv(1:xsize5) )
   ALLOCATE( g_east_rcv(1:xsize5) )

! Allocation for data exchange in the y direction
   ysize  = (xu - xl + 1)*(zu - zl + 1)
   ysize2 = ysize + ysize
   ysize3 = ysize2 + ysize
   ysize4 = ysize3 + ysize
   ysize5 = ysize4 + ysize

   ALLOCATE( phi_south_snd(1:ysize) )
   ALLOCATE( phi_north_snd(1:ysize) )
   ALLOCATE( phi_south_rcv(1:ysize) )
   ALLOCATE( phi_north_rcv(1:ysize) )

   ALLOCATE( f_south_snd(1:ysize) )
   ALLOCATE( f_north_snd(1:ysize) )
   ALLOCATE( f_south_rcv(1:ysize) )
   ALLOCATE( f_north_rcv(1:ysize) )

   ALLOCATE( g_south_snd(1:ysize5) )
   ALLOCATE( g_north_snd(1:ysize5) )
   ALLOCATE( g_south_rcv(1:ysize5) )
   ALLOCATE( g_north_rcv(1:ysize5) )

! Allocation for data exchanged in the y direction
   zsize  = (xu - xl + 1)*(yu - yl + 1)
   zsize2 = zsize + zsize
   zsize3 = zsize2 + zsize
   zsize4 = zsize3 + zsize
   zsize5 = zsize4 + zsize

   ALLOCATE( phi_bot_snd(1:zsize) )
   ALLOCATE( phi_top_snd(1:zsize) )
   ALLOCATE( phi_bot_rcv(1:zsize) )
   ALLOCATE( phi_top_rcv(1:zsize) )

   ALLOCATE( f_bot_snd(1:zsize) )
   ALLOCATE( f_bot_rcv(1:zsize) )
   ALLOCATE( f_top_snd(1:zsize) )
   ALLOCATE( f_top_rcv(1:zsize) )

   ALLOCATE( g_bot_snd(1:zsize5) )
   ALLOCATE( g_bot_rcv(1:zsize5) )
   ALLOCATE( g_top_snd(1:zsize5) )
   ALLOCATE( g_top_rcv(1:zsize5) )

! Allocation for data exchange in edges along the x direction
   xedge = xu - xl + 1

   ALLOCATE( phi_tn_snd(1:xedge) )
   ALLOCATE( phi_ts_snd(1:xedge) )
   ALLOCATE( phi_bn_snd(1:xedge) )
   ALLOCATE( phi_bs_snd(1:xedge) )
   ALLOCATE( phi_tn_rcv(1:xedge) )
   ALLOCATE( phi_ts_rcv(1:xedge) )
   ALLOCATE( phi_bn_rcv(1:xedge) )
   ALLOCATE( phi_bs_rcv(1:xedge) )
   
   ALLOCATE( g_tn_snd(1:xedge) )
   ALLOCATE( g_ts_snd(1:xedge) )
   ALLOCATE( g_bn_snd(1:xedge) )
   ALLOCATE( g_bs_snd(1:xedge) )
   ALLOCATE( g_tn_rcv(1:xedge) )
   ALLOCATE( g_ts_rcv(1:xedge) )
   ALLOCATE( g_bn_rcv(1:xedge) )
   ALLOCATE( g_bs_rcv(1:xedge) )

! Allocation for data exchange in edges along the y direction
   yedge = yu - yl + 1

   ALLOCATE( phi_te_snd(1:yedge) )
   ALLOCATE( phi_tw_snd(1:yedge) )
   ALLOCATE( phi_be_snd(1:yedge) )
   ALLOCATE( phi_bw_snd(1:yedge) )
   ALLOCATE( phi_te_rcv(1:yedge) )
   ALLOCATE( phi_tw_rcv(1:yedge) )
   ALLOCATE( phi_be_rcv(1:yedge) )
   ALLOCATE( phi_bw_rcv(1:yedge) )

   ALLOCATE( g_te_snd(1:yedge) )
   ALLOCATE( g_tw_snd(1:yedge) )
   ALLOCATE( g_be_snd(1:yedge) )
   ALLOCATE( g_bw_snd(1:yedge) )
   ALLOCATE( g_te_rcv(1:yedge) )
   ALLOCATE( g_tw_rcv(1:yedge) )
   ALLOCATE( g_be_rcv(1:yedge) )
   ALLOCATE( g_bw_rcv(1:yedge) )

! Allocation for data exchange in edges along the z direction
   zedge = zu - zl + 1

   ALLOCATE( phi_ne_snd(1:zedge) )
   ALLOCATE( phi_nw_snd(1:zedge) )
   ALLOCATE( phi_se_snd(1:zedge) )
   ALLOCATE( phi_sw_snd(1:zedge) )
   ALLOCATE( phi_ne_rcv(1:zedge) )
   ALLOCATE( phi_nw_rcv(1:zedge) )
   ALLOCATE( phi_se_rcv(1:zedge) )
   ALLOCATE( phi_sw_rcv(1:zedge) )

   ALLOCATE( g_ne_snd(1:zedge) )
   ALLOCATE( g_nw_snd(1:zedge) )
   ALLOCATE( g_se_snd(1:zedge) )
   ALLOCATE( g_sw_snd(1:zedge) )
   ALLOCATE( g_ne_rcv(1:zedge) )
   ALLOCATE( g_nw_rcv(1:zedge) )
   ALLOCATE( g_se_rcv(1:zedge) )
   ALLOCATE( g_sw_rcv(1:zedge) )

 ELSE
   DEALLOCATE( ni )
   DEALLOCATE( phi )
   DEALLOCATE( f )
   DEALLOCATE( g )
   DEALLOCATE( f_west_snd )
   DEALLOCATE( f_west_rcv )
   DEALLOCATE( f_east_snd )
   DEALLOCATE( f_east_rcv )
   DEALLOCATE( g_west_snd )
   DEALLOCATE( g_east_snd )
   DEALLOCATE( g_west_rcv )
   DEALLOCATE( g_east_rcv )
   DEALLOCATE( phi_west_snd )
   DEALLOCATE( phi_west_rcv )
   DEALLOCATE( phi_east_snd )
   DEALLOCATE( phi_east_rcv )
   DEALLOCATE( f_south_snd )
   DEALLOCATE( f_south_rcv )
   DEALLOCATE( f_north_snd )
   DEALLOCATE( f_north_rcv )
   DEALLOCATE( g_south_snd )
   DEALLOCATE( g_north_snd )
   DEALLOCATE( g_south_rcv )
   DEALLOCATE( g_north_rcv )
   DEALLOCATE( phi_south_snd )
   DEALLOCATE( phi_south_rcv )
   DEALLOCATE( phi_north_snd )
   DEALLOCATE( phi_north_rcv )
   DEALLOCATE( f_bot_snd )
   DEALLOCATE( f_bot_rcv )
   DEALLOCATE( f_top_snd )
   DEALLOCATE( f_top_rcv )
   DEALLOCATE( g_bot_snd )
   DEALLOCATE( g_bot_rcv )
   DEALLOCATE( g_top_snd )
   DEALLOCATE( g_top_rcv )
   DEALLOCATE( phi_bot_snd )
   DEALLOCATE( phi_bot_rcv )
   DEALLOCATE( phi_top_snd )
   DEALLOCATE( phi_top_rcv )
   DEALLOCATE( phi_tn_snd )
   DEALLOCATE( phi_ts_snd )
   DEALLOCATE( phi_te_snd )
   DEALLOCATE( phi_tw_snd )
   DEALLOCATE( phi_bn_snd )
   DEALLOCATE( phi_bs_snd )
   DEALLOCATE( phi_be_snd )
   DEALLOCATE( phi_bw_snd )
   DEALLOCATE( phi_ne_snd )
   DEALLOCATE( phi_nw_snd )
   DEALLOCATE( phi_se_snd )
   DEALLOCATE( phi_sw_snd ) 
   DEALLOCATE( phi_tn_rcv )
   DEALLOCATE( phi_ts_rcv )
   DEALLOCATE( phi_te_rcv )
   DEALLOCATE( phi_tw_rcv )
   DEALLOCATE( phi_bn_rcv )
   DEALLOCATE( phi_bs_rcv )
   DEALLOCATE( phi_be_rcv )
   DEALLOCATE( phi_bw_rcv )
   DEALLOCATE( phi_ne_rcv )
   DEALLOCATE( phi_nw_rcv )
   DEALLOCATE( phi_se_rcv )
   DEALLOCATE( phi_sw_rcv )

 END IF

 RETURN
 END SUBROUTINE MemAlloc
