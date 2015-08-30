!-------------------------------------------------------------------------------
! Program    : LBS3D-DEVEL-OFFLOAD
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Driver for the offload OpenMP implementation of the Zheng-Shu-Chew
! multiphase LBM using D3Q7/D3Q19 discretization and periodic boundary
! conditions, including the gravitational force. For details:
!
! Journal of Computational Physics 218: 353-371, 2006.
!
! The average velocity, mass conservation factor, effective radius of the drop,
! pressure difference between the inside and the outside of the drop and the
! error with respect to the analytical value given by Laplace's equation are
! written to file "stats.out"
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

#define ALLOC   alloc_if(.TRUE.)
#define REUSE   alloc_if(.FALSE.)
#define FREE    free_if(.TRUE.)
#define KEEP    free_if(.FALSE.)

 PROGRAM main

! Common Variables
 USE Timers
 USE Domain
 USE FluidParams
 USE LBMParams
 USE MIC_LIB
 USE OMP_LIB
 IMPLICIT NONE

! Define subroutines that will be offloaded
!DIR$ ATTRIBUTES OFFLOAD: mic :: InitMIC, CollisionMIC, StreamMIC, PostCollisionMIC, PostStreamMIC, PackMIC, PackMIC_f, UpdatePhiMIC

 INTEGER :: SYNC

! Read input parameters, broadcast, and create virtual CPU grid
 tStart = omp_get_wtime()

 micNum = OFFLOAD_NUMBER_OF_DEVICES()
 micID  = OFFLOAD_GET_DEVICE_NUMBER()

!DIR$ OMP OFFLOAD TARGET(mic:0) INOUT( micThreads )
!$OMP PARALLEL
!$OMP MASTER
 micThreads = omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL
 WRITE(*,*)'Number of MIC devices  :',micNUM
 WRITE(*,*)'Executing on device ID :',micID
 WRITE(*,*)'Number of OMP Threads  :',micThreads

 CALL Parameters
 now = 0
 nxt = NG
 now_mic = 0
 nxt_mic = NG_MIC

! Allocate memory for common arrays and initialize
 CALL MemAlloc(1)

!-------------------------------------------------------------------------------
! Let's move all this data to the MIC and make sure it is not cleared on exit
! Transfer is synchronous unless the SIGNAL(tag) specifier is added
! Scalars do not require an alloc/free statement
!-------------------------------------------------------------------------------
tDataStart = omp_get_wtime()
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( f0_mic, f1_mic, f2_mic : ALLOC KEEP )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( f3_mic,f4_mic,f5_mic, f6_mic : ALLOC KEEP )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( g0_mic, g1_mic, g2_mic, g3_mic : ALLOC KEEP )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( g4_mic, g5_mic, g6_mic : ALLOC KEEP )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( g7_mic, g8_mic, g9_mic, g10_mic : ALLOC KEEP )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( g11_mic, g12_mic, g13_mic, g14_mic : ALLOC KEEP )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( g15_mic, g16_mic, g17_mic, g18_mic : ALLOC KEEP )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( phi_mic, lapPhi_mic, gradPhiX_mic : ALLOC KEEP )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( gradPhiY_mic,gradPhiZ_mic,bubbles : ALLOC KEEP )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( f_buff_mic, f_buff_cpu : ALLOC KEEP )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( buff_mic, buff_cpu : ALLOC KEEP )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( buff_mic_edge_x, buff_cpu_edge_x : ALLOC KEEP )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( buff_mic_edge_z, buff_cpu_edge_z : ALLOC KEEP )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( buff_mic_phi, buff_cpu_phi : ALLOC KEEP )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( buff_mic_phi_edge_x, buff_cpu_phi_edge_x : ALLOC KEEP )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( buff_mic_phi_edge_z, buff_cpu_phi_edge_z : ALLOC KEEP )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( inv6, inv12, now_mic, nxt_mic, iStep, RelaxStep )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( NX, NXG, NY, NYG, NZ, NZG, NG, MaxStep )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( NY_MIC, NYG_MIC, NXYG_MIC, NG_MIC )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( alpha4, phiStar2, kappa, gamma, invTauPhi )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( invEta2, invTauPhi1, Eg0n, Eg1n, Eg2n, invCs_sq )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( eta, eta2, invTauRhoOne, EgC0n, EgC1n, EgC2n )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( pConv, eps, iStep, nBubbles, IntWidth )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( rhoH, rhoL, phiStar, tCall, xl, yl, zl )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( NPX, NPY, NPZ, xmin, xmax, ymin, ymax, zmin, zmax )
tDataEnd = omp_get_wtime()
tData = tDataEnd - tDataStart

! Initialize data
!DIR$ OFFLOAD BEGIN TARGET(MIC:0) SIGNAL(SYNC)
 CALL InitMIC
!DIR$ END OFFLOAD
 CALL Init
!DIR$ OFFLOAD_WAIT TARGET(MIC:0) WAIT(SYNC)

! ***** WARNING *****
! RelaxStats Needs to be run on both the mic and the host...
! Save initialized data
! ***** WARNING *****
 CALL RelaxStats
#ifdef LBMIO
 CALL VtkSave
#endif

! Interface relaxation loop
 grav = 0.D0
 tCurrent = omp_get_wtime()
 tSetup = tCurrent - tStart
! Main Iteration Loop
!-------------------------------------------------------------------------------
! Let's offload. We need to get the f values back for the final log save
!-------------------------------------------------------------------------------
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( grav )
 DO iStep = 1, RelaxStep

!DIR$ OFFLOAD BEGIN TARGET(MIC:0) IN(now_mic,nxt_mic)           &
                                  OUT(f_buff_mic : REUSE KEEP ) &
                                  SIGNAL(SYNC)
   CALL CollisionMIC
   CALL PackMIC_F
!DIR$ END OFFLOAD

   CALL Collision
   CALL PackCPU_F
!DIR$ OFFLOAD_WAIT TARGET(mic:0) WAIT(SYNC)

!DIR$ OFFLOAD BEGIN TARGET(MIC:0) IN(now_mic,nxt_mic)                &
                                  IN(f_buff_cpu       : REUSE KEEP ) &
                                  OUT(buff_mic        : REUSE KEEP ) &
                                  OUT(buff_mic_edge_x : REUSE KEEP ) &
                                  OUT(buff_mic_edge_z : REUSE KEEP ) &
                                  SIGNAL(SYNC)
   CALL PostCollisionMIC
   CALL StreamMIC
   CALL PackMIC
!DIR$ END OFFLOAD

   CALL PostCollision
   CALL Stream
   CALL PackCPU
!DIR$ OFFLOAD_WAIT TARGET(mic:0) WAIT(SYNC)

!DIR$ OFFLOAD BEGIN TARGET(mic:0) IN(nxt_mic)                             &
                                  OUT(buff_mic_phi        : REUSE KEEP )  &
                                  OUT(buff_mic_phi_edge_x : REUSE KEEP )  &
                                  OUT(buff_mic_phi_edge_z : REUSE KEEP )  &
                                  SIGNAL(SYNC)
    CALL UpdatePhiMIC
!DIR$ END OFFLOAD

!DIR$ OFFLOAD_WAIT TARGET(mic:0) WAIT(SYNC)
    CALL UpdatePhi

!DIR$ OFFLOAD BEGIN TARGET(mic:0) IN(nxt_mic)                            &
                                  IN(buff_cpu            : REUSE KEEP )  &
                                  IN(buff_cpu_edge_x     : REUSE KEEP )  &
                                  IN(buff_cpu_edge_z     : REUSE KEEP )  &
                                  IN(buff_cpu_phi        : REUSE KEEP )  &
                                  IN(buff_cpu_phi_edge_x : REUSE KEEP )  &
                                  IN(buff_cpu_phi_edge_z : REUSE KEEP )  &
                                  SIGNAL(SYNC)
   CALL PostStreamMIC
!DIR$ END OFFLOAD

   CALL PostStream
!DIR$ OFFLOAD_WAIT TARGET(mic:0) WAIT(SYNC)

   now = NG - now
   nxt = NG - nxt

   now_mic = NG_MIC - now_mic
   nxt_mic = NG_MIC - nxt_mic

#ifdef LBMIO
   IF( MOD(iStep,tStat) == 0 ) CALL RelaxStats
   IF( MOD(iStep,tSave) == 0 ) CALL VtkSave
#endif

   IF( eps < pConv ) EXIT

 END DO

 tCurrent = omp_get_wtime()
 tRelax = tCurrent - tSetup - tStart
 CompletedRelaxSteps = iStep

! Save relaxed configuration data
 STAGE = 1
#ifdef LBMIO
 CALL Stats
 CALL VtkSave
 tCurrent = omp_get_wtime()
 tIO = tCurrent - tRelax - tSetup - tStart
#else
 tIO = 0.D0
#endif

! Main Iteration Loop
 grav = 2.D0*phiStar*gravity
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) IN( grav )
 DO iStep = 1, Maxstep

!DIR$ OFFLOAD BEGIN TARGET(MIC:0) IN(now_mic,nxt_mic)           &
                                  OUT(f_buff_mic : REUSE KEEP ) &
                                  SIGNAL(SYNC)
   CALL CollisionMIC
   CALL PackMIC_F
!DIR$ END OFFLOAD

   CALL Collision
   CALL PackCPU_F
!DIR$ OFFLOAD_WAIT TARGET(mic:0) WAIT(SYNC)

!DIR$ OFFLOAD BEGIN TARGET(MIC:0) IN(now_mic,nxt_mic)                &
                                  IN(f_buff_cpu       : REUSE KEEP ) &
                                  OUT(buff_mic        : REUSE KEEP ) &
                                  OUT(buff_mic_edge_x : REUSE KEEP ) &
                                  OUT(buff_mic_edge_z : REUSE KEEP ) &
                                  SIGNAL(SYNC)
   CALL PostCollisionMIC
   CALL StreamMIC
   CALL PackMIC
!DIR$ END OFFLOAD

   CALL PostCollision
   CALL Stream
   CALL PackCPU
!DIR$ OFFLOAD_WAIT TARGET(mic:0) WAIT(SYNC)

!DIR$ OFFLOAD BEGIN TARGET(mic:0) IN(nxt_mic)                             &
                                  OUT(buff_mic_phi        : REUSE KEEP )  &
                                  OUT(buff_mic_phi_edge_x : REUSE KEEP )  &
                                  OUT(buff_mic_phi_edge_z : REUSE KEEP )  &
                                  SIGNAL(SYNC)
    CALL UpdatePhiMIC
!DIR$ END OFFLOAD

!DIR$ OFFLOAD_WAIT TARGET(mic:0) WAIT(SYNC)
    CALL UpdatePhi

!DIR$ OFFLOAD BEGIN TARGET(mic:0) IN(nxt_mic)                            &
                                  IN(buff_cpu            : REUSE KEEP )  &
                                  IN(buff_cpu_edge_x     : REUSE KEEP )  &
                                  IN(buff_cpu_edge_z     : REUSE KEEP )  &
                                  IN(buff_cpu_phi        : REUSE KEEP )  &
                                  IN(buff_cpu_phi_edge_x : REUSE KEEP )  &
                                  IN(buff_cpu_phi_edge_z : REUSE KEEP )  &
                                  SIGNAL(SYNC)
   CALL PostStreamMIC
!DIR$ END OFFLOAD
   CALL PostStream
!DIR$ OFFLOAD_WAIT TARGET(mic:0) WAIT(SYNC)

   now = NG - now
   nxt = NG - nxt

   now_mic = NG_MIC - now_mic
   nxt_mic = NG_MIC - nxt_mic

#ifdef LBMIO
   IF( MOD(iStep,tStat) == 0 ) CALL Stats
   IF( MOD(iStep,tSave) == 0 ) CALL VtkSave
#endif

 END DO
 tCurrent = omp_get_wtime()
! nIO is the approximate number of times we save the full domain data to file
! In most sane implementations this dominates the IO cost.
 nIO   = (CompletedRelaxSteps + MaxStep) / tSave
 tEvol = tCurrent - tIO*nIO - tRelax - tSetup - tStart
 tRun  = tCurrent - tStart

!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) OUT( f0_mic, f1_mic, f2_mic : REUSE FREE )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) OUT( f3_mic,f4_mic,f5_mic, f6_mic : REUSE FREE )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) OUT( g0_mic, g1_mic, g2_mic, g3_mic : REUSE FREE )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) OUT( g4_mic, g5_mic, g6_mic : REUSE FREE )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) OUT( g7_mic, g8_mic, g9_mic, g10_mic : REUSE FREE )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) OUT( g11_mic, g12_mic, g13_mic, g14_mic : REUSE FREE )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) OUT( g15_mic, g16_mic, g17_mic, g18_mic : REUSE FREE )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) OUT( phi_mic, lapPhi_mic, gradPhiX_mic : REUSE FREE )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) OUT( gradPhiY_mic,gradPhiZ_mic,bubbles : REUSE FREE )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) OUT( f_buff_mic, f_buff_cpu : REUSE FREE )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) OUT( buff_mic, buff_cpu : REUSE FREE )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) OUT( buff_mic_edge_x, buff_cpu_edge_x : REUSE FREE )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) OUT( buff_mic_edge_z, buff_cpu_edge_z : REUSE FREE )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) OUT( buff_mic_phi, buff_cpu_phi : REUSE FREE )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) OUT( buff_mic_phi_edge_x, buff_cpu_phi_edge_x : REUSE FREE )
!DIR$ OFFLOAD_TRANSFER TARGET(MIC:0) OUT( buff_mic_phi_edge_z, buff_cpu_phi_edge_z : REUSE FREE )

! Save final data
 CALL FinalSave
#ifdef LBMIO
 CALL VtkSave
#endif

! Free memory, finalize MPI and end program
 CALL MemAlloc(2)

 END PROGRAM
