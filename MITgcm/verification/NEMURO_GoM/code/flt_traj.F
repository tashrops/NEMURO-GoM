C $Header: /u/gcmpack/MITgcm/pkg/flt/flt_traj.F,v 1.14 2012/03/30 18:25:03 jmc Exp $
C $Name:  $

#include "FLT_OPTIONS.h"

CBOP 0
C !ROUTINE: FLT_TRAJ

C !INTERFACE:
      SUBROUTINE FLT_TRAJ (
     I                      myTime, myIter, myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE FLT_TRAJ
C     | o This routine samples the model state at float position
C     |   every flt_int_traj time steps and writes output.
C     *==========================================================*

C     !USES:
      IMPLICIT NONE
C     == global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FLT_SIZE.h"
#include "FLT.h"
#include "FLT_BUFF.h"
#ifdef ALLOW_EXCH2
#include "W2_EXCH2_SIZE.h"
#include "W2_EXCH2_TOPOLOGY.h"
#endif
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
C#include "GCHEM.h"
C#include "GCHEM_FIELDS.h"

C     !INPUT PARAMETERS:
C     myTime :: current time in simulation
C     myIter :: current iteration number
C     myThid :: my Thread Id number
      _RL myTime
      INTEGER myIter, myThid

C     !FUNCTIONS:
      _RL FLT_MAP_K2R
      EXTERNAL FLT_MAP_K2R

C     !LOCAL VARIABLES:
      INTEGER bi, bj, nFlds
      INTEGER ip, kp, ii
      _RL write_flag(max_npart_tile)
      _RL ix, jy, i0x, j0y, xx, yy, zz
      _RL uu, vv, tt, ss, pp
      _RL NO, NH, SI, DON, PON, OP, sp, lp, sz, lz, pz
      _RL chl, l_lim_sp, l_lim_lp, temp, cm, sp_chl2c, lp_chl2c

      INTEGER imax
      PARAMETER (imax=30)
      _RL tmp(imax)
      _RL npart_read, npart_times
      _RS dummyRS(1)
      INTEGER fp, ioUnit, irecord, npart_count
      CHARACTER*(MAX_LEN_FNAM) fn
      CHARACTER*(MAX_LEN_MBUF) msgBuf
#ifdef ALLOW_EXCH2
      INTEGER nT
#endif
CEOP

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C--   set number of fields to write
      nFlds = 0
      IF ( flt_selectTrajOutp.GE.1 ) nFlds = nFlds + 8
      IF ( flt_selectTrajOutp.GE.2 ) nFlds = nFlds + 22

C--   check buffer size
      IF ( nFlds.GT.fltBufDim ) THEN
         _BEGIN_MASTER(myThid)
         WRITE(msgBuf,'(3(A,I4))') ' FLT_TRAJ: fltBufDim=', fltBufDim,
     &                             ' too small (<', nFlds, ' )'
         CALL PRINT_ERROR( msgBuf, myThid )
         WRITE(msgBuf,'(2A)')     ' FLT_TRAJ: => increase fltBufDim',
     &                            ' in "FLT_SIZE.h" & recompile'
         CALL PRINT_ERROR( msgBuf, myThid )
         _END_MASTER(myThid)
         CALL ALL_PROC_DIE( myThid )
         STOP 'ABNORMAL END: S/R FLT_TRAJ'
      ENDIF

      IF ( myIter.EQ.nIter0 .OR. flt_selectTrajOutp.LE.0 ) RETURN

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Calculate position + other fields at float position and fill up IO-buffer

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

#ifdef ALLOW_EXCH2
         nT = W2_myTileList(bi,bj)
         i0x = DFLOAT( exch2_txGlobalo(nT) - 1 )
         j0y = DFLOAT( exch2_tyGlobalo(nT) - 1 )
#else
         i0x = DFLOAT( myXGlobalLo-1 + (bi-1)*sNx )
         j0y = DFLOAT( myYGlobalLo-1 + (bj-1)*sNy )
#endif
         DO ip=1,npart_tile(bi,bj)

            ix = ipart(ip,bi,bj)
            jy = jpart(ip,bi,bj)
            CALL FLT_MAP_IJLOCAL2XY( xx, yy,
     I                               ix, jy, bi,bj, myThid )
            zz = FLT_MAP_K2R( kpart(ip,bi,bj),bi,bj,myThid )
            kp = NINT(kpart(ip,bi,bj))
            tmp(1) = npart(ip,bi,bj)
            tmp(2) = myTime
            tmp(3) = xx
            tmp(4) = yy
            tmp(5) = zz
            tmp(6) = ix + i0x
            tmp(7) = jy + j0y
            tmp(8) = kpart(ip,bi,bj)

            IF ( ( flt_selectTrajOutp.GE.2 )   .AND.
     &           ( myTime.GE.tstart(ip,bi,bj)) .AND.
     &           ( tend(ip,bi,bj).EQ.-1. .OR. myTime.LE.tend(ip,bi,bj))
     &         ) THEN
              IF ( kp.LT.1 .OR. kp.GT.Nr ) THEN
                WRITE(msgBuf,'(2A,I8)') '** WARNING ** FLT_TRAJ: ',
     &            ' illegal value for kp=',kp
                CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                              SQUEEZE_RIGHT, myThid )
                WRITE(msgBuf,'(A,1P5E20.13)')
     &            ' FLT_TRAJ: ', (flt_io_buff(ii,ip,bi,bj),ii=1,5)
                CALL PRINT_MESSAGE( msgBuf, errorMessageUnit,
     &                              SQUEEZE_RIGHT, myThid )
c               CALL PRINT_ERROR( msgBuf, myThid )
c               STOP 'ABNORMAL END: S/R FLT_TRAJ'
C-- jmc: not sure if this is right but added to avoid Pb in FLT_BILINEAR:
                kp = MIN( MAX(kp,1), Nr)
              ENDIF
              CALL FLT_BILINEAR  (ix,jy,uu,uVel,  kp,1,bi,bj,myThid)
              CALL FLT_BILINEAR  (ix,jy,vv,vVel,  kp,2,bi,bj,myThid)
              CALL FLT_BILINEAR2D(ix,jy,pp,depth,     0,bi,bj,myThid)
              CALL FLT_BILINEAR  (ix,jy,ss,salt,  kp,0,bi,bj,myThid) 

              tmp( 9) = pp
              tmp(10) = uu
              tmp(11) = vv
              tmp(12) = ss

        CALL FLT_BILINEAR  (ix,jy,NO,pTracer_tavg_field2(:,:,:,:,:,1),
     &  kp,0,bi,bj,myThid)
        CALL FLT_BILINEAR  (ix,jy,NH,pTracer_tavg_field2(:,:,:,:,:,2),
     &  kp,0,bi,bj,myThid)
        CALL FLT_BILINEAR  (ix,jy,SI,pTracer_tavg_field2(:,:,:,:,:,3),
     &  kp,0,bi,bj,myThid)
        CALL FLT_BILINEAR  (ix,jy,DON,pTracer_tavg_field2(:,:,:,:,:,4),
     &  kp,0,bi,bj,myThid)
        CALL FLT_BILINEAR  (ix,jy,PON,pTracer_tavg_field2(:,:,:,:,:,5),
     &  kp,0,bi,bj,myThid)
        CALL FLT_BILINEAR  (ix,jy,OP,pTracer_tavg_field2(:,:,:,:,:,6),
     &  kp,0,bi,bj,myThid)
        CALL FLT_BILINEAR  (ix,jy,sp,pTracer_tavg_field2(:,:,:,:,:,7),
     &  kp,0,bi,bj,myThid)
        CALL FLT_BILINEAR  (ix,jy,lp,pTracer_tavg_field2(:,:,:,:,:,8),
     &  kp,0,bi,bj,myThid)
        CALL FLT_BILINEAR  (ix,jy,sz,pTracer_tavg_field2(:,:,:,:,:,9),
     &  kp,0,bi,bj,myThid)
        CALL FLT_BILINEAR  (ix,jy,lz,pTracer_tavg_field2(:,:,:,:,:,10),
     &  kp,0,bi,bj,myThid)
        CALL FLT_BILINEAR  (ix,jy,pz,pTracer_tavg_field2(:,:,:,:,:,11),
     &  kp,0,bi,bj,myThid)

              tmp(13) = NO
              tmp(14) = NH
              tmp(15) = SI
              tmp(16) = DON
              tmp(17) = PON
              tmp(18) = OP
              tmp(19) = sp
              tmp(20) = lp
              tmp(21) = sz
              tmp(22) = lz
              tmp(23) = pz

              CALL FLT_BILINEAR  (ix,jy,chl,phy_chl2,
     &  kp,0,bi,bj,myThid)
              CALL FLT_BILINEAR  (ix,jy,l_lim_sp,l_lim_sp_field2,
     &  kp,0,bi,bj,myThid)
              CALL FLT_BILINEAR  (ix,jy,l_lim_lp,l_lim_lp_field2,
     &  kp,0,bi,bj,myThid)
              CALL FLT_BILINEAR  (ix,jy,temp,tmp_field2,
     &  kp,0,bi,bj,myThid)
              CALL FLT_BILINEAR  (ix,jy,cm,hFacC,
     &  kp,0,bi,bj,myThid)
              CALL FLT_BILINEAR  (ix,jy,sp_chl2c,chl2c_sp_field2,
     &  kp,0,bi,bj,myThid)
              CALL FLT_BILINEAR  (ix,jy,lp_chl2c,chl2c_lp_field2,
     &  kp,0,bi,bj,myThid)
              
              tmp(24) = chl
              tmp(25) = l_lim_sp
              tmp(26) = l_lim_lp
              tmp(27) = temp
              tmp(28) = sp_chl2c
              tmp(29) = lp_chl2c
              tmp(30) = cm

                write_flag(ip)=1.0

            ELSEIF ( flt_selectTrajOutp.GE.2 ) THEN
              tmp( 9) = flt_nan
              tmp(10) = flt_nan
              tmp(11) = flt_nan
              tmp(12) = flt_nan
              tmp(13) = flt_nan
              tmp(14) = flt_nan
              tmp(15) = flt_nan
              tmp(16) = flt_nan
              tmp(17) = flt_nan
              tmp(18) = flt_nan
              tmp(19) = flt_nan
              tmp(20) = flt_nan
              tmp(21) = flt_nan
              tmp(22) = flt_nan
              tmp(23) = flt_nan
              tmp(24) = flt_nan
              tmp(25) = flt_nan
              tmp(26) = flt_nan
              tmp(27) = flt_nan
              tmp(28) = flt_nan
              tmp(29) = flt_nan
              tmp(30) = flt_nan
              
              write_flag(ip)=0.0

            ENDIF

            DO ii=1,nFlds
              flt_io_buff(ii,ip,bi,bj) = tmp(ii)
            ENDDO

         ENDDO

       ENDDO
      ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C--   Write shared buffer to file
        
      _BARRIER
      _BEGIN_MASTER(myThid)

      fn = 'float_trajectories'
      fp = writeBinaryPrec

      DO bj=1,nSy
       DO bi=1,nSx

CC (1) read actual number floats from file (if exists)
C         ioUnit = -2
C         CALL MDS_READVEC_LOC(  fn, fp, ioUnit, 'RL', nFlds,
C     O                          tmp, dummyRS,
C     I                          bi, bj, 1, myThid )
C         IF ( ioUnit.GT. 0 ) THEN
C            npart_read  = tmp(1)
C            npart_times = tmp(5)
C            ii = NINT(tmp(7))
CC-       for backward compatibility with old trajectory files:
C            IF ( ii.EQ.0 ) ii = 13
C            IF ( ii.NE.nFlds ) THEN
C              WRITE(msgBuf,'(A,I4,A)')
C     &            'FLT_TRAJ: nFlds=', nFlds,' different from'
C              CALL PRINT_ERROR( msgBuf, myThid )
C              WRITE(msgBuf,'(3A,I4,A)')
C     &            'previous file (',fn(1:18),') value =',ii
C              CALL PRINT_ERROR( msgBuf, myThid )
C              CALL ALL_PROC_DIE( 0 )
C              STOP 'ABNORMAL END: S/R FLT_TRAJ'
C            ENDIF
CC-       close the read-unit (safer to use a different unit for writing)
C            CLOSE( ioUnit )
C         ELSE
C            npart_read  = 0.
C            npart_times = 0.
C            tmp(2)      = myTime
C         ENDIF

C (2) write new actual number floats and time axis into file
C-    the standard routine mds_writevec_loc can be used here

CC     total number of records in this file
CC         tmp(1) = DBLE(npart_tile(bi,bj))+npart_read
C         tmp(1) = npart_read+DBLE(npart_tile(bi,bj))+npart_read
CC     first time of writing floats (do not change when written)
Cc        tmp(2) = tmp(2)
CC     current time
C         tmp(3) = myTime
CC     timestep
C         tmp(4) = flt_int_traj
CC     total number of timesteps
C         tmp(5) = npart_times + 1.
CC     total number of floats
C         tmp(6) = max_npart
CC     total number of fields
C         tmp(7) = nFlds
C         DO ii=8,nFlds
C           tmp(ii) = 0.
C         ENDDO
C        ioUnit = -1
C         CALL MDS_WRITEVEC_LOC( fn, fp, ioUnit, 'RL', nFlds,
C     &                          tmp, dummyRS,
C     &                          bi, bj, -1, myIter, myThid )

C ********** Write FLT File ********** C
        ioUnit = -1
        DO ip=1,npart_tile(bi,bj)
        if (write_flag(ip)==1.0) then 
C        if (flt_io_buff(9,ip,bi,bj)/=-999) then 

        part_rec = part_rec+1
        irecord = part_rec*-1

            CALL MDS_WRITEVEC_LOC( fn, fp, ioUnit, 'RL', nFlds,
     I                             flt_io_buff(1,ip,bi,bj), dummyRS,
     I                             bi, bj, irecord, myIter, myThid )

        end if 
         ENDDO
         CLOSE( ioUnit )
C ************************************** C

       ENDDO
      ENDDO

      _END_MASTER(myThid)
      _BARRIER

      RETURN
      END
