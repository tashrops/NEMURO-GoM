C $Header: /u/gcmpack/MITgcm/pkg/ptracers/ptracers_calc_wsurf_tr.F,v 1.4 2016/08/23 18:29:49 jahn Exp $
C $Name:  $

#include "PTRACERS_OPTIONS.h"

CBOP
C     !ROUTINE: PTRACERS_CALC_WSURF_TR
C     !INTERFACE:
      SUBROUTINE PTRACERS_CALC_WSURF_TR(wVelFld,
     I                                  myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE CALC_WSURF_TR
C     | o Compute a correction for the source/sink of tracer
C     |   due to the linear free surface.
C     | o The source/sink results from W*Tr not summing to
C     |   zero at the free surface.
C     | o Here, we compute an area-integrated correction
C     |   to be applied in biogeochemistry packages.
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global variables
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SURFACE.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"

C     !INPUT PARAMETERS:
C     wvelFld  :: vertical velocity field
C     myTime   :: Current time in simulation
C     myIter   :: Current iteration number in simulation
C     myThid   :: Thread number for this instance of the routine.
      _RL wVelFld (1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RL myTime
      INTEGER myIter
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
      INTEGER i,j,bi,bj,ks,itr
      _RL wT_Tile(nSx,nSy,PTRACERS_num)
      _RL wT_Mean(PTRACERS_num)

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      DO itr=1,PTRACERS_numInUse
       totSurfCorPTr(iTr)  = 0.0 _d 0
       meanSurfCorPTr(iTr) = 0.0 _d 0
       wT_mean(iTr)        = 0.0 _d 0
      ENDDO

C-- compute per-tile sums
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO itr=1,PTRACERS_numInUse
         wT_Tile(bi,bj,itr) = 0.
         DO j=1,sNy
          DO i=1,sNx
           ks = ksurfC(i,j,bi,bj)
           IF (ks.LE.Nr) THEN
            wT_Tile(bi,bj,itr) = wT_Tile(bi,bj,itr)
     &       + rA(i,j,bi,bj)*wVelFld(i,j,ks,bi,bj)
     &                      *pTracer(i,j,ks,bi,bj,iTr)
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C NOTE: 03/19/19 found that line 70 had "ptracer" instead of pTracer?
C Dont know if this had any impact.

C-- sum over tiles
      DO itr=1,PTRACERS_numInUse
       CALL GLOBAL_SUM_TILE_RL( wT_Tile(1,1,itr), wT_Mean(iTr),
     &                          myThid )
       totSurfCorPTr(iTr) = wT_Mean(iTr)
      ENDDO
        
C      IF ( globalArea.GT.0. ) THEN
C       _BEGIN_MASTER( myThid )
C       DO itr=1,PTRACERS_numInUse
C        meanSurfCorPTr(iTr) = wT_Mean(iTr) / globalArea
C        ENDDO
        
C       _END_MASTER( myThid )
C      ENDIF

      _BARRIER

      RETURN
      END
