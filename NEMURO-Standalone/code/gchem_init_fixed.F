C $Header: /u/gcmpack/MITgcm/pkg/gchem/gchem_init_fixed.F,v 1.12 2016/05/23 13:09:31 jmc Exp $
C $Name:  $

#include "GCHEM_OPTIONS.h"

C !INTERFACE: ==========================================================
      SUBROUTINE GCHEM_INIT_FIXED(myThid )

C !DESCRIPTION:
C calls subroutines that initialized fixed variables for any tracer
C experiment

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GCHEM.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid
      INTEGER        diagNum
      INTEGER        diagMate
      CHARACTER*8    diagName
      CHARACTER*16   diagCode
      CHARACTER*16   diagUnits
      CHARACTER*(80) diagTitle
CEOP

#ifdef ALLOW_GCHEM

C     !LOCAL VARIABLES:

C jmc: simpler to put MASTER thread directive within individual S/R :
C      that leave the choice, for each S/R, to do either single (master)
C      thread tasks or Multi-threaded tasks.
c        _BEGIN_MASTER( myThid )
ccccccccccccccccccccccccccc
c initialize parameters   c
ccccccccccccccccccccccccccc
#ifdef ALLOW_CFC
      IF ( useCFC ) THEN
         print*,'QQ start cfc params'
         CALL CFC_PARAM( myThid )
         print*,'QQ read cfc atmos data'
         CALL CFC_ATMOS( myThid )
      ENDIF
#endif /* ALLOW_CFC */
#ifdef ALLOW_DIC
      IF ( useDIC ) THEN
         print*,'QQ load dic parameters, initial fixed'
         CALL DIC_INIT_FIXED( myThid )
      ENDIF
#endif /* ALLOW_DIC */
#ifdef ALLOW_BLING
      IF ( useBLING ) THEN
         CALL BLING_INIT_FIXED( myThid )
      ENDIF
#endif /* ALLOW_BLING */
c         _END_MASTER( myThid )

#ifdef ALLOW_DARWIN
      IF ( useDARWIN ) THEN
         CALL DARWIN_INIT_FIXED(myThid)
      ENDIF
#endif

#ifdef ALLOW_DIAGNOSTICS
C     Define diagnostics Names :
      IF ( useDiagnostics ) THEN
        CALL GCHEM_DIAGNOSTICS_INIT( myThid )

C     ***  Edited by Taylor Shropshire ***

      diagName  = 'phy_chl     '
      diagTitle = 'Chlorophyll-a Concentration'
      diagUnits = 'mg Chl-a/m3'
      diagCode  = 'SMR     MR     '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'l_lim_sp_field     '
      diagTitle = 'Small Phytoplankton Light Limitation'
      diagUnits = 'value 0-1'
      diagCode  = 'SMR     MR     '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'l_lim_lp_field     '
      diagTitle = 'Large Phytoplankton Light Limitation'
      diagUnits = 'value 0-1'
      diagCode  = 'SMR     MR     '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'tmp_field     '
      diagTitle = 'Water Temperature'
      diagUnits = 'deg C'
      diagCode  = 'SMR     MR     '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'surf_chl     '
      diagTitle = 'Chlorophyll-a Concentration'
      diagUnits = 'mg Chl-a/m3'
      diagCode  = 'SM      M1      '      
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'sp2sz_ifield     '
      diagTitle = 'SZ Grazing on SP'
      diagUnits = 'mmol N/m2 per day'
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'lp2lz_ifield     '
      diagTitle = 'lZ Grazing on LP'
      diagUnits = 'mmol N/m2 per day'
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'sz2lz_ifield     '
      diagTitle = 'LZ Grazing on SZ'
      diagUnits = 'mmol N/m2 per day'
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'lp2pz_ifield     '
      diagTitle = 'PZ Grazing on LP'
      diagUnits = 'mmol N/m2 per day'
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'sz2pz_ifield     '
      diagTitle = 'PZ Grazing on SZ'
      diagUnits = 'mmol N/m2 per day'
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'lz2pz_ifield     '
      diagTitle = 'PZ Grazing on LZ'
      diagUnits = 'mmol N/m2 per day'
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'lz_diet_field     '
      diagTitle = 'LZ diet'
      diagUnits = '1-LP, 2-SZ, 3-LZ'
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'pz_diet_field     '
      diagTitle = 'PZ diet'
      diagUnits = '1-LP, 2-SZ, 3-LZ'
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'lz_trophic_field     '
      diagTitle = 'LZ trophic level'
      diagUnits = 'trophic level'
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'pz_trophic_field     '
      diagTitle = 'PZ trophic level'
      diagUnits = 'tropohic level'
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'sp_npp_ifield     '
      diagTitle = 'SP NPP'
      diagUnits = 'mmol N/m2 per day'
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'lp_npp_ifield     '
      diagTitle = 'LP NPP'
      diagUnits = 'mmol N/m2 per day'
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'sp_ifield     '
      diagTitle = 'SP Biomass'
      diagUnits = 'mmol N/m2'
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'lp_ifield     '
      diagTitle = 'LP Biomass'
      diagUnits = 'mmol N/m2'
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'sz_ifield     '
      diagTitle = 'SZ Biomass'
      diagUnits = 'mmol N/m2'
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'lz_ifield     '
      diagTitle = 'LZ Biomass'
      diagUnits = 'mmol N/m2'
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'pz_ifield     '
      diagTitle = 'PZ Biomass'
      diagUnits = 'mmol N/m2'
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )


C       End New Diagnostics
C       ******************************

      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_GCHEM */

      RETURN
      END
