C $Header: /u/gcmpack/MITgcm/pkg/gchem/gchem_readparms.F,v 1.20 2014/05/27 23:41:28 jmc Exp $
C $Name:  $

#include "GCHEM_OPTIONS.h"

CBOP
C !ROUTINE: GCHEM_READPARMS

C !INTERFACE: ==========================================================
      SUBROUTINE GCHEM_READPARMS( myThid )

C !DESCRIPTION:
C     Initialize GCHEM parameters, read in data.gchem

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GCHEM.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_GCHEM

C !LOCAL VARIABLES: ====================================================
C  tIter0               :: retired parameter
C  iTracer              :: loop indices
C  iUnit                :: unit number for I/O
C  msgBuf               :: message buffer
      INTEGER tIter0
      INTEGER iUnit
      CHARACTER*(MAX_LEN_MBUF) msgBuf
C      integer:: dd
C      real:: temp
C      real,dimension(365):: light_forcing
C      character(len=46):: part1
C      character(len=29):: part2
C      character(len=200):: fname_l
CEOP

C- Sub-package on/off flags: not fully implemented, requires
C  to test the flag before any corresponding pkg S/R call

      NAMELIST /GCHEM_PARM01/
     &                   nsubtime,     
     &                   useDIC,
     &                   useCFC,
     &                   useDARWIN,
     &                   Filename1,
     &                   Filename2,
     &                   Filename3,
     &                   Filename4,
     &                   Filename5,
     &           day_count, bio_count, year_count,
     &           day_count_tot, bad_cell_count, bio_count2,
     &           dumpFreq2,
C
C Small Phytoplankton Variables:
     &           vmax_sp, k_NO_sp, k_NH_sp, inh_NH_NO_sp,
     &           ref_resp_sp, ref_mort_sp, ext_excr_sp,
     &           tc_v_sp, tc_m_sp, tc_r_sp, photo_chem_sp,
     &           photo_inh_sp, ext_sp,
C Large Phytoplankton 
     &           vmax_lp, k_NO_lp, k_NH_lp, inh_NH_NO_lp,
     &           ref_resp_lp, ref_mort_lp, ext_excr_lp,
     &           tc_v_lp, tc_m_lp, tc_r_lp, photo_chem_lp,
     &           photo_inh_lp, ext_lp, k_SI_lp,
C
C Small Zooplankton:
     &           gmax_sz_sp, tc_g_sz_sp, iv_sz_sp, thresh_sz_sp,
     &           ref_mort_sz, tc_m_sz, ae_sz, gge_sz,
C
C Large Zooplankton:
     &           gmax_lz_sp, tc_g_lz_sp, iv_lz_sp, thresh_lz_sp,
     &           gmax_lz_lp, tc_g_lz_lp, iv_lz_lp, thresh_lz_lp,
     &           gmax_lz_sz, tc_g_lz_sz, iv_lz_sz, thresh_lz_sz,
     &           ref_mort_lz, tc_m_lz, ae_lz, gge_lz,
C
C Predatory Zooplankton:
     &           gmax_pz_lp, tc_g_pz_lp, iv_pz_lp, thresh_pz_lp,
     &           gmax_pz_sz, tc_g_pz_sz, iv_pz_sz, thresh_pz_sz,
     &           gmax_pz_lz, tc_g_pz_lz, iv_pz_lz, thresh_pz_lz,
     &           ref_mort_pz, tc_m_pz, ae_pz, gge_pz, inh_szlz_lp,
     &           inh_lz_sz,
C
C Nitrogen:
     &           ref_nitr, tc_nitr, ref_dec_PON_NH, tc_dec_PON_NH,
     &           ref_dec_PON_DON, tc_dec_PON_DON, ref_dec_DON_NH,
     &           tc_dec_DON_NH,
C
C Silica:
     &           ref_dec_OP_SI, tc_dec_OP_SI, r_SI_N, r_SI_N_riv,
C
C Other:
     &           min_val, dt_b, dt_b_sec, PAR_frac, ext_w, 
     &           scale_tmp_phyto, scale_tmp_phyto2, scale_tmp_zoo,
     &           scale_tmp_decomp, diag_tscale, sed_min, sed_max,
     &           sink_min, sink_max, PON_k, PON_flux_k,
     &           sig_slope1, sig_mid1, sig_slope2, sig_mid2,

C
C *** End NEMURO Variables Names ***
C       
C Chl:C Model
     &           sp_chl2c_min, sp_chl2c_max,
     &           lp_chl2c_min, lp_chl2c_max,
     &           alpha_chl,    
C
C Diagnostic Fields 
C     &          sp_npp, lp_npp, phy_chl, np_frac,
C     &          sz_graze, lz_graze, pz_graze, 
C
C Other:
     &           gchem_ForcingPeriod, gchem_ForcingCycle,
     &           tIter0

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      IF ( .NOT.useGCHEM ) THEN
C-    pkg GCHEM is not used
        _BEGIN_MASTER(myThid)
C-    Track pkg activation status:
C     print a (weak) warning if data.gchem is found
         CALL PACKAGES_UNUSED_MSG( 'useGCHEM', ' ', ' ' )
        _END_MASTER(myThid)
        RETURN
      ENDIF

      _BEGIN_MASTER(myThid)

C Set defaults values for parameters in GCHEM.h
       useDIC = .FALSE.
       useCFC = .FALSE.
       useDARWIN = .FALSE.
       tIter0 = UNSET_I
       nsubtime = 1
       day_count = 0
       bio_count = 0
       bio_count2 = 0
       dumpFreq2= 0
       year_count = 1
       day_count_tot = 0
       bad_cell_count = 0 
       Filename1=' '
       Filename2=' '
       Filename3=' '
       Filename4=' '
       Filename5=' '
C
C Small Phytoplankton Parameters
       vmax_sp = 0.0
       k_NO_sp = 0.0                 
       k_NH_sp = 0.0
       inh_NH_NO_sp = 0.0
       ref_resp_sp = 0.0
       ref_mort_sp = 0.0
       ext_excr_sp = 0.0
       tc_v_sp = 0.0
       tc_m_sp = 0.0
       tc_r_sp = 0.0
       photo_chem_sp = 0.0
       photo_inh_sp = 0.0
       ext_sp = 0.0
C
C Large Phytoplankton Parameters       
       vmax_lp = 0.0
       k_NO_lp = 0.0
       k_NH_lp = 0.0
       k_SI_lp = 0.0
       inh_NH_NO_lp = 0.0
       ref_resp_lp = 0.0
       ref_mort_lp = 0.0
       ext_excr_lp = 0.0
       tc_v_lp = 0.0
       tc_m_lp = 0.0
       tc_r_lp = 0.0
       photo_chem_lp = 0.0
       photo_inh_lp = 0.0
       ext_lp = 0.0
C
C Small Zooplankton 
       gmax_sz_sp = 0.0
       tc_g_sz_sp = 0.0
       iv_sz_sp = 0.0
       thresh_sz_sp = 0.0
       ref_mort_sz = 0.0
       tc_m_sz = 0.0
       ae_sz = 0.0
       gge_sz = 0.0
C
C Large Zooplankton 
       gmax_lz_sp = 0.0
       tc_g_lz_sp = 0.0
       iv_lz_sp = 0.0
       thresh_lz_sp = 0.0
       gmax_lz_lp = 0.0
       tc_g_lz_lp = 0.0
       iv_lz_lp = 0.0
       thresh_lz_lp = 0.0
       gmax_lz_sz = 0.0
       tc_g_lz_sz = 0.0
       iv_lz_sz = 0.0
       thresh_lz_sz = 0.0
       ref_mort_lz = 0.0
       tc_m_lz = 0.0
       ae_lz = 0.0
       gge_lz = 0.0
C
C Predatory Zooplankton 
       gmax_pz_lp = 0.0 
       tc_g_pz_lp = 0.0
       iv_pz_lp = 0.0
       thresh_pz_lp = 0.0
       gmax_pz_sz = 0.0
       tc_g_pz_sz = 0.0
       iv_pz_sz = 0.0
       thresh_pz_sz = 0.0
       gmax_pz_lz = 0.0
       tc_g_pz_lz = 0.0
       iv_pz_lz = 0.0
       thresh_pz_lz = 0.0
       ref_mort_pz = 0.0
       tc_m_pz = 0.0
       ae_pz = 0.0
       gge_pz = 0.0
       inh_szlz_lp = 0.0
       inh_lz_sz = 0.0
C
C Nitrogen: 
       ref_nitr = 0.0
       tc_nitr = 0.0 
       ref_dec_PON_NH = 0.0 
       tc_dec_PON_NH = 0.0
       ref_dec_PON_DON = 0.0
       tc_dec_PON_DON = 0.0
       ref_dec_DON_NH = 0.0
       tc_dec_DON_NH = 0.0
C 
C Silica: 
       ref_dec_OP_SI = 0.0
       tc_dec_OP_SI = 0.0
       r_SI_N = 0.0
       r_SI_N_riv= 0.0
C
C Other:
        min_val = 0.0 
        dt_b = 0.0
        PAR_frac = 0.0
        ext_w = 0.0
        dt_b_sec = 0.0
        scale_tmp_phyto = 0.0
        scale_tmp_phyto2 = 0.0
        scale_tmp_zoo = 0.0
        scale_tmp_decomp = 0.0
        diag_tscale = 0.0
        sed_min = 0.0
        sed_max = 0.0
        sink_min = 0.0
        sink_max = 0.0
        PON_k = 0.0
        PON_flux_k = 0.0
        sig_slope1 = 0.0
        sig_mid1 = 0.0
        sig_slope2 = 0.0
        sig_mid2 = 0.0

C       
C Chl:C Model 
       sp_chl2c_min=0.0
       sp_chl2c_max=0.0
       lp_chl2c_min=0.0
       lp_chl2c_max=0.0
       alpha_chl=0.0      
        
C
C default periodic forcing to same as for physics
       gchem_ForcingPeriod=externForcingPeriod
       gchem_ForcingCycle=externForcingCycle

C Open and read the data.gchem file
      WRITE(msgBuf,'(A)') ' GCHEM_READPARMS: opening data.gchem'
      CALL PRINT_MESSAGE(msgBuf, standardMessageUnit,
     &                   SQUEEZE_RIGHT , 1)
      CALL OPEN_COPY_DATA_FILE(
     I                   'data.gchem', 'GCHEM_PARM01',
     O                   iUnit,
     I                   myThid )
      READ(UNIT=iUnit,NML=GCHEM_PARM01)
      WRITE(msgBuf,'(A)')
     &  ' GCHEM_READPARMS: finished reading data.gchem'
      CALL PRINT_MESSAGE(msgBuf, standardMessageUnit,
     &                   SQUEEZE_RIGHT , 1)

C Close the open data file
      CLOSE(iUnit)

C- Check for retired parameters:
      IF ( tIter0 .NE. UNSET_I ) THEN
c      nRetired = nRetired+1
       WRITE(msgBuf,'(A,A)')
     &  'S/R GCHEM_READPARMS: Paramater "tIter0" is',
     &  ' no longer allowed in file "data.gchem"'
       CALL PRINT_ERROR( msgBuf , myThid)
       WRITE(msgBuf,'(A,A)')
     &  'S/R GCHEM_READPARMS: "tIter0" has been moved to',
     &  ' PTRACERS_Iter0 in file "data.ptracers".'
       CALL PRINT_ERROR( msgBuf , myThid)
       STOP 'ABNORMAL END: S/R GCHEM_READPARMS'
      ENDIF

      _END_MASTER(myThid)

C Everyone else must wait for the parameters to be loaded
      _BARRIER

#ifdef ALLOW_DIC
      IF ( useDIC ) THEN
        CALL DIC_READPARMS(myThid)
      ENDIF
#endif

#ifdef ALLOW_CFC
      IF ( useCFC ) THEN
        CALL CFC_READPARMS(myThid)
      ENDIF
#endif

#ifdef ALLOW_DARWIN
      IF ( useDARWIN ) THEN
        CALL DARWIN_READPARMS(myThid)
      ENDIF
#endif

#endif /* ALLOW_GCHEM */

      RETURN
      END
