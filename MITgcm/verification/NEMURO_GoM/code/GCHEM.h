C $Header: /u/gcmpack/MITgcm/pkg/gchem/GCHEM.h,v 1.13 2009/06/30 16:43:00 jahn Exp $
C $Name:  $

#ifdef ALLOW_GCHEM

CBOP
C    !ROUTINE: GCHEM.h
C    !INTERFACE:

C    !DESCRIPTION:
C Contains tracer parameters and input files for chemical tracers.
C These can be read in from data.gchem
C
C--   COMMON /GCHEM_PARM_L/ Logical valued parameters used by GCHEM pkg.
C     useDIC    :: flag to turn on/off DIC pkg
C     useCFC    :: flag to turn on/off CFC pkg
C     useDARWIN :: flag to turn on/off darwin pkg
C
C--   COMMON /GCHEM_PARAMS/
C  nsubtime    :: number of chemistry timesteps per deltaTtracer
C                 (default 1)
C  Filename*   :: various spare filenames
C  gchem_int*  :: place holder to read in a integer number, set at run time
C  gchem_rl*   :: place holder to read in a real number, set at run time
C  gchem_ForcingPeriod :: periodic forcing parameter specific for gchem (seconds)
C  gchem_ForcingCycle  :: periodic forcing parameter specific for gchem (seconds)

CEOP

      COMMON /GCHEM_PARM_L/
     &              useDIC,
     &              useCFC,
     &              useDARWIN

      LOGICAL useDIC, useCFC, useDARWIN

      COMMON /GCHEM_PARAMS/
     &                   Filename1,
     &                   Filename2,
     &                   Filename3,
     &                   Filename4,
     &                   Filename5,
     &                   nsubtime, 
     &                   bio_count, day_count, year_count, 
     &                   day_count_tot, bad_cell_count, bio_count2,
     &			 dumpFreq2,
C
C *** Define Variable Names For NEMURO *** 
C
C Small Phytoplankton Variables:
     &           vmax_sp, k_NO_sp, k_NH_sp, inh_NH_NO_sp,
     & 		 ref_resp_sp, ref_mort_sp, ext_excr_sp,
     & 		 tc_v_sp, tc_m_sp, tc_r_sp, photo_chem_sp,
     & 		 photo_inh_sp, ext_sp, 
C
C Large Phtoplankton Variables: 
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
C Diagnostics 
C     &          sp_npp, lp_npp, np_frac,
C     & 	 sz_graze, lz_graze, pz_graze, 
C
C Chl:C	
     &		 sp_chl2c_min, sp_chl2c_max,
     & 	 	 lp_chl2c_min, lp_chl2c_max,
     &		 alpha_chl,
C
C Other:
     &	         min_val, dt_b, dt_b_sec, PAR_frac, ext_w, 
     &           scale_tmp_phyto, scale_tmp_phyto2, 
     &		 scale_tmp_zoo, scale_tmp_decomp, 
     &		 diag_tscale, sed_min, sed_max, sink_min, sink_max,
     &           PON_k, PON_flux_k, sig_slope1, sig_mid1,
     &           sig_slope2, sig_mid2,

C
C *** End NEMURO Variables Names ***
C
     &           gchem_ForcingPeriod, gchem_ForcingCycle

      INTEGER nsubtime
      CHARACTER*(MAX_LEN_FNAM) Filename1
      CHARACTER*(MAX_LEN_FNAM) Filename2
      CHARACTER*(MAX_LEN_FNAM) Filename3
      CHARACTER*(MAX_LEN_FNAM) Filename4
      CHARACTER*(MAX_LEN_FNAM) Filename5
      INTEGER day_count
      INTEGER day_count_tot
      INTEGER year_count
      INTEGER bio_count
      INTEGER bio_count2
      INTEGER dumpFreq2      
      INTEGER bad_cell_count

C *** Define NEMURO Variables Type ***
 
C Small Phytoplankton Variables (#=13):
      _RL     vmax_sp
      _RL     k_NO_sp
      _RL     k_NH_sp
      _RL     inh_NH_NO_sp
      _RL     ref_resp_sp
      _RL     ref_mort_sp
      _RL     ext_excr_sp
      _RL     tc_v_sp
      _RL     tc_m_sp
      _RL     tc_r_sp
      _RL     photo_chem_sp
      _RL     photo_inh_sp
      _RL     ext_sp

C Large Phtoplankton Variables (#=14):
      _RL     vmax_lp
      _RL     k_NO_lp
      _RL     k_NH_lp
      _RL     k_SI_lp
      _RL     inh_NH_NO_lp
      _RL     ref_resp_lp
      _RL     ref_mort_lp
      _RL     ext_excr_lp
      _RL     tc_v_lp
      _RL     tc_m_lp
      _RL     tc_r_lp
      _RL     photo_chem_lp
      _RL     photo_inh_lp
      _RL     ext_lp

C Small Zooplankton (#=8):
      _RL     gmax_sz_sp
      _RL     tc_g_sz_sp
      _RL     iv_sz_sp
      _RL     thresh_sz_sp
      _RL     ref_mort_sz
      _RL     tc_m_sz
      _RL     ae_sz
      _RL     gge_sz

C Large Zooplakton (#=16):
      _RL     gmax_lz_sp
      _RL     tc_g_lz_sp
      _RL     iv_lz_sp
      _RL     thresh_lz_sp
      _RL     gmax_lz_lp
      _RL     tc_g_lz_lp
      _RL     iv_lz_lp
      _RL     thresh_lz_lp
      _RL     gmax_lz_sz
      _RL     tc_g_lz_sz
      _RL     iv_lz_sz
      _RL     thresh_lz_sz
      _RL     ref_mort_lz
      _RL     tc_m_lz
      _RL     ae_lz
      _RL     gge_lz

C Predatory Zooplankton (#=18):
      _RL     gmax_pz_lp
      _RL     tc_g_pz_lp
      _RL     iv_pz_lp
      _RL     thresh_pz_lp
      _RL     gmax_pz_sz
      _RL     tc_g_pz_sz
      _RL     iv_pz_sz
      _RL     thresh_pz_sz
      _RL     gmax_pz_lz
      _RL     tc_g_pz_lz
      _RL     iv_pz_lz
      _RL     thresh_pz_lz
      _RL     ref_mort_pz
      _RL     tc_m_pz
      _RL     ae_pz
      _RL     gge_pz
      _RL     inh_szlz_lp
      _RL     inh_lz_sz

C Nitrogen (#=8): 
      _RL     ref_nitr
      _RL     tc_nitr
      _RL     ref_dec_PON_NH
      _RL     tc_dec_PON_NH
      _RL     ref_dec_PON_DON
      _RL     tc_dec_PON_DON
      _RL     ref_dec_DON_NH
      _RL     tc_dec_DON_NH

C Silica (#=3):
      _RL     ref_dec_OP_SI
      _RL     tc_dec_OP_SI
      _RL     r_SI_N
      _RL     r_SI_N_riv

C Other
      _RL     min_val
      _RL     dt_b
      _RL     PAR_frac
      _RL     ext_w 
      _RL     dt_b_sec
      _RL     scale_tmp_phyto
      _RL     scale_tmp_phyto2
      _RL     scale_tmp_zoo
      _RL     scale_tmp_decomp
      _RL     diag_tscale
      _RL     sed_min
      _RL     sed_max
      _RL     sink_min
      _RL     sink_max
      _RL     PON_k
      _RL     PON_flux_k
      _RL     sig_slope1
      _RL     sig_mid1
      _RL     sig_slope2
      _RL     sig_mid2

C *** END NEMURO Variables Type ***

C Diagnostics
C Diagnostics fields i.e. l_lim and phy_chl are defined in
C GCHEM_FIELDS.h

C Chl:C
      _RL sp_chl2c_min
      _RL sp_chl2c_max
      _RL lp_chl2c_min
      _RL lp_chl2c_max
      _RL alpha_chl

C Other 
      _RL     gchem_ForcingPeriod
      _RL     gchem_ForcingCycle

#endif /* ALLOW_GCHEM */
