C $Header: /u/gcmpack/MITgcm/pkg/gchem/GCHEM_FIELDS.h,v 1.1 2004/11/28 23:48:31 mlosch Exp $
C $Name:  $

#ifdef ALLOW_GCHEM
CBOP
C    !ROUTINE: GCHEM_FIELDS.h
C    !INTERFACE:
 
C    !DESCRIPTION:
C Contains tracer fields specifically for chemical tracers.
C
C  gchemTendency :: 3DxPTRACER_num field that store the tendencies due
C                   to the bio-geochemical model

#ifndef GCHEM_SEPARATE_FORCING
      _RL gchemTendency(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,
     &                  PTRACERS_num)

      COMMON /GCHEM_FIELDS/ 
     &     gchemTendency
#endif /* GCHEM_SEPARATE_FORCING */

      _RL phy_chl(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL l_lim_sp_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL l_lim_lp_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL tmp_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      _RL phy_chl_day(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL l_lim_sp_field_day(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL l_lim_lp_field_day(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL tmp_field_day(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      _RL sz_graz_sp_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL lz_graz_lp_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL lz_graz_sz_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL pz_graz_lp_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL pz_graz_sz_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL pz_graz_lz_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      _RL chl2c_sp_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL chl2c_lp_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      _RL sp_NPP_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL lp_NPP_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
	
      _RL sp_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL lp_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL sz_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL lz_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL pz_field(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      _RL pTracer_tavg_field (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,
     &              PTRACERS_num)

      COMMON /GCHEM_FIELDS/
     &     l_lim_sp_field, l_lim_lp_field, phy_chl, tmp_field,
     &     sz_graz_sp_field, lz_graz_lp_field, lz_graz_sz_field, 
     &     pz_graz_lp_field, pz_graz_sz_field, pz_graz_lz_field,
     &     chl2c_sp_field,chl2c_lp_field, phy_chl_day, 
     &     pTracer_tavg_field, sp_NPP_field, lp_NPP_field,
     &     l_lim_sp_field_day, l_lim_lp_field_day, tmp_field_day, 
     &     sp_field,  lp_field,  sz_field,  lz_field,  pz_field

CEOP
#endif /* ALLOW_GCHEM */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
