#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=icx
CCC=icx
CXX=icx
FC=ifort
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/c_src/compfg.o \
	${OBJECTDIR}/c_src/fock2.o \
	${OBJECTDIR}/c_src/fockdorbs.o \
	${OBJECTDIR}/c_src/getdat.o \
	${OBJECTDIR}/c_src/gettxt.o \
	${OBJECTDIR}/c_src/helect.o \
	${OBJECTDIR}/c_src/mopac.o \
	${OBJECTDIR}/c_src/readmo.o \
	${OBJECTDIR}/c_src/run_mopac.o \
	${OBJECTDIR}/f_src/Common_arrays_C.o \
	${OBJECTDIR}/f_src/H_bond_correction_EC_plus_ER.o \
	${OBJECTDIR}/f_src/H_bond_correction_EH_plus.o \
	${OBJECTDIR}/f_src/H_bond_correction_PM6_DH_Dispersion.o \
	${OBJECTDIR}/f_src/H_bond_correction_PM6_DH_type.o \
	${OBJECTDIR}/f_src/H_bond_correction_bits.o \
	${OBJECTDIR}/f_src/H_bonds4.o \
	${OBJECTDIR}/f_src/MOZYME_C.o \
	${OBJECTDIR}/f_src/aababc.o \
	${OBJECTDIR}/f_src/aababc_I.o \
	${OBJECTDIR}/f_src/aabacd.o \
	${OBJECTDIR}/f_src/aabacd_I.o \
	${OBJECTDIR}/f_src/aabbcd.o \
	${OBJECTDIR}/f_src/aabbcd_I.o \
	${OBJECTDIR}/f_src/add_more_interactions.o \
	${OBJECTDIR}/f_src/addfck_I.o \
	${OBJECTDIR}/f_src/addhb.o \
	${OBJECTDIR}/f_src/addhcr_I.o \
	${OBJECTDIR}/f_src/addnuc_I.o \
	${OBJECTDIR}/f_src/adjvec.o \
	${OBJECTDIR}/f_src/afmm_mod.o \
	${OBJECTDIR}/f_src/aintgs_I.o \
	${OBJECTDIR}/f_src/alphaf_I.o \
	${OBJECTDIR}/f_src/analyt.o \
	${OBJECTDIR}/f_src/analyt_C.o \
	${OBJECTDIR}/f_src/anavib.o \
	${OBJECTDIR}/f_src/asum_I.o \
	${OBJECTDIR}/f_src/atomrs.o \
	${OBJECTDIR}/f_src/aval_I.o \
	${OBJECTDIR}/f_src/babbbc.o \
	${OBJECTDIR}/f_src/babbbc_I.o \
	${OBJECTDIR}/f_src/babbcd.o \
	${OBJECTDIR}/f_src/babbcd_I.o \
	${OBJECTDIR}/f_src/bangle.o \
	${OBJECTDIR}/f_src/bangle_I.o \
	${OBJECTDIR}/f_src/bfn.o \
	${OBJECTDIR}/f_src/bfn_I.o \
	${OBJECTDIR}/f_src/bintgs_I.o \
	${OBJECTDIR}/f_src/blas_I.o \
	${OBJECTDIR}/f_src/buildf.o \
	${OBJECTDIR}/f_src/calpar.o \
	${OBJECTDIR}/f_src/calpar_I.o \
	${OBJECTDIR}/f_src/ccprod_I.o \
	${OBJECTDIR}/f_src/ccrep.o \
	${OBJECTDIR}/f_src/cdiag_I.o \
	${OBJECTDIR}/f_src/chanel_C.o \
	${OBJECTDIR}/f_src/charg_I.o \
	${OBJECTDIR}/f_src/check.o \
	${OBJECTDIR}/f_src/chkion.o \
	${OBJECTDIR}/f_src/chklew.o \
	${OBJECTDIR}/f_src/chrge.o \
	${OBJECTDIR}/f_src/chrge_I.o \
	${OBJECTDIR}/f_src/chrge_for_MOZYME.o \
	${OBJECTDIR}/f_src/cnvg.o \
	${OBJECTDIR}/f_src/cnvg_I.o \
	${OBJECTDIR}/f_src/cnvgz.o \
	${OBJECTDIR}/f_src/coe.o \
	${OBJECTDIR}/f_src/coe_I.o \
	${OBJECTDIR}/f_src/collid_I.o \
	${OBJECTDIR}/f_src/collis_I.o \
	${OBJECTDIR}/f_src/compct.o \
	${OBJECTDIR}/f_src/compfg_I.o \
	${OBJECTDIR}/f_src/conref_C.o \
	${OBJECTDIR}/f_src/convert_storage.o \
	${OBJECTDIR}/f_src/cosmo.o \
	${OBJECTDIR}/f_src/cosmo_C.o \
	${OBJECTDIR}/f_src/dang.o \
	${OBJECTDIR}/f_src/dang_I.o \
	${OBJECTDIR}/f_src/datin.o \
	${OBJECTDIR}/f_src/datin_I.o \
	${OBJECTDIR}/f_src/daxpy_I.o \
	${OBJECTDIR}/f_src/dcart.o \
	${OBJECTDIR}/f_src/dcart_I.o \
	${OBJECTDIR}/f_src/ddot_I.o \
	${OBJECTDIR}/f_src/ddpo_I.o \
	${OBJECTDIR}/f_src/delmol_I.o \
	${OBJECTDIR}/f_src/delri_I.o \
	${OBJECTDIR}/f_src/delsta.o \
	${OBJECTDIR}/f_src/densit.o \
	${OBJECTDIR}/f_src/density_cuda_i.o \
	${OBJECTDIR}/f_src/density_for_GPU.o \
	${OBJECTDIR}/f_src/density_for_MOZYME.o \
	${OBJECTDIR}/f_src/deri0.o \
	${OBJECTDIR}/f_src/deri0_I.o \
	${OBJECTDIR}/f_src/deri1.o \
	${OBJECTDIR}/f_src/deri1_I.o \
	${OBJECTDIR}/f_src/deri2.o \
	${OBJECTDIR}/f_src/deri21.o \
	${OBJECTDIR}/f_src/deri21_I.o \
	${OBJECTDIR}/f_src/deri22.o \
	${OBJECTDIR}/f_src/deri22_I.o \
	${OBJECTDIR}/f_src/deri23.o \
	${OBJECTDIR}/f_src/deri23_I.o \
	${OBJECTDIR}/f_src/deri2_I.o \
	${OBJECTDIR}/f_src/deritr.o \
	${OBJECTDIR}/f_src/deritr_I.o \
	${OBJECTDIR}/f_src/deriv.o \
	${OBJECTDIR}/f_src/deriv_I.o \
	${OBJECTDIR}/f_src/dernvo.o \
	${OBJECTDIR}/f_src/dernvo_I.o \
	${OBJECTDIR}/f_src/ders_I.o \
	${OBJECTDIR}/f_src/dex2_I.o \
	${OBJECTDIR}/f_src/dfield.o \
	${OBJECTDIR}/f_src/dfield_I.o \
	${OBJECTDIR}/f_src/dfock2.o \
	${OBJECTDIR}/f_src/dfock2_I.o \
	${OBJECTDIR}/f_src/dfport.o \
	${OBJECTDIR}/f_src/dfpsav.o \
	${OBJECTDIR}/f_src/dfpsav_I.o \
	${OBJECTDIR}/f_src/dftd3_bits.o \
	${OBJECTDIR}/f_src/dgefa_I.o \
	${OBJECTDIR}/f_src/dgetri_I.o \
	${OBJECTDIR}/f_src/dhc.o \
	${OBJECTDIR}/f_src/dhc_I.o \
	${OBJECTDIR}/f_src/dhcore.o \
	${OBJECTDIR}/f_src/dhcore_I.o \
	${OBJECTDIR}/f_src/diag.o \
	${OBJECTDIR}/f_src/diag_I.o \
	${OBJECTDIR}/f_src/diag_for_GPU.o \
	${OBJECTDIR}/f_src/diagg.o \
	${OBJECTDIR}/f_src/diagg1.o \
	${OBJECTDIR}/f_src/diagg2.o \
	${OBJECTDIR}/f_src/diagi.o \
	${OBJECTDIR}/f_src/diagi_I.o \
	${OBJECTDIR}/f_src/diat.o \
	${OBJECTDIR}/f_src/diat2_I.o \
	${OBJECTDIR}/f_src/diat_I.o \
	${OBJECTDIR}/f_src/diegrd_I.o \
	${OBJECTDIR}/f_src/dielen_I.o \
	${OBJECTDIR}/f_src/digit.o \
	${OBJECTDIR}/f_src/digit_I.o \
	${OBJECTDIR}/f_src/dihed.o \
	${OBJECTDIR}/f_src/dihed_I.o \
	${OBJECTDIR}/f_src/dijkl1.o \
	${OBJECTDIR}/f_src/dijkl1_I.o \
	${OBJECTDIR}/f_src/dijkl2.o \
	${OBJECTDIR}/f_src/dijkl2_I.o \
	${OBJECTDIR}/f_src/dijkld_I.o \
	${OBJECTDIR}/f_src/dimens.o \
	${OBJECTDIR}/f_src/dimens_I.o \
	${OBJECTDIR}/f_src/dipole_for_MOZYME.o \
	${OBJECTDIR}/f_src/disp_DnX.o \
	${OBJECTDIR}/f_src/dist2.o \
	${OBJECTDIR}/f_src/dist2_I.o \
	${OBJECTDIR}/f_src/dmecip_I.o \
	${OBJECTDIR}/f_src/dnrm2_I.o \
	${OBJECTDIR}/f_src/dofs.o \
	${OBJECTDIR}/f_src/dofs_I.o \
	${OBJECTDIR}/f_src/dopen_I.o \
	${OBJECTDIR}/f_src/dot.o \
	${OBJECTDIR}/f_src/dot_I.o \
	${OBJECTDIR}/f_src/drc.o \
	${OBJECTDIR}/f_src/drc_I.o \
	${OBJECTDIR}/f_src/drcout.o \
	${OBJECTDIR}/f_src/drcout_I.o \
	${OBJECTDIR}/f_src/drepp2_I.o \
	${OBJECTDIR}/f_src/drotat_I.o \
	${OBJECTDIR}/f_src/dsum_I.o \
	${OBJECTDIR}/f_src/dtran2.o \
	${OBJECTDIR}/f_src/dtran2_I.o \
	${OBJECTDIR}/f_src/dtrans.o \
	${OBJECTDIR}/f_src/dtrmm_I.o \
	${OBJECTDIR}/f_src/dtrmv_I.o \
	${OBJECTDIR}/f_src/dtrti2_I.o \
	${OBJECTDIR}/f_src/dvfill_I.o \
	${OBJECTDIR}/f_src/ef_C.o \
	${OBJECTDIR}/f_src/eigen.o \
	${OBJECTDIR}/f_src/eigenvectors_LAPACK.o \
	${OBJECTDIR}/f_src/eimp.o \
	${OBJECTDIR}/f_src/einvit_I.o \
	${OBJECTDIR}/f_src/eiscor_I.o \
	${OBJECTDIR}/f_src/elau_I.o \
	${OBJECTDIR}/f_src/elemts_C.o \
	${OBJECTDIR}/f_src/elenuc_I.o \
	${OBJECTDIR}/f_src/elesn_I.o \
	${OBJECTDIR}/f_src/en_I.o \
	${OBJECTDIR}/f_src/enpart.o \
	${OBJECTDIR}/f_src/enpart_I.o \
	${OBJECTDIR}/f_src/epsab_I.o \
	${OBJECTDIR}/f_src/epseta.o \
	${OBJECTDIR}/f_src/epseta_I.o \
	${OBJECTDIR}/f_src/epslon_I.o \
	${OBJECTDIR}/f_src/eqlrat_I.o \
	${OBJECTDIR}/f_src/esn_I.o \
	${OBJECTDIR}/f_src/esp.o \
	${OBJECTDIR}/f_src/esp1_I.o \
	${OBJECTDIR}/f_src/esp_C.o \
	${OBJECTDIR}/f_src/esp_utilities.o \
	${OBJECTDIR}/f_src/espfit_I.o \
	${OBJECTDIR}/f_src/estpi1_I.o \
	${OBJECTDIR}/f_src/etime_I.o \
	${OBJECTDIR}/f_src/etrbk3_I.o \
	${OBJECTDIR}/f_src/etred3_I.o \
	${OBJECTDIR}/f_src/euler_C.o \
	${OBJECTDIR}/f_src/evvrsp_I.o \
	${OBJECTDIR}/f_src/exchng.o \
	${OBJECTDIR}/f_src/exchng_I.o \
	${OBJECTDIR}/f_src/fbx_I.o \
	${OBJECTDIR}/f_src/fcnpp_I.o \
	${OBJECTDIR}/f_src/ffreq1_I.o \
	${OBJECTDIR}/f_src/ffreq2_I.o \
	${OBJECTDIR}/f_src/fhpatn_I.o \
	${OBJECTDIR}/f_src/fillij.o \
	${OBJECTDIR}/f_src/findn1.o \
	${OBJECTDIR}/f_src/finish.o \
	${OBJECTDIR}/f_src/flepo.o \
	${OBJECTDIR}/f_src/flepo_I.o \
	${OBJECTDIR}/f_src/fmat.o \
	${OBJECTDIR}/f_src/fock1.o \
	${OBJECTDIR}/f_src/fock1_I.o \
	${OBJECTDIR}/f_src/fock1_for_MOZYME.o \
	${OBJECTDIR}/f_src/fock2.o \
	${OBJECTDIR}/f_src/fock2_I.o \
	${OBJECTDIR}/f_src/fock2z.o \
	${OBJECTDIR}/f_src/fockd2_I.o \
	${OBJECTDIR}/f_src/fordd_I.o \
	${OBJECTDIR}/f_src/formxy.o \
	${OBJECTDIR}/f_src/formxy_I.o \
	${OBJECTDIR}/f_src/forsav.o \
	${OBJECTDIR}/f_src/frame.o \
	${OBJECTDIR}/f_src/frame_I.o \
	${OBJECTDIR}/f_src/freda_I.o \
	${OBJECTDIR}/f_src/fsub_I.o \
	${OBJECTDIR}/f_src/funcon_C.o \
	${OBJECTDIR}/f_src/genun.o \
	${OBJECTDIR}/f_src/genun_I.o \
	${OBJECTDIR}/f_src/genvec_I.o \
	${OBJECTDIR}/f_src/geochk.o \
	${OBJECTDIR}/f_src/getdat.o \
	${OBJECTDIR}/f_src/getgeg.o \
	${OBJECTDIR}/f_src/getgeg_I.o \
	${OBJECTDIR}/f_src/getgeo.o \
	${OBJECTDIR}/f_src/getgeo_I.o \
	${OBJECTDIR}/f_src/getpdb.o \
	${OBJECTDIR}/f_src/getsym.o \
	${OBJECTDIR}/f_src/getsym_I.o \
	${OBJECTDIR}/f_src/gettxt.o \
	${OBJECTDIR}/f_src/gettxt_I.o \
	${OBJECTDIR}/f_src/getval.o \
	${OBJECTDIR}/f_src/getval_I.o \
	${OBJECTDIR}/f_src/gmetry.o \
	${OBJECTDIR}/f_src/gmetry_I.o \
	${OBJECTDIR}/f_src/gover_I.o \
	${OBJECTDIR}/f_src/greek.o \
	${OBJECTDIR}/f_src/grids_I.o \
	${OBJECTDIR}/f_src/gstore_I.o \
	${OBJECTDIR}/f_src/h1elec.o \
	${OBJECTDIR}/f_src/h1elec_I.o \
	${OBJECTDIR}/f_src/haddon_I.o \
	${OBJECTDIR}/f_src/hbonds.o \
	${OBJECTDIR}/f_src/hcore.o \
	${OBJECTDIR}/f_src/hcore_I.o \
	${OBJECTDIR}/f_src/hcore_for_MOZYME.o \
	${OBJECTDIR}/f_src/hcored_I.o \
	${OBJECTDIR}/f_src/helect.o \
	${OBJECTDIR}/f_src/helect_I.o \
	${OBJECTDIR}/f_src/helecz.o \
	${OBJECTDIR}/f_src/hmuf_I.o \
	${OBJECTDIR}/f_src/hplusf_I.o \
	${OBJECTDIR}/f_src/hybrid.o \
	${OBJECTDIR}/f_src/ijbo.o \
	${OBJECTDIR}/f_src/inid_I.o \
	${OBJECTDIR}/f_src/inighd_I.o \
	${OBJECTDIR}/f_src/init_filenames.o \
	${OBJECTDIR}/f_src/initsn_I.o \
	${OBJECTDIR}/f_src/initsv_I.o \
	${OBJECTDIR}/f_src/insymc_I.o \
	${OBJECTDIR}/f_src/interp.o \
	${OBJECTDIR}/f_src/interp_I.o \
	${OBJECTDIR}/f_src/ionout.o \
	${OBJECTDIR}/f_src/ird_I.o \
	${OBJECTDIR}/f_src/isitsc.o \
	${OBJECTDIR}/f_src/iten_I.o \
	${OBJECTDIR}/f_src/iter_C.o \
	${OBJECTDIR}/f_src/iter_for_MOZYME.o \
	${OBJECTDIR}/f_src/jab.o \
	${OBJECTDIR}/f_src/jab_I.o \
	${OBJECTDIR}/f_src/jab_for_MOZYME.o \
	${OBJECTDIR}/f_src/jcarin.o \
	${OBJECTDIR}/f_src/jcarin_I.o \
	${OBJECTDIR}/f_src/jdate.o \
	${OBJECTDIR}/f_src/journal_references_C.o \
	${OBJECTDIR}/f_src/kab.o \
	${OBJECTDIR}/f_src/kab_I.o \
	${OBJECTDIR}/f_src/kab_for_MOZYME.o \
	${OBJECTDIR}/f_src/lbfgs.o \
	${OBJECTDIR}/f_src/lewis.o \
	${OBJECTDIR}/f_src/ligand.o \
	${OBJECTDIR}/f_src/linear_cosmo.o \
	${OBJECTDIR}/f_src/linmin.o \
	${OBJECTDIR}/f_src/linmin_I.o \
	${OBJECTDIR}/f_src/linpack.o \
	${OBJECTDIR}/f_src/local.o \
	${OBJECTDIR}/f_src/local2.o \
	${OBJECTDIR}/f_src/local_for_MOZYME.o \
	${OBJECTDIR}/f_src/lsame_I.o \
	${OBJECTDIR}/f_src/lyse.o \
	${OBJECTDIR}/f_src/makeuf_I.o \
	${OBJECTDIR}/f_src/makopr_I.o \
	${OBJECTDIR}/f_src/maksym.o \
	${OBJECTDIR}/f_src/maksym_I.o \
	${OBJECTDIR}/f_src/makvec.o \
	${OBJECTDIR}/f_src/mamult.o \
	${OBJECTDIR}/f_src/mamult_I.o \
	${OBJECTDIR}/f_src/mamult_cuda_i.o \
	${OBJECTDIR}/f_src/maps_C.o \
	${OBJECTDIR}/f_src/mat33.o \
	${OBJECTDIR}/f_src/mat33_I.o \
	${OBJECTDIR}/f_src/matout.o \
	${OBJECTDIR}/f_src/matout_I.o \
	${OBJECTDIR}/f_src/mbonds.o \
	${OBJECTDIR}/f_src/me08a_I.o \
	${OBJECTDIR}/f_src/meci_C.o \
	${OBJECTDIR}/f_src/mecid.o \
	${OBJECTDIR}/f_src/mecid_I.o \
	${OBJECTDIR}/f_src/mecih.o \
	${OBJECTDIR}/f_src/mecih_I.o \
	${OBJECTDIR}/f_src/mecip.o \
	${OBJECTDIR}/f_src/mecip_I.o \
	${OBJECTDIR}/f_src/mepchg_I.o \
	${OBJECTDIR}/f_src/mepmap_I.o \
	${OBJECTDIR}/f_src/meprot_I.o \
	${OBJECTDIR}/f_src/minv.o \
	${OBJECTDIR}/f_src/minv_I.o \
	${OBJECTDIR}/f_src/mlmo.o \
	${OBJECTDIR}/f_src/mndod.o \
	${OBJECTDIR}/f_src/mndod_C.o \
	${OBJECTDIR}/f_src/mod_atomradii.o \
	${OBJECTDIR}/f_src/mod_calls_cublas.o \
	${OBJECTDIR}/f_src/mod_gpu_info.o \
	${OBJECTDIR}/f_src/mod_vars_cuda.o \
	${OBJECTDIR}/f_src/modchg.o \
	${OBJECTDIR}/f_src/modgra.o \
	${OBJECTDIR}/f_src/moldat.o \
	${OBJECTDIR}/f_src/molkst_C.o \
	${OBJECTDIR}/f_src/molmec_C.o \
	${OBJECTDIR}/f_src/molval.o \
	${OBJECTDIR}/f_src/molval_I.o \
	${OBJECTDIR}/f_src/mopend.o \
	${OBJECTDIR}/f_src/mopend_I.o \
	${OBJECTDIR}/f_src/mpcbds_I.o \
	${OBJECTDIR}/f_src/mpcpop_I.o \
	${OBJECTDIR}/f_src/mpcsyb.o \
	${OBJECTDIR}/f_src/mtxm.o \
	${OBJECTDIR}/f_src/mtxm_I.o \
	${OBJECTDIR}/f_src/mtxmc.o \
	${OBJECTDIR}/f_src/mtxmc_I.o \
	${OBJECTDIR}/f_src/mullik.o \
	${OBJECTDIR}/f_src/mullik_I.o \
	${OBJECTDIR}/f_src/mult.o \
	${OBJECTDIR}/f_src/mult33.o \
	${OBJECTDIR}/f_src/mult33_I.o \
	${OBJECTDIR}/f_src/mult_I.o \
	${OBJECTDIR}/f_src/mxm.o \
	${OBJECTDIR}/f_src/mxm_I.o \
	${OBJECTDIR}/f_src/mxmt.o \
	${OBJECTDIR}/f_src/mxmt_I.o \
	${OBJECTDIR}/f_src/mxv.o \
	${OBJECTDIR}/f_src/mxv_I.o \
	${OBJECTDIR}/f_src/myword.o \
	${OBJECTDIR}/f_src/myword_I.o \
	${OBJECTDIR}/f_src/naican_I.o \
	${OBJECTDIR}/f_src/naicas_I.o \
	${OBJECTDIR}/f_src/names.o \
	${OBJECTDIR}/f_src/new_esp.o \
	${OBJECTDIR}/f_src/newflg.o \
	${OBJECTDIR}/f_src/ngamtg_I.o \
	${OBJECTDIR}/f_src/ngefis_I.o \
	${OBJECTDIR}/f_src/ngidri_I.o \
	${OBJECTDIR}/f_src/ngoke_I.o \
	${OBJECTDIR}/f_src/nllsn_I.o \
	${OBJECTDIR}/f_src/nonbet_I.o \
	${OBJECTDIR}/f_src/nonope_I.o \
	${OBJECTDIR}/f_src/nonor_I.o \
	${OBJECTDIR}/f_src/nuchar.o \
	${OBJECTDIR}/f_src/nuchar_I.o \
	${OBJECTDIR}/f_src/nxtmer.o \
	${OBJECTDIR}/f_src/openda_I.o \
	${OBJECTDIR}/f_src/orient_I.o \
	${OBJECTDIR}/f_src/osinv.o \
	${OBJECTDIR}/f_src/osinv_I.o \
	${OBJECTDIR}/f_src/outer1.o \
	${OBJECTDIR}/f_src/outer2.o \
	${OBJECTDIR}/f_src/overlaps_C.o \
	${OBJECTDIR}/f_src/ovlp_I.o \
	${OBJECTDIR}/f_src/packp_I.o \
	${OBJECTDIR}/f_src/parameters_C.o \
	${OBJECTDIR}/f_src/parameters_for_PM7_C.o \
	${OBJECTDIR}/f_src/parameters_for_PM7_Sparkles_C.o \
	${OBJECTDIR}/f_src/parameters_for_PM7_TS_C.o \
	${OBJECTDIR}/f_src/partxy.o \
	${OBJECTDIR}/f_src/partxy_I.o \
	${OBJECTDIR}/f_src/pdgrid_I.o \
	${OBJECTDIR}/f_src/perm.o \
	${OBJECTDIR}/f_src/perm_I.o \
	${OBJECTDIR}/f_src/picopt.o \
	${OBJECTDIR}/f_src/pinout.o \
	${OBJECTDIR}/f_src/plato_I.o \
	${OBJECTDIR}/f_src/pmep.o \
	${OBJECTDIR}/f_src/pmep1_I.o \
	${OBJECTDIR}/f_src/pmep_I.o \
	${OBJECTDIR}/f_src/pmepco_I.o \
	${OBJECTDIR}/f_src/poij_I.o \
	${OBJECTDIR}/f_src/pol_vol_I.o \
	${OBJECTDIR}/f_src/post_scf_corrections.o \
	${OBJECTDIR}/f_src/potcal_I.o \
	${OBJECTDIR}/f_src/powsav_I.o \
	${OBJECTDIR}/f_src/printp_I.o \
	${OBJECTDIR}/f_src/prtdrc.o \
	${OBJECTDIR}/f_src/prtdrc_I.o \
	${OBJECTDIR}/f_src/prtgra.o \
	${OBJECTDIR}/f_src/prthco_I.o \
	${OBJECTDIR}/f_src/prthes_I.o \
	${OBJECTDIR}/f_src/prtlmo.o \
	${OBJECTDIR}/f_src/prtpar_I.o \
	${OBJECTDIR}/f_src/prttim.o \
	${OBJECTDIR}/f_src/prttim_I.o \
	${OBJECTDIR}/f_src/pulay.o \
	${OBJECTDIR}/f_src/pulay_I.o \
	${OBJECTDIR}/f_src/quadr.o \
	${OBJECTDIR}/f_src/quadr_I.o \
	${OBJECTDIR}/f_src/react1.o \
	${OBJECTDIR}/f_src/react1_I.o \
	${OBJECTDIR}/f_src/reada.o \
	${OBJECTDIR}/f_src/reada_I.o \
	${OBJECTDIR}/f_src/redatm_I.o \
	${OBJECTDIR}/f_src/refer.o \
	${OBJECTDIR}/f_src/refer_I.o \
	${OBJECTDIR}/f_src/refkey_C.o \
	${OBJECTDIR}/f_src/reorth.o \
	${OBJECTDIR}/f_src/repp_I.o \
	${OBJECTDIR}/f_src/reppd2_I.o \
	${OBJECTDIR}/f_src/reppd_I.o \
	${OBJECTDIR}/f_src/reseq.o \
	${OBJECTDIR}/f_src/resolv.o \
	${OBJECTDIR}/f_src/resolv_I.o \
	${OBJECTDIR}/f_src/rijkl_I.o \
	${OBJECTDIR}/f_src/rotat_I.o \
	${OBJECTDIR}/f_src/rotatd_I.o \
	${OBJECTDIR}/f_src/rotate.o \
	${OBJECTDIR}/f_src/rotate_C.o \
	${OBJECTDIR}/f_src/rotate_I.o \
	${OBJECTDIR}/f_src/rotlmo.o \
	${OBJECTDIR}/f_src/rotmat_I.o \
	${OBJECTDIR}/f_src/rotmol.o \
	${OBJECTDIR}/f_src/rotmol_I.o \
	${OBJECTDIR}/f_src/rsc_I.o \
	${OBJECTDIR}/f_src/rsp.o \
	${OBJECTDIR}/f_src/scfcri.o \
	${OBJECTDIR}/f_src/schmib.o \
	${OBJECTDIR}/f_src/schmib_I.o \
	${OBJECTDIR}/f_src/schmit.o \
	${OBJECTDIR}/f_src/schmit_I.o \
	${OBJECTDIR}/f_src/scprm_I.o \
	${OBJECTDIR}/f_src/search_I.o \
	${OBJECTDIR}/f_src/second.o \
	${OBJECTDIR}/f_src/second_I.o \
	${OBJECTDIR}/f_src/selmos.o \
	${OBJECTDIR}/f_src/set.o \
	${OBJECTDIR}/f_src/set_I.o \
	${OBJECTDIR}/f_src/set_up_MOZYME_arrays.o \
	${OBJECTDIR}/f_src/set_up_RAPID.o \
	${OBJECTDIR}/f_src/set_up_dentate.o \
	${OBJECTDIR}/f_src/setup3_I.o \
	${OBJECTDIR}/f_src/setup_mopac_arrays.o \
	${OBJECTDIR}/f_src/setupg.o \
	${OBJECTDIR}/f_src/setupk.o \
	${OBJECTDIR}/f_src/solrot.o \
	${OBJECTDIR}/f_src/solrot_I.o \
	${OBJECTDIR}/f_src/sort.o \
	${OBJECTDIR}/f_src/sort_I.o \
	${OBJECTDIR}/f_src/sp_two_electron.o \
	${OBJECTDIR}/f_src/spcore_I.o \
	${OBJECTDIR}/f_src/spline_I.o \
	${OBJECTDIR}/f_src/ss_I.o \
	${OBJECTDIR}/f_src/suma2_I.o \
	${OBJECTDIR}/f_src/supdot.o \
	${OBJECTDIR}/f_src/supdot_I.o \
	${OBJECTDIR}/f_src/superd.o \
	${OBJECTDIR}/f_src/superd_I.o \
	${OBJECTDIR}/f_src/surfa_I.o \
	${OBJECTDIR}/f_src/surfac_I.o \
	${OBJECTDIR}/f_src/swap.o \
	${OBJECTDIR}/f_src/swap_I.o \
	${OBJECTDIR}/f_src/switch.o \
	${OBJECTDIR}/f_src/symdec_I.o \
	${OBJECTDIR}/f_src/symh.o \
	${OBJECTDIR}/f_src/symh_I.o \
	${OBJECTDIR}/f_src/symmetry_C.o \
	${OBJECTDIR}/f_src/symopr.o \
	${OBJECTDIR}/f_src/symopr_I.o \
	${OBJECTDIR}/f_src/symp_I.o \
	${OBJECTDIR}/f_src/sympop.o \
	${OBJECTDIR}/f_src/sympop_I.o \
	${OBJECTDIR}/f_src/symr.o \
	${OBJECTDIR}/f_src/symr_I.o \
	${OBJECTDIR}/f_src/symt.o \
	${OBJECTDIR}/f_src/symt_I.o \
	${OBJECTDIR}/f_src/symtry.o \
	${OBJECTDIR}/f_src/symtry_I.o \
	${OBJECTDIR}/f_src/tf_I.o \
	${OBJECTDIR}/f_src/tidy.o \
	${OBJECTDIR}/f_src/time_I.o \
	${OBJECTDIR}/f_src/timer.o \
	${OBJECTDIR}/f_src/timer_I.o \
	${OBJECTDIR}/f_src/timout.o \
	${OBJECTDIR}/f_src/transf_I.o \
	${OBJECTDIR}/f_src/trsub_I.o \
	${OBJECTDIR}/f_src/trudgu_I.o \
	${OBJECTDIR}/f_src/trugdu_I.o \
	${OBJECTDIR}/f_src/trugud_I.o \
	${OBJECTDIR}/f_src/tx_I.o \
	${OBJECTDIR}/f_src/txtype.o \
	${OBJECTDIR}/f_src/upcase.o \
	${OBJECTDIR}/f_src/upcase_I.o \
	${OBJECTDIR}/f_src/update.o \
	${OBJECTDIR}/f_src/update_I.o \
	${OBJECTDIR}/f_src/values.o \
	${OBJECTDIR}/f_src/vastkind.o \
	${OBJECTDIR}/f_src/vecprt.o \
	${OBJECTDIR}/f_src/vecprt_I.o \
	${OBJECTDIR}/f_src/vecprt_for_MOZYME.o \
	${OBJECTDIR}/f_src/volume.o \
	${OBJECTDIR}/f_src/volume_I.o \
	${OBJECTDIR}/f_src/w2mat_I.o \
	${OBJECTDIR}/f_src/worder_I.o \
	${OBJECTDIR}/f_src/wrdkey_I.o \
	${OBJECTDIR}/f_src/wrtkey.o \
	${OBJECTDIR}/f_src/wrtkey_I.o \
	${OBJECTDIR}/f_src/wrttxt.o \
	${OBJECTDIR}/f_src/wrttxt_I.o \
	${OBJECTDIR}/f_src/wstore_I.o \
	${OBJECTDIR}/f_src/xyzcry.o \
	${OBJECTDIR}/f_src/xyzcry_I.o \
	${OBJECTDIR}/f_src/xyzint.o \
	${OBJECTDIR}/f_src/xyzint_I.o \
	${OBJECTDIR}/f_src/zerom_I.o


# C Compiler Flags
CFLAGS=-lifcore -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread -lm

# CC Compiler Flags
CCFLAGS=-lifcore -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread -lm
CXXFLAGS=-lifcore -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread -lm

# Fortran Compiler Flags
FFLAGS=-lpthread -lstdc++ -I/${OBJECTDIR} -module ${OBJECTDIR} -fPIE

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libcmopac.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libcmopac.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libcmopac.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libcmopac.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libcmopac.a

${OBJECTDIR}/c_src/compfg.o: c_src/compfg.c ${OBJECTDIR}/f_src/chanel_C.o ${OBJECTDIR}/f_src/Common_arrays_C.o ${OBJECTDIR}/f_src/cosmo_C.o ${OBJECTDIR}/f_src/deriv_I.o ${OBJECTDIR}/f_src/dihed_I.o ${OBJECTDIR}/f_src/dot_I.o ${OBJECTDIR}/f_src/elemts_C.o ${OBJECTDIR}/f_src/funcon_C.o ${OBJECTDIR}/f_src/gmetry_I.o ${OBJECTDIR}/f_src/hcore_I.o ${OBJECTDIR}/f_src/mecip_I.o ${OBJECTDIR}/f_src/molkst_C.o ${OBJECTDIR}/f_src/molmec_C.o ${OBJECTDIR}/f_src/prtpar_I.o ${OBJECTDIR}/f_src/symtry_I.o ${OBJECTDIR}/f_src/timer_I.o ${OBJECTDIR}/f_src/vastkind.o ${OBJECTDIR}/f_src/volume_I.o ${OBJECTDIR}/f_src/linear_cosmo.o
	${MKDIR} -p ${OBJECTDIR}/c_src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/c_src/compfg.o c_src/compfg.c

${OBJECTDIR}/c_src/fock2.o: c_src/fock2.c
	${MKDIR} -p ${OBJECTDIR}/c_src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/c_src/fock2.o c_src/fock2.c

${OBJECTDIR}/c_src/fockdorbs.o: c_src/fockdorbs.c
	${MKDIR} -p ${OBJECTDIR}/c_src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/c_src/fockdorbs.o c_src/fockdorbs.c

${OBJECTDIR}/c_src/getdat.o: c_src/getdat.c
	${MKDIR} -p ${OBJECTDIR}/c_src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/c_src/getdat.o c_src/getdat.c

${OBJECTDIR}/c_src/gettxt.o: c_src/gettxt.c
	${MKDIR} -p ${OBJECTDIR}/c_src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/c_src/gettxt.o c_src/gettxt.c

${OBJECTDIR}/c_src/helect.o: c_src/helect.c
	${MKDIR} -p ${OBJECTDIR}/c_src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/c_src/helect.o c_src/helect.c

${OBJECTDIR}/c_src/mopac.o: c_src/mopac.c
	${MKDIR} -p ${OBJECTDIR}/c_src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/c_src/mopac.o c_src/mopac.c

${OBJECTDIR}/c_src/readmo.o: c_src/readmo.c
	${MKDIR} -p ${OBJECTDIR}/c_src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/c_src/readmo.o c_src/readmo.c

${OBJECTDIR}/c_src/run_mopac.o: c_src/run_mopac.c ${OBJECTDIR}/f_src/calpar_I.o ${OBJECTDIR}/f_src/moldat.o ${OBJECTDIR}/f_src/chanel_C.o ${OBJECTDIR}/f_src/geochk.o ${OBJECTDIR}/f_src/Common_arrays_C.o ${OBJECTDIR}/f_src/compfg_I.o ${OBJECTDIR}/f_src/cosmo_C.o ${OBJECTDIR}/f_src/datin_I.o ${OBJECTDIR}/f_src/drc_I.o ${OBJECTDIR}/f_src/fbx_I.o ${OBJECTDIR}/f_src/flepo_I.o ${OBJECTDIR}/f_src/fordd_I.o ${OBJECTDIR}/f_src/funcon_C.o ${OBJECTDIR}/f_src/maps_C.o ${OBJECTDIR}/f_src/molkst_C.o ${OBJECTDIR}/f_src/parameters_C.o ${OBJECTDIR}/f_src/pmep_I.o ${OBJECTDIR}/f_src/react1_I.o ${OBJECTDIR}/f_src/second_I.o ${OBJECTDIR}/f_src/symmetry_C.o ${OBJECTDIR}/f_src/vastkind.o ${OBJECTDIR}/f_src/parameters_for_PM7_Sparkles_C.o ${OBJECTDIR}/f_src/wrttxt_I.o
	${MKDIR} -p ${OBJECTDIR}/c_src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/c_src/run_mopac.o c_src/run_mopac.c

${OBJECTDIR}/f_src/Common_arrays_C.o: f_src/Common_arrays_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/Common_arrays_C.o f_src/Common_arrays_C.F90

${OBJECTDIR}/f_src/H_bond_correction_EC_plus_ER.o: f_src/H_bond_correction_EC_plus_ER.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/H_bond_correction_EC_plus_ER.o f_src/H_bond_correction_EC_plus_ER.F90

${OBJECTDIR}/f_src/H_bond_correction_EH_plus.o: f_src/H_bond_correction_EH_plus.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/H_bond_correction_EH_plus.o f_src/H_bond_correction_EH_plus.F90

${OBJECTDIR}/f_src/H_bond_correction_PM6_DH_Dispersion.o: f_src/H_bond_correction_PM6_DH_Dispersion.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/H_bond_correction_PM6_DH_Dispersion.o f_src/H_bond_correction_PM6_DH_Dispersion.F90

${OBJECTDIR}/f_src/H_bond_correction_PM6_DH_type.o: f_src/H_bond_correction_PM6_DH_type.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/H_bond_correction_PM6_DH_type.o f_src/H_bond_correction_PM6_DH_type.F90

${OBJECTDIR}/f_src/H_bond_correction_bits.o: f_src/H_bond_correction_bits.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/H_bond_correction_bits.o f_src/H_bond_correction_bits.F90

${OBJECTDIR}/f_src/H_bonds4.o: f_src/H_bonds4.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/H_bonds4.o f_src/H_bonds4.F90

${OBJECTDIR}/f_src/MOZYME_C.o: f_src/MOZYME_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/MOZYME_C.o f_src/MOZYME_C.F90

${OBJECTDIR}/f_src/aababc.o: f_src/aababc.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/aababc.o f_src/aababc.F90

${OBJECTDIR}/f_src/aababc_I.o: f_src/aababc_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/aababc_I.o f_src/aababc_I.F90

${OBJECTDIR}/f_src/aabacd.o: f_src/aabacd.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/aabacd.o f_src/aabacd.F90

${OBJECTDIR}/f_src/aabacd_I.o: f_src/aabacd_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/aabacd_I.o f_src/aabacd_I.F90

${OBJECTDIR}/f_src/aabbcd.o: f_src/aabbcd.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/aabbcd.o f_src/aabbcd.F90

${OBJECTDIR}/f_src/aabbcd_I.o: f_src/aabbcd_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/aabbcd_I.o f_src/aabbcd_I.F90

${OBJECTDIR}/f_src/add_more_interactions.o: f_src/add_more_interactions.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/add_more_interactions.o f_src/add_more_interactions.F90

${OBJECTDIR}/f_src/addfck_I.o: f_src/addfck_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/addfck_I.o f_src/addfck_I.F90

${OBJECTDIR}/f_src/addhb.o: f_src/addhb.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/addhb.o f_src/addhb.F90

${OBJECTDIR}/f_src/addhcr_I.o: f_src/addhcr_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/addhcr_I.o f_src/addhcr_I.F90

${OBJECTDIR}/f_src/addnuc_I.o: f_src/addnuc_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/addnuc_I.o f_src/addnuc_I.F90

${OBJECTDIR}/f_src/adjvec.o: f_src/adjvec.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/adjvec.o f_src/adjvec.F90

${OBJECTDIR}/f_src/afmm_mod.o: f_src/afmm_mod.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/afmm_mod.o f_src/afmm_mod.F90

${OBJECTDIR}/f_src/aintgs_I.o: f_src/aintgs_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/aintgs_I.o f_src/aintgs_I.F90

${OBJECTDIR}/f_src/alphaf_I.o: f_src/alphaf_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/alphaf_I.o f_src/alphaf_I.F90

${OBJECTDIR}/f_src/analyt.o: f_src/analyt.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/analyt.o f_src/analyt.F90

${OBJECTDIR}/f_src/analyt_C.o: f_src/analyt_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/analyt_C.o f_src/analyt_C.F90

${OBJECTDIR}/f_src/anavib.o: f_src/anavib.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/anavib.o f_src/anavib.F90

${OBJECTDIR}/f_src/asum_I.o: f_src/asum_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/asum_I.o f_src/asum_I.F90

${OBJECTDIR}/f_src/atomrs.o: f_src/atomrs.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/atomrs.o f_src/atomrs.F90

${OBJECTDIR}/f_src/aval_I.o: f_src/aval_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/aval_I.o f_src/aval_I.F90

${OBJECTDIR}/f_src/babbbc.o: f_src/babbbc.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/babbbc.o f_src/babbbc.F90

${OBJECTDIR}/f_src/babbbc_I.o: f_src/babbbc_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/babbbc_I.o f_src/babbbc_I.F90

${OBJECTDIR}/f_src/babbcd.o: f_src/babbcd.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/babbcd.o f_src/babbcd.F90

${OBJECTDIR}/f_src/babbcd_I.o: f_src/babbcd_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/babbcd_I.o f_src/babbcd_I.F90

${OBJECTDIR}/f_src/bangle.o: f_src/bangle.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/bangle.o f_src/bangle.F90

${OBJECTDIR}/f_src/bangle_I.o: f_src/bangle_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/bangle_I.o f_src/bangle_I.F90

${OBJECTDIR}/f_src/bfn.o: f_src/bfn.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/bfn.o f_src/bfn.F90

${OBJECTDIR}/f_src/bfn_I.o: f_src/bfn_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/bfn_I.o f_src/bfn_I.F90

${OBJECTDIR}/f_src/bintgs_I.o: f_src/bintgs_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/bintgs_I.o f_src/bintgs_I.F90

${OBJECTDIR}/f_src/blas_I.o: f_src/blas_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/blas_I.o f_src/blas_I.F90

${OBJECTDIR}/f_src/buildf.o: f_src/buildf.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/buildf.o f_src/buildf.F90

${OBJECTDIR}/f_src/calpar.o: f_src/calpar.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/calpar.o f_src/calpar.F90

${OBJECTDIR}/f_src/calpar_I.o: f_src/calpar_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/calpar_I.o f_src/calpar_I.F90

${OBJECTDIR}/f_src/ccprod_I.o: f_src/ccprod_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ccprod_I.o f_src/ccprod_I.F90

${OBJECTDIR}/f_src/ccrep.o: f_src/ccrep.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ccrep.o f_src/ccrep.F90

${OBJECTDIR}/f_src/cdiag_I.o: f_src/cdiag_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/cdiag_I.o f_src/cdiag_I.F90

${OBJECTDIR}/f_src/chanel_C.o: f_src/chanel_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/chanel_C.o f_src/chanel_C.F90

${OBJECTDIR}/f_src/charg_I.o: f_src/charg_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/charg_I.o f_src/charg_I.F90

${OBJECTDIR}/f_src/check.o: f_src/check.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/check.o f_src/check.F90

${OBJECTDIR}/f_src/chkion.o: f_src/chkion.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/chkion.o f_src/chkion.F90

${OBJECTDIR}/f_src/chklew.o: f_src/chklew.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/chklew.o f_src/chklew.F90

${OBJECTDIR}/f_src/chrge.o: f_src/chrge.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/chrge.o f_src/chrge.F90

${OBJECTDIR}/f_src/chrge_I.o: f_src/chrge_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/chrge_I.o f_src/chrge_I.F90

${OBJECTDIR}/f_src/chrge_for_MOZYME.o: f_src/chrge_for_MOZYME.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/chrge_for_MOZYME.o f_src/chrge_for_MOZYME.F90

${OBJECTDIR}/f_src/cnvg.o: f_src/cnvg.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/cnvg.o f_src/cnvg.F90

${OBJECTDIR}/f_src/cnvg_I.o: f_src/cnvg_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/cnvg_I.o f_src/cnvg_I.F90

${OBJECTDIR}/f_src/cnvgz.o: f_src/cnvgz.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/cnvgz.o f_src/cnvgz.F90

${OBJECTDIR}/f_src/coe.o: f_src/coe.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/coe.o f_src/coe.F90

${OBJECTDIR}/f_src/coe_I.o: f_src/coe_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/coe_I.o f_src/coe_I.F90

${OBJECTDIR}/f_src/collid_I.o: f_src/collid_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/collid_I.o f_src/collid_I.F90

${OBJECTDIR}/f_src/collis_I.o: f_src/collis_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/collis_I.o f_src/collis_I.F90

${OBJECTDIR}/f_src/compct.o: f_src/compct.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/compct.o f_src/compct.F90

${OBJECTDIR}/f_src/compfg_I.o: f_src/compfg_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/compfg_I.o f_src/compfg_I.F90

${OBJECTDIR}/f_src/conref_C.o: f_src/conref_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/conref_C.o f_src/conref_C.F90

${OBJECTDIR}/f_src/convert_storage.o: f_src/convert_storage.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/convert_storage.o f_src/convert_storage.F90

${OBJECTDIR}/f_src/cosmo.o: f_src/cosmo.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/cosmo.o f_src/cosmo.F90

${OBJECTDIR}/f_src/cosmo_C.o: f_src/cosmo_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/cosmo_C.o f_src/cosmo_C.F90

${OBJECTDIR}/f_src/dang.o: f_src/dang.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dang.o f_src/dang.F90

${OBJECTDIR}/f_src/dang_I.o: f_src/dang_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dang_I.o f_src/dang_I.F90

${OBJECTDIR}/f_src/datin.o: f_src/datin.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/datin.o f_src/datin.F90

${OBJECTDIR}/f_src/datin_I.o: f_src/datin_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/datin_I.o f_src/datin_I.F90

${OBJECTDIR}/f_src/daxpy_I.o: f_src/daxpy_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/daxpy_I.o f_src/daxpy_I.F90

${OBJECTDIR}/f_src/dcart.o: f_src/dcart.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dcart.o f_src/dcart.F90

${OBJECTDIR}/f_src/dcart_I.o: f_src/dcart_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dcart_I.o f_src/dcart_I.F90

${OBJECTDIR}/f_src/ddot_I.o: f_src/ddot_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ddot_I.o f_src/ddot_I.F90

${OBJECTDIR}/f_src/ddpo_I.o: f_src/ddpo_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ddpo_I.o f_src/ddpo_I.F90

${OBJECTDIR}/f_src/delmol_I.o: f_src/delmol_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/delmol_I.o f_src/delmol_I.F90

${OBJECTDIR}/f_src/delri_I.o: f_src/delri_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/delri_I.o f_src/delri_I.F90

${OBJECTDIR}/f_src/delsta.o: f_src/delsta.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/delsta.o f_src/delsta.F90

${OBJECTDIR}/f_src/densit.o: f_src/densit.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/densit.o f_src/densit.F90

${OBJECTDIR}/f_src/density_cuda_i.o: f_src/density_cuda_i.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/density_cuda_i.o f_src/density_cuda_i.F90

${OBJECTDIR}/f_src/density_for_GPU.o: f_src/density_for_GPU.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/density_for_GPU.o f_src/density_for_GPU.F90

${OBJECTDIR}/f_src/density_for_MOZYME.o: f_src/density_for_MOZYME.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/density_for_MOZYME.o f_src/density_for_MOZYME.F90

${OBJECTDIR}/f_src/deri0.o: f_src/deri0.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/deri0.o f_src/deri0.F90

${OBJECTDIR}/f_src/deri0_I.o: f_src/deri0_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/deri0_I.o f_src/deri0_I.F90

${OBJECTDIR}/f_src/deri1.o: f_src/deri1.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/deri1.o f_src/deri1.F90

${OBJECTDIR}/f_src/deri1_I.o: f_src/deri1_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/deri1_I.o f_src/deri1_I.F90

${OBJECTDIR}/f_src/deri2.o: f_src/deri2.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/deri2.o f_src/deri2.F90

${OBJECTDIR}/f_src/deri21.o: f_src/deri21.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/deri21.o f_src/deri21.F90

${OBJECTDIR}/f_src/deri21_I.o: f_src/deri21_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/deri21_I.o f_src/deri21_I.F90

${OBJECTDIR}/f_src/deri22.o: f_src/deri22.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/deri22.o f_src/deri22.F90

${OBJECTDIR}/f_src/deri22_I.o: f_src/deri22_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/deri22_I.o f_src/deri22_I.F90

${OBJECTDIR}/f_src/deri23.o: f_src/deri23.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/deri23.o f_src/deri23.F90

${OBJECTDIR}/f_src/deri23_I.o: f_src/deri23_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/deri23_I.o f_src/deri23_I.F90

${OBJECTDIR}/f_src/deri2_I.o: f_src/deri2_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/deri2_I.o f_src/deri2_I.F90

${OBJECTDIR}/f_src/deritr.o: f_src/deritr.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/deritr.o f_src/deritr.F90

${OBJECTDIR}/f_src/deritr_I.o: f_src/deritr_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/deritr_I.o f_src/deritr_I.F90

${OBJECTDIR}/f_src/deriv.o: f_src/deriv.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/deriv.o f_src/deriv.F90

${OBJECTDIR}/f_src/deriv_I.o: f_src/deriv_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/deriv_I.o f_src/deriv_I.F90

${OBJECTDIR}/f_src/dernvo.o: f_src/dernvo.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dernvo.o f_src/dernvo.F90

${OBJECTDIR}/f_src/dernvo_I.o: f_src/dernvo_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dernvo_I.o f_src/dernvo_I.F90

${OBJECTDIR}/f_src/ders_I.o: f_src/ders_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ders_I.o f_src/ders_I.F90

${OBJECTDIR}/f_src/dex2_I.o: f_src/dex2_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dex2_I.o f_src/dex2_I.F90

${OBJECTDIR}/f_src/dfield.o: f_src/dfield.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dfield.o f_src/dfield.F90

${OBJECTDIR}/f_src/dfield_I.o: f_src/dfield_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dfield_I.o f_src/dfield_I.F90

${OBJECTDIR}/f_src/dfock2.o: f_src/dfock2.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dfock2.o f_src/dfock2.F90

${OBJECTDIR}/f_src/dfock2_I.o: f_src/dfock2_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dfock2_I.o f_src/dfock2_I.F90

${OBJECTDIR}/f_src/dfport.o: f_src/dfport.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dfport.o f_src/dfport.F90

${OBJECTDIR}/f_src/dfpsav.o: f_src/dfpsav.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dfpsav.o f_src/dfpsav.F90

${OBJECTDIR}/f_src/dfpsav_I.o: f_src/dfpsav_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dfpsav_I.o f_src/dfpsav_I.F90

${OBJECTDIR}/f_src/dftd3_bits.o: f_src/dftd3_bits.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dftd3_bits.o f_src/dftd3_bits.F90

${OBJECTDIR}/f_src/dgefa_I.o: f_src/dgefa_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dgefa_I.o f_src/dgefa_I.F90

${OBJECTDIR}/f_src/dgetri_I.o: f_src/dgetri_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dgetri_I.o f_src/dgetri_I.F90

${OBJECTDIR}/f_src/dhc.o: f_src/dhc.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dhc.o f_src/dhc.F90

${OBJECTDIR}/f_src/dhc_I.o: f_src/dhc_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dhc_I.o f_src/dhc_I.F90

${OBJECTDIR}/f_src/dhcore.o: f_src/dhcore.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dhcore.o f_src/dhcore.F90

${OBJECTDIR}/f_src/dhcore_I.o: f_src/dhcore_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dhcore_I.o f_src/dhcore_I.F90

${OBJECTDIR}/f_src/diag.o: f_src/diag.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/diag.o f_src/diag.F90

${OBJECTDIR}/f_src/diag_I.o: f_src/diag_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/diag_I.o f_src/diag_I.F90

${OBJECTDIR}/f_src/diag_for_GPU.o: f_src/diag_for_GPU.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/diag_for_GPU.o f_src/diag_for_GPU.F90

${OBJECTDIR}/f_src/diagg.o: f_src/diagg.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/diagg.o f_src/diagg.F90

${OBJECTDIR}/f_src/diagg1.o: f_src/diagg1.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/diagg1.o f_src/diagg1.F90

${OBJECTDIR}/f_src/diagg2.o: f_src/diagg2.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/diagg2.o f_src/diagg2.F90

${OBJECTDIR}/f_src/diagi.o: f_src/diagi.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/diagi.o f_src/diagi.F90

${OBJECTDIR}/f_src/diagi_I.o: f_src/diagi_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/diagi_I.o f_src/diagi_I.F90

${OBJECTDIR}/f_src/diat.o: f_src/diat.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/diat.o f_src/diat.F90

${OBJECTDIR}/f_src/diat2_I.o: f_src/diat2_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/diat2_I.o f_src/diat2_I.F90

${OBJECTDIR}/f_src/diat_I.o: f_src/diat_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/diat_I.o f_src/diat_I.F90

${OBJECTDIR}/f_src/diegrd_I.o: f_src/diegrd_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/diegrd_I.o f_src/diegrd_I.F90

${OBJECTDIR}/f_src/dielen_I.o: f_src/dielen_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dielen_I.o f_src/dielen_I.F90

${OBJECTDIR}/f_src/digit.o: f_src/digit.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/digit.o f_src/digit.F90

${OBJECTDIR}/f_src/digit_I.o: f_src/digit_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/digit_I.o f_src/digit_I.F90

${OBJECTDIR}/f_src/dihed.o: f_src/dihed.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dihed.o f_src/dihed.F90

${OBJECTDIR}/f_src/dihed_I.o: f_src/dihed_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dihed_I.o f_src/dihed_I.F90

${OBJECTDIR}/f_src/dijkl1.o: f_src/dijkl1.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dijkl1.o f_src/dijkl1.F90

${OBJECTDIR}/f_src/dijkl1_I.o: f_src/dijkl1_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dijkl1_I.o f_src/dijkl1_I.F90

${OBJECTDIR}/f_src/dijkl2.o: f_src/dijkl2.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dijkl2.o f_src/dijkl2.F90

${OBJECTDIR}/f_src/dijkl2_I.o: f_src/dijkl2_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dijkl2_I.o f_src/dijkl2_I.F90

${OBJECTDIR}/f_src/dijkld_I.o: f_src/dijkld_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dijkld_I.o f_src/dijkld_I.F90

${OBJECTDIR}/f_src/dimens.o: f_src/dimens.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dimens.o f_src/dimens.F90

${OBJECTDIR}/f_src/dimens_I.o: f_src/dimens_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dimens_I.o f_src/dimens_I.F90

${OBJECTDIR}/f_src/dipole_for_MOZYME.o: f_src/dipole_for_MOZYME.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dipole_for_MOZYME.o f_src/dipole_for_MOZYME.F90

${OBJECTDIR}/f_src/disp_DnX.o: f_src/disp_DnX.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/disp_DnX.o f_src/disp_DnX.F90

${OBJECTDIR}/f_src/dist2.o: f_src/dist2.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dist2.o f_src/dist2.F90

${OBJECTDIR}/f_src/dist2_I.o: f_src/dist2_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dist2_I.o f_src/dist2_I.F90

${OBJECTDIR}/f_src/dmecip_I.o: f_src/dmecip_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dmecip_I.o f_src/dmecip_I.F90

${OBJECTDIR}/f_src/dnrm2_I.o: f_src/dnrm2_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dnrm2_I.o f_src/dnrm2_I.F90

${OBJECTDIR}/f_src/dofs.o: f_src/dofs.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dofs.o f_src/dofs.F90

${OBJECTDIR}/f_src/dofs_I.o: f_src/dofs_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dofs_I.o f_src/dofs_I.F90

${OBJECTDIR}/f_src/dopen_I.o: f_src/dopen_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dopen_I.o f_src/dopen_I.F90

${OBJECTDIR}/f_src/dot.o: f_src/dot.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dot.o f_src/dot.F90

${OBJECTDIR}/f_src/dot_I.o: f_src/dot_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dot_I.o f_src/dot_I.F90

${OBJECTDIR}/f_src/drc.o: f_src/drc.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/drc.o f_src/drc.F90

${OBJECTDIR}/f_src/drc_I.o: f_src/drc_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/drc_I.o f_src/drc_I.F90

${OBJECTDIR}/f_src/drcout.o: f_src/drcout.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/drcout.o f_src/drcout.F90

${OBJECTDIR}/f_src/drcout_I.o: f_src/drcout_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/drcout_I.o f_src/drcout_I.F90

${OBJECTDIR}/f_src/drepp2_I.o: f_src/drepp2_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/drepp2_I.o f_src/drepp2_I.F90

${OBJECTDIR}/f_src/drotat_I.o: f_src/drotat_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/drotat_I.o f_src/drotat_I.F90

${OBJECTDIR}/f_src/dsum_I.o: f_src/dsum_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dsum_I.o f_src/dsum_I.F90

${OBJECTDIR}/f_src/dtran2.o: f_src/dtran2.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dtran2.o f_src/dtran2.F90

${OBJECTDIR}/f_src/dtran2_I.o: f_src/dtran2_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dtran2_I.o f_src/dtran2_I.F90

${OBJECTDIR}/f_src/dtrans.o: f_src/dtrans.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dtrans.o f_src/dtrans.F90

${OBJECTDIR}/f_src/dtrmm_I.o: f_src/dtrmm_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dtrmm_I.o f_src/dtrmm_I.F90

${OBJECTDIR}/f_src/dtrmv_I.o: f_src/dtrmv_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dtrmv_I.o f_src/dtrmv_I.F90

${OBJECTDIR}/f_src/dtrti2_I.o: f_src/dtrti2_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dtrti2_I.o f_src/dtrti2_I.F90

${OBJECTDIR}/f_src/dvfill_I.o: f_src/dvfill_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/dvfill_I.o f_src/dvfill_I.F90

${OBJECTDIR}/f_src/ef_C.o: f_src/ef_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ef_C.o f_src/ef_C.F90

${OBJECTDIR}/f_src/eigen.o: f_src/eigen.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/eigen.o f_src/eigen.F90

${OBJECTDIR}/f_src/eigenvectors_LAPACK.o: f_src/eigenvectors_LAPACK.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/eigenvectors_LAPACK.o f_src/eigenvectors_LAPACK.F90

${OBJECTDIR}/f_src/eimp.o: f_src/eimp.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/eimp.o f_src/eimp.F90

${OBJECTDIR}/f_src/einvit_I.o: f_src/einvit_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/einvit_I.o f_src/einvit_I.F90

${OBJECTDIR}/f_src/eiscor_I.o: f_src/eiscor_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/eiscor_I.o f_src/eiscor_I.F90

${OBJECTDIR}/f_src/elau_I.o: f_src/elau_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/elau_I.o f_src/elau_I.F90

${OBJECTDIR}/f_src/elemts_C.o: f_src/elemts_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/elemts_C.o f_src/elemts_C.F90

${OBJECTDIR}/f_src/elenuc_I.o: f_src/elenuc_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/elenuc_I.o f_src/elenuc_I.F90

${OBJECTDIR}/f_src/elesn_I.o: f_src/elesn_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/elesn_I.o f_src/elesn_I.F90

${OBJECTDIR}/f_src/en_I.o: f_src/en_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/en_I.o f_src/en_I.F90

${OBJECTDIR}/f_src/enpart.o: f_src/enpart.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/enpart.o f_src/enpart.F90

${OBJECTDIR}/f_src/enpart_I.o: f_src/enpart_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/enpart_I.o f_src/enpart_I.F90

${OBJECTDIR}/f_src/epsab_I.o: f_src/epsab_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/epsab_I.o f_src/epsab_I.F90

${OBJECTDIR}/f_src/epseta.o: f_src/epseta.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/epseta.o f_src/epseta.F90

${OBJECTDIR}/f_src/epseta_I.o: f_src/epseta_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/epseta_I.o f_src/epseta_I.F90

${OBJECTDIR}/f_src/epslon_I.o: f_src/epslon_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/epslon_I.o f_src/epslon_I.F90

${OBJECTDIR}/f_src/eqlrat_I.o: f_src/eqlrat_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/eqlrat_I.o f_src/eqlrat_I.F90

${OBJECTDIR}/f_src/esn_I.o: f_src/esn_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/esn_I.o f_src/esn_I.F90

${OBJECTDIR}/f_src/esp.o: f_src/esp.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/esp.o f_src/esp.F90

${OBJECTDIR}/f_src/esp1_I.o: f_src/esp1_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/esp1_I.o f_src/esp1_I.F90

${OBJECTDIR}/f_src/esp_C.o: f_src/esp_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/esp_C.o f_src/esp_C.F90

${OBJECTDIR}/f_src/esp_utilities.o: f_src/esp_utilities.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/esp_utilities.o f_src/esp_utilities.F90

${OBJECTDIR}/f_src/espfit_I.o: f_src/espfit_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/espfit_I.o f_src/espfit_I.F90

${OBJECTDIR}/f_src/estpi1_I.o: f_src/estpi1_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/estpi1_I.o f_src/estpi1_I.F90

${OBJECTDIR}/f_src/etime_I.o: f_src/etime_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/etime_I.o f_src/etime_I.F90

${OBJECTDIR}/f_src/etrbk3_I.o: f_src/etrbk3_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/etrbk3_I.o f_src/etrbk3_I.F90

${OBJECTDIR}/f_src/etred3_I.o: f_src/etred3_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/etred3_I.o f_src/etred3_I.F90

${OBJECTDIR}/f_src/euler_C.o: f_src/euler_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/euler_C.o f_src/euler_C.F90

${OBJECTDIR}/f_src/evvrsp_I.o: f_src/evvrsp_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/evvrsp_I.o f_src/evvrsp_I.F90

${OBJECTDIR}/f_src/exchng.o: f_src/exchng.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/exchng.o f_src/exchng.F90

${OBJECTDIR}/f_src/exchng_I.o: f_src/exchng_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/exchng_I.o f_src/exchng_I.F90

${OBJECTDIR}/f_src/fbx_I.o: f_src/fbx_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/fbx_I.o f_src/fbx_I.F90

${OBJECTDIR}/f_src/fcnpp_I.o: f_src/fcnpp_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/fcnpp_I.o f_src/fcnpp_I.F90

${OBJECTDIR}/f_src/ffreq1_I.o: f_src/ffreq1_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ffreq1_I.o f_src/ffreq1_I.F90

${OBJECTDIR}/f_src/ffreq2_I.o: f_src/ffreq2_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ffreq2_I.o f_src/ffreq2_I.F90

${OBJECTDIR}/f_src/fhpatn_I.o: f_src/fhpatn_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/fhpatn_I.o f_src/fhpatn_I.F90

${OBJECTDIR}/f_src/fillij.o: f_src/fillij.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/fillij.o f_src/fillij.F90

${OBJECTDIR}/f_src/findn1.o: f_src/findn1.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/findn1.o f_src/findn1.F90

${OBJECTDIR}/f_src/finish.o: f_src/finish.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/finish.o f_src/finish.F90

${OBJECTDIR}/f_src/flepo.o: f_src/flepo.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/flepo.o f_src/flepo.F90

${OBJECTDIR}/f_src/flepo_I.o: f_src/flepo_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/flepo_I.o f_src/flepo_I.F90

${OBJECTDIR}/f_src/fmat.o: f_src/fmat.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/fmat.o f_src/fmat.F90

${OBJECTDIR}/f_src/fock1.o: f_src/fock1.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/fock1.o f_src/fock1.F90

${OBJECTDIR}/f_src/fock1_I.o: f_src/fock1_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/fock1_I.o f_src/fock1_I.F90

${OBJECTDIR}/f_src/fock1_for_MOZYME.o: f_src/fock1_for_MOZYME.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/fock1_for_MOZYME.o f_src/fock1_for_MOZYME.F90

${OBJECTDIR}/f_src/fock2.o: f_src/fock2.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/fock2.o f_src/fock2.F90

${OBJECTDIR}/f_src/fock2_I.o: f_src/fock2_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/fock2_I.o f_src/fock2_I.F90

${OBJECTDIR}/f_src/fock2z.o: f_src/fock2z.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/fock2z.o f_src/fock2z.F90

${OBJECTDIR}/f_src/fockd2_I.o: f_src/fockd2_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/fockd2_I.o f_src/fockd2_I.F90

${OBJECTDIR}/f_src/fordd_I.o: f_src/fordd_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/fordd_I.o f_src/fordd_I.F90

${OBJECTDIR}/f_src/formxy.o: f_src/formxy.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/formxy.o f_src/formxy.F90

${OBJECTDIR}/f_src/formxy_I.o: f_src/formxy_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/formxy_I.o f_src/formxy_I.F90

${OBJECTDIR}/f_src/forsav.o: f_src/forsav.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/forsav.o f_src/forsav.F90

${OBJECTDIR}/f_src/frame.o: f_src/frame.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/frame.o f_src/frame.F90

${OBJECTDIR}/f_src/frame_I.o: f_src/frame_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/frame_I.o f_src/frame_I.F90

${OBJECTDIR}/f_src/freda_I.o: f_src/freda_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/freda_I.o f_src/freda_I.F90

${OBJECTDIR}/f_src/fsub_I.o: f_src/fsub_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/fsub_I.o f_src/fsub_I.F90

${OBJECTDIR}/f_src/funcon_C.o: f_src/funcon_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/funcon_C.o f_src/funcon_C.F90

${OBJECTDIR}/f_src/genun.o: f_src/genun.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/genun.o f_src/genun.F90

${OBJECTDIR}/f_src/genun_I.o: f_src/genun_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/genun_I.o f_src/genun_I.F90

${OBJECTDIR}/f_src/genvec_I.o: f_src/genvec_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/genvec_I.o f_src/genvec_I.F90

${OBJECTDIR}/f_src/geochk.o: f_src/geochk.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/geochk.o f_src/geochk.F90

${OBJECTDIR}/f_src/getdat.o: f_src/getdat.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/getdat.o f_src/getdat.F90

${OBJECTDIR}/f_src/getgeg.o: f_src/getgeg.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/getgeg.o f_src/getgeg.F90

${OBJECTDIR}/f_src/getgeg_I.o: f_src/getgeg_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/getgeg_I.o f_src/getgeg_I.F90

${OBJECTDIR}/f_src/getgeo.o: f_src/getgeo.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/getgeo.o f_src/getgeo.F90

${OBJECTDIR}/f_src/getgeo_I.o: f_src/getgeo_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/getgeo_I.o f_src/getgeo_I.F90

${OBJECTDIR}/f_src/getpdb.o: f_src/getpdb.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/getpdb.o f_src/getpdb.F90

${OBJECTDIR}/f_src/getsym.o: f_src/getsym.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/getsym.o f_src/getsym.F90

${OBJECTDIR}/f_src/getsym_I.o: f_src/getsym_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/getsym_I.o f_src/getsym_I.F90

${OBJECTDIR}/f_src/gettxt.o: f_src/gettxt.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/gettxt.o f_src/gettxt.F90

${OBJECTDIR}/f_src/gettxt_I.o: f_src/gettxt_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/gettxt_I.o f_src/gettxt_I.F90

${OBJECTDIR}/f_src/getval.o: f_src/getval.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/getval.o f_src/getval.F90

${OBJECTDIR}/f_src/getval_I.o: f_src/getval_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/getval_I.o f_src/getval_I.F90

${OBJECTDIR}/f_src/gmetry.o: f_src/gmetry.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/gmetry.o f_src/gmetry.F90

${OBJECTDIR}/f_src/gmetry_I.o: f_src/gmetry_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/gmetry_I.o f_src/gmetry_I.F90

${OBJECTDIR}/f_src/gover_I.o: f_src/gover_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/gover_I.o f_src/gover_I.F90

${OBJECTDIR}/f_src/greek.o: f_src/greek.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/greek.o f_src/greek.F90

${OBJECTDIR}/f_src/grids_I.o: f_src/grids_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/grids_I.o f_src/grids_I.F90

${OBJECTDIR}/f_src/gstore_I.o: f_src/gstore_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/gstore_I.o f_src/gstore_I.F90

${OBJECTDIR}/f_src/h1elec.o: f_src/h1elec.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/h1elec.o f_src/h1elec.F90

${OBJECTDIR}/f_src/h1elec_I.o: f_src/h1elec_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/h1elec_I.o f_src/h1elec_I.F90

${OBJECTDIR}/f_src/haddon_I.o: f_src/haddon_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/haddon_I.o f_src/haddon_I.F90

${OBJECTDIR}/f_src/hbonds.o: f_src/hbonds.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/hbonds.o f_src/hbonds.F90

${OBJECTDIR}/f_src/hcore.o: f_src/hcore.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/hcore.o f_src/hcore.F90

${OBJECTDIR}/f_src/hcore_I.o: f_src/hcore_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/hcore_I.o f_src/hcore_I.F90

${OBJECTDIR}/f_src/hcore_for_MOZYME.o: f_src/hcore_for_MOZYME.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/hcore_for_MOZYME.o f_src/hcore_for_MOZYME.F90

${OBJECTDIR}/f_src/hcored_I.o: f_src/hcored_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/hcored_I.o f_src/hcored_I.F90

${OBJECTDIR}/f_src/helect.o: f_src/helect.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/helect.o f_src/helect.F90

${OBJECTDIR}/f_src/helect_I.o: f_src/helect_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/helect_I.o f_src/helect_I.F90

${OBJECTDIR}/f_src/helecz.o: f_src/helecz.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/helecz.o f_src/helecz.F90

${OBJECTDIR}/f_src/hmuf_I.o: f_src/hmuf_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/hmuf_I.o f_src/hmuf_I.F90

${OBJECTDIR}/f_src/hplusf_I.o: f_src/hplusf_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/hplusf_I.o f_src/hplusf_I.F90

${OBJECTDIR}/f_src/hybrid.o: f_src/hybrid.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/hybrid.o f_src/hybrid.F90

${OBJECTDIR}/f_src/ijbo.o: f_src/ijbo.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ijbo.o f_src/ijbo.F90

${OBJECTDIR}/f_src/inid_I.o: f_src/inid_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/inid_I.o f_src/inid_I.F90

${OBJECTDIR}/f_src/inighd_I.o: f_src/inighd_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/inighd_I.o f_src/inighd_I.F90

${OBJECTDIR}/f_src/init_filenames.o: f_src/init_filenames.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/init_filenames.o f_src/init_filenames.F90

${OBJECTDIR}/f_src/initsn_I.o: f_src/initsn_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/initsn_I.o f_src/initsn_I.F90

${OBJECTDIR}/f_src/initsv_I.o: f_src/initsv_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/initsv_I.o f_src/initsv_I.F90

${OBJECTDIR}/f_src/insymc_I.o: f_src/insymc_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/insymc_I.o f_src/insymc_I.F90

${OBJECTDIR}/f_src/interp.o: f_src/interp.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/interp.o f_src/interp.F90

${OBJECTDIR}/f_src/interp_I.o: f_src/interp_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/interp_I.o f_src/interp_I.F90

${OBJECTDIR}/f_src/ionout.o: f_src/ionout.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ionout.o f_src/ionout.F90

${OBJECTDIR}/f_src/ird_I.o: f_src/ird_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ird_I.o f_src/ird_I.F90

${OBJECTDIR}/f_src/isitsc.o: f_src/isitsc.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/isitsc.o f_src/isitsc.F90

${OBJECTDIR}/f_src/iten_I.o: f_src/iten_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/iten_I.o f_src/iten_I.F90

${OBJECTDIR}/f_src/iter_C.o: f_src/iter_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/iter_C.o f_src/iter_C.F90

${OBJECTDIR}/f_src/iter_for_MOZYME.o: f_src/iter_for_MOZYME.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/iter_for_MOZYME.o f_src/iter_for_MOZYME.F90

${OBJECTDIR}/f_src/jab.o: f_src/jab.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/jab.o f_src/jab.F90

${OBJECTDIR}/f_src/jab_I.o: f_src/jab_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/jab_I.o f_src/jab_I.F90

${OBJECTDIR}/f_src/jab_for_MOZYME.o: f_src/jab_for_MOZYME.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/jab_for_MOZYME.o f_src/jab_for_MOZYME.F90

${OBJECTDIR}/f_src/jcarin.o: f_src/jcarin.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/jcarin.o f_src/jcarin.F90

${OBJECTDIR}/f_src/jcarin_I.o: f_src/jcarin_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/jcarin_I.o f_src/jcarin_I.F90

${OBJECTDIR}/f_src/jdate.o: f_src/jdate.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/jdate.o f_src/jdate.F90

${OBJECTDIR}/f_src/journal_references_C.o: f_src/journal_references_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/journal_references_C.o f_src/journal_references_C.F90

${OBJECTDIR}/f_src/kab.o: f_src/kab.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/kab.o f_src/kab.F90

${OBJECTDIR}/f_src/kab_I.o: f_src/kab_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/kab_I.o f_src/kab_I.F90

${OBJECTDIR}/f_src/kab_for_MOZYME.o: f_src/kab_for_MOZYME.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/kab_for_MOZYME.o f_src/kab_for_MOZYME.F90

${OBJECTDIR}/f_src/lbfgs.o: f_src/lbfgs.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/lbfgs.o f_src/lbfgs.F90

${OBJECTDIR}/f_src/lewis.o: f_src/lewis.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/lewis.o f_src/lewis.F90

${OBJECTDIR}/f_src/ligand.o: f_src/ligand.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ligand.o f_src/ligand.F90

${OBJECTDIR}/f_src/linear_cosmo.o: f_src/linear_cosmo.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/linear_cosmo.o f_src/linear_cosmo.F90

${OBJECTDIR}/f_src/linmin.o: f_src/linmin.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/linmin.o f_src/linmin.F90

${OBJECTDIR}/f_src/linmin_I.o: f_src/linmin_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/linmin_I.o f_src/linmin_I.F90

${OBJECTDIR}/f_src/linpack.o: f_src/linpack.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/linpack.o f_src/linpack.F90

${OBJECTDIR}/f_src/local.o: f_src/local.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/local.o f_src/local.F90

${OBJECTDIR}/f_src/local2.o: f_src/local2.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/local2.o f_src/local2.F90

${OBJECTDIR}/f_src/local_for_MOZYME.o: f_src/local_for_MOZYME.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/local_for_MOZYME.o f_src/local_for_MOZYME.F90

${OBJECTDIR}/f_src/lsame_I.o: f_src/lsame_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/lsame_I.o f_src/lsame_I.F90

${OBJECTDIR}/f_src/lyse.o: f_src/lyse.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/lyse.o f_src/lyse.F90

${OBJECTDIR}/f_src/makeuf_I.o: f_src/makeuf_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/makeuf_I.o f_src/makeuf_I.F90

${OBJECTDIR}/f_src/makopr_I.o: f_src/makopr_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/makopr_I.o f_src/makopr_I.F90

${OBJECTDIR}/f_src/maksym.o: f_src/maksym.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/maksym.o f_src/maksym.F90

${OBJECTDIR}/f_src/maksym_I.o: f_src/maksym_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/maksym_I.o f_src/maksym_I.F90

${OBJECTDIR}/f_src/makvec.o: f_src/makvec.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/makvec.o f_src/makvec.F90

${OBJECTDIR}/f_src/mamult.o: f_src/mamult.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mamult.o f_src/mamult.F90

${OBJECTDIR}/f_src/mamult_I.o: f_src/mamult_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mamult_I.o f_src/mamult_I.F90

${OBJECTDIR}/f_src/mamult_cuda_i.o: f_src/mamult_cuda_i.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mamult_cuda_i.o f_src/mamult_cuda_i.F90

${OBJECTDIR}/f_src/maps_C.o: f_src/maps_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/maps_C.o f_src/maps_C.F90

${OBJECTDIR}/f_src/mat33.o: f_src/mat33.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mat33.o f_src/mat33.F90

${OBJECTDIR}/f_src/mat33_I.o: f_src/mat33_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mat33_I.o f_src/mat33_I.F90

${OBJECTDIR}/f_src/matout.o: f_src/matout.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/matout.o f_src/matout.F90

${OBJECTDIR}/f_src/matout_I.o: f_src/matout_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/matout_I.o f_src/matout_I.F90

${OBJECTDIR}/f_src/mbonds.o: f_src/mbonds.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mbonds.o f_src/mbonds.F90

${OBJECTDIR}/f_src/me08a_I.o: f_src/me08a_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/me08a_I.o f_src/me08a_I.F90

${OBJECTDIR}/f_src/meci_C.o: f_src/meci_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/meci_C.o f_src/meci_C.F90

${OBJECTDIR}/f_src/mecid.o: f_src/mecid.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mecid.o f_src/mecid.F90

${OBJECTDIR}/f_src/mecid_I.o: f_src/mecid_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mecid_I.o f_src/mecid_I.F90

${OBJECTDIR}/f_src/mecih.o: f_src/mecih.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mecih.o f_src/mecih.F90

${OBJECTDIR}/f_src/mecih_I.o: f_src/mecih_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mecih_I.o f_src/mecih_I.F90

${OBJECTDIR}/f_src/mecip.o: f_src/mecip.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mecip.o f_src/mecip.F90

${OBJECTDIR}/f_src/mecip_I.o: f_src/mecip_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mecip_I.o f_src/mecip_I.F90

${OBJECTDIR}/f_src/mepchg_I.o: f_src/mepchg_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mepchg_I.o f_src/mepchg_I.F90

${OBJECTDIR}/f_src/mepmap_I.o: f_src/mepmap_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mepmap_I.o f_src/mepmap_I.F90

${OBJECTDIR}/f_src/meprot_I.o: f_src/meprot_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/meprot_I.o f_src/meprot_I.F90

${OBJECTDIR}/f_src/minv.o: f_src/minv.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/minv.o f_src/minv.F90

${OBJECTDIR}/f_src/minv_I.o: f_src/minv_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/minv_I.o f_src/minv_I.F90

${OBJECTDIR}/f_src/mlmo.o: f_src/mlmo.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mlmo.o f_src/mlmo.F90

${OBJECTDIR}/f_src/mndod.o: f_src/mndod.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mndod.o f_src/mndod.F90

${OBJECTDIR}/f_src/mndod_C.o: f_src/mndod_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mndod_C.o f_src/mndod_C.F90

${OBJECTDIR}/f_src/mod_atomradii.o: f_src/mod_atomradii.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mod_atomradii.o f_src/mod_atomradii.F90

${OBJECTDIR}/f_src/mod_calls_cublas.o: f_src/mod_calls_cublas.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mod_calls_cublas.o f_src/mod_calls_cublas.F90

${OBJECTDIR}/f_src/mod_gpu_info.o: f_src/mod_gpu_info.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mod_gpu_info.o f_src/mod_gpu_info.F90

${OBJECTDIR}/f_src/mod_vars_cuda.o: f_src/mod_vars_cuda.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mod_vars_cuda.o f_src/mod_vars_cuda.F90

${OBJECTDIR}/f_src/modchg.o: f_src/modchg.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/modchg.o f_src/modchg.F90

${OBJECTDIR}/f_src/modgra.o: f_src/modgra.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/modgra.o f_src/modgra.F90

${OBJECTDIR}/f_src/moldat.o: f_src/moldat.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/moldat.o f_src/moldat.F90

${OBJECTDIR}/f_src/molkst_C.o: f_src/molkst_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/molkst_C.o f_src/molkst_C.F90

${OBJECTDIR}/f_src/molmec_C.o: f_src/molmec_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/molmec_C.o f_src/molmec_C.F90

${OBJECTDIR}/f_src/molval.o: f_src/molval.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/molval.o f_src/molval.F90

${OBJECTDIR}/f_src/molval_I.o: f_src/molval_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/molval_I.o f_src/molval_I.F90

${OBJECTDIR}/f_src/mopend.o: f_src/mopend.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mopend.o f_src/mopend.F90

${OBJECTDIR}/f_src/mopend_I.o: f_src/mopend_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mopend_I.o f_src/mopend_I.F90

${OBJECTDIR}/f_src/mpcbds_I.o: f_src/mpcbds_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mpcbds_I.o f_src/mpcbds_I.F90

${OBJECTDIR}/f_src/mpcpop_I.o: f_src/mpcpop_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mpcpop_I.o f_src/mpcpop_I.F90

${OBJECTDIR}/f_src/mpcsyb.o: f_src/mpcsyb.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mpcsyb.o f_src/mpcsyb.F90

${OBJECTDIR}/f_src/mtxm.o: f_src/mtxm.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mtxm.o f_src/mtxm.F90

${OBJECTDIR}/f_src/mtxm_I.o: f_src/mtxm_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mtxm_I.o f_src/mtxm_I.F90

${OBJECTDIR}/f_src/mtxmc.o: f_src/mtxmc.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mtxmc.o f_src/mtxmc.F90

${OBJECTDIR}/f_src/mtxmc_I.o: f_src/mtxmc_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mtxmc_I.o f_src/mtxmc_I.F90

${OBJECTDIR}/f_src/mullik.o: f_src/mullik.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mullik.o f_src/mullik.F90

${OBJECTDIR}/f_src/mullik_I.o: f_src/mullik_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mullik_I.o f_src/mullik_I.F90

${OBJECTDIR}/f_src/mult.o: f_src/mult.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mult.o f_src/mult.F90

${OBJECTDIR}/f_src/mult33.o: f_src/mult33.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mult33.o f_src/mult33.F90

${OBJECTDIR}/f_src/mult33_I.o: f_src/mult33_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mult33_I.o f_src/mult33_I.F90

${OBJECTDIR}/f_src/mult_I.o: f_src/mult_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mult_I.o f_src/mult_I.F90

${OBJECTDIR}/f_src/mxm.o: f_src/mxm.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mxm.o f_src/mxm.F90

${OBJECTDIR}/f_src/mxm_I.o: f_src/mxm_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mxm_I.o f_src/mxm_I.F90

${OBJECTDIR}/f_src/mxmt.o: f_src/mxmt.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mxmt.o f_src/mxmt.F90

${OBJECTDIR}/f_src/mxmt_I.o: f_src/mxmt_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mxmt_I.o f_src/mxmt_I.F90

${OBJECTDIR}/f_src/mxv.o: f_src/mxv.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mxv.o f_src/mxv.F90

${OBJECTDIR}/f_src/mxv_I.o: f_src/mxv_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/mxv_I.o f_src/mxv_I.F90

${OBJECTDIR}/f_src/myword.o: f_src/myword.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/myword.o f_src/myword.F90

${OBJECTDIR}/f_src/myword_I.o: f_src/myword_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/myword_I.o f_src/myword_I.F90

${OBJECTDIR}/f_src/naican_I.o: f_src/naican_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/naican_I.o f_src/naican_I.F90

${OBJECTDIR}/f_src/naicas_I.o: f_src/naicas_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/naicas_I.o f_src/naicas_I.F90

${OBJECTDIR}/f_src/names.o: f_src/names.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/names.o f_src/names.F90

${OBJECTDIR}/f_src/new_esp.o: f_src/new_esp.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/new_esp.o f_src/new_esp.F90

${OBJECTDIR}/f_src/newflg.o: f_src/newflg.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/newflg.o f_src/newflg.F90

${OBJECTDIR}/f_src/ngamtg_I.o: f_src/ngamtg_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ngamtg_I.o f_src/ngamtg_I.F90

${OBJECTDIR}/f_src/ngefis_I.o: f_src/ngefis_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ngefis_I.o f_src/ngefis_I.F90

${OBJECTDIR}/f_src/ngidri_I.o: f_src/ngidri_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ngidri_I.o f_src/ngidri_I.F90

${OBJECTDIR}/f_src/ngoke_I.o: f_src/ngoke_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ngoke_I.o f_src/ngoke_I.F90

${OBJECTDIR}/f_src/nllsn_I.o: f_src/nllsn_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/nllsn_I.o f_src/nllsn_I.F90

${OBJECTDIR}/f_src/nonbet_I.o: f_src/nonbet_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/nonbet_I.o f_src/nonbet_I.F90

${OBJECTDIR}/f_src/nonope_I.o: f_src/nonope_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/nonope_I.o f_src/nonope_I.F90

${OBJECTDIR}/f_src/nonor_I.o: f_src/nonor_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/nonor_I.o f_src/nonor_I.F90

${OBJECTDIR}/f_src/nuchar.o: f_src/nuchar.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/nuchar.o f_src/nuchar.F90

${OBJECTDIR}/f_src/nuchar_I.o: f_src/nuchar_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/nuchar_I.o f_src/nuchar_I.F90

${OBJECTDIR}/f_src/nxtmer.o: f_src/nxtmer.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/nxtmer.o f_src/nxtmer.F90

${OBJECTDIR}/f_src/openda_I.o: f_src/openda_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/openda_I.o f_src/openda_I.F90

${OBJECTDIR}/f_src/orient_I.o: f_src/orient_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/orient_I.o f_src/orient_I.F90

${OBJECTDIR}/f_src/osinv.o: f_src/osinv.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/osinv.o f_src/osinv.F90

${OBJECTDIR}/f_src/osinv_I.o: f_src/osinv_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/osinv_I.o f_src/osinv_I.F90

${OBJECTDIR}/f_src/outer1.o: f_src/outer1.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/outer1.o f_src/outer1.F90

${OBJECTDIR}/f_src/outer2.o: f_src/outer2.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/outer2.o f_src/outer2.F90

${OBJECTDIR}/f_src/overlaps_C.o: f_src/overlaps_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/overlaps_C.o f_src/overlaps_C.F90

${OBJECTDIR}/f_src/ovlp_I.o: f_src/ovlp_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ovlp_I.o f_src/ovlp_I.F90

${OBJECTDIR}/f_src/packp_I.o: f_src/packp_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/packp_I.o f_src/packp_I.F90

${OBJECTDIR}/f_src/parameters_C.o: f_src/parameters_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/parameters_C.o f_src/parameters_C.F90

${OBJECTDIR}/f_src/parameters_for_PM7_C.o: f_src/parameters_for_PM7_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/parameters_for_PM7_C.o f_src/parameters_for_PM7_C.F90

${OBJECTDIR}/f_src/parameters_for_PM7_Sparkles_C.o: f_src/parameters_for_PM7_Sparkles_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/parameters_for_PM7_Sparkles_C.o f_src/parameters_for_PM7_Sparkles_C.F90

${OBJECTDIR}/f_src/parameters_for_PM7_TS_C.o: f_src/parameters_for_PM7_TS_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/parameters_for_PM7_TS_C.o f_src/parameters_for_PM7_TS_C.F90

${OBJECTDIR}/f_src/partxy.o: f_src/partxy.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/partxy.o f_src/partxy.F90

${OBJECTDIR}/f_src/partxy_I.o: f_src/partxy_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/partxy_I.o f_src/partxy_I.F90

${OBJECTDIR}/f_src/pdgrid_I.o: f_src/pdgrid_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/pdgrid_I.o f_src/pdgrid_I.F90

${OBJECTDIR}/f_src/perm.o: f_src/perm.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/perm.o f_src/perm.F90

${OBJECTDIR}/f_src/perm_I.o: f_src/perm_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/perm_I.o f_src/perm_I.F90

${OBJECTDIR}/f_src/picopt.o: f_src/picopt.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/picopt.o f_src/picopt.F90

${OBJECTDIR}/f_src/pinout.o: f_src/pinout.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/pinout.o f_src/pinout.F90

${OBJECTDIR}/f_src/plato_I.o: f_src/plato_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/plato_I.o f_src/plato_I.F90

${OBJECTDIR}/f_src/pmep.o: f_src/pmep.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/pmep.o f_src/pmep.F90

${OBJECTDIR}/f_src/pmep1_I.o: f_src/pmep1_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/pmep1_I.o f_src/pmep1_I.F90

${OBJECTDIR}/f_src/pmep_I.o: f_src/pmep_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/pmep_I.o f_src/pmep_I.F90

${OBJECTDIR}/f_src/pmepco_I.o: f_src/pmepco_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/pmepco_I.o f_src/pmepco_I.F90

${OBJECTDIR}/f_src/poij_I.o: f_src/poij_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/poij_I.o f_src/poij_I.F90

${OBJECTDIR}/f_src/pol_vol_I.o: f_src/pol_vol_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/pol_vol_I.o f_src/pol_vol_I.F90

${OBJECTDIR}/f_src/post_scf_corrections.o: f_src/post_scf_corrections.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/post_scf_corrections.o f_src/post_scf_corrections.F90

${OBJECTDIR}/f_src/potcal_I.o: f_src/potcal_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/potcal_I.o f_src/potcal_I.F90

${OBJECTDIR}/f_src/powsav_I.o: f_src/powsav_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/powsav_I.o f_src/powsav_I.F90

${OBJECTDIR}/f_src/printp_I.o: f_src/printp_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/printp_I.o f_src/printp_I.F90

${OBJECTDIR}/f_src/prtdrc.o: f_src/prtdrc.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/prtdrc.o f_src/prtdrc.F90

${OBJECTDIR}/f_src/prtdrc_I.o: f_src/prtdrc_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/prtdrc_I.o f_src/prtdrc_I.F90

${OBJECTDIR}/f_src/prtgra.o: f_src/prtgra.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/prtgra.o f_src/prtgra.F90

${OBJECTDIR}/f_src/prthco_I.o: f_src/prthco_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/prthco_I.o f_src/prthco_I.F90

${OBJECTDIR}/f_src/prthes_I.o: f_src/prthes_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/prthes_I.o f_src/prthes_I.F90

${OBJECTDIR}/f_src/prtlmo.o: f_src/prtlmo.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/prtlmo.o f_src/prtlmo.F90

${OBJECTDIR}/f_src/prtpar_I.o: f_src/prtpar_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/prtpar_I.o f_src/prtpar_I.F90

${OBJECTDIR}/f_src/prttim.o: f_src/prttim.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/prttim.o f_src/prttim.F90

${OBJECTDIR}/f_src/prttim_I.o: f_src/prttim_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/prttim_I.o f_src/prttim_I.F90

${OBJECTDIR}/f_src/pulay.o: f_src/pulay.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/pulay.o f_src/pulay.F90

${OBJECTDIR}/f_src/pulay_I.o: f_src/pulay_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/pulay_I.o f_src/pulay_I.F90

${OBJECTDIR}/f_src/quadr.o: f_src/quadr.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/quadr.o f_src/quadr.F90

${OBJECTDIR}/f_src/quadr_I.o: f_src/quadr_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/quadr_I.o f_src/quadr_I.F90

${OBJECTDIR}/f_src/react1.o: f_src/react1.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/react1.o f_src/react1.F90

${OBJECTDIR}/f_src/react1_I.o: f_src/react1_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/react1_I.o f_src/react1_I.F90

${OBJECTDIR}/f_src/reada.o: f_src/reada.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/reada.o f_src/reada.F90

${OBJECTDIR}/f_src/reada_I.o: f_src/reada_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/reada_I.o f_src/reada_I.F90

${OBJECTDIR}/f_src/redatm_I.o: f_src/redatm_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/redatm_I.o f_src/redatm_I.F90

${OBJECTDIR}/f_src/refer.o: f_src/refer.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/refer.o f_src/refer.F90

${OBJECTDIR}/f_src/refer_I.o: f_src/refer_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/refer_I.o f_src/refer_I.F90

${OBJECTDIR}/f_src/refkey_C.o: f_src/refkey_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/refkey_C.o f_src/refkey_C.F90

${OBJECTDIR}/f_src/reorth.o: f_src/reorth.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/reorth.o f_src/reorth.F90

${OBJECTDIR}/f_src/repp_I.o: f_src/repp_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/repp_I.o f_src/repp_I.F90

${OBJECTDIR}/f_src/reppd2_I.o: f_src/reppd2_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/reppd2_I.o f_src/reppd2_I.F90

${OBJECTDIR}/f_src/reppd_I.o: f_src/reppd_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/reppd_I.o f_src/reppd_I.F90

${OBJECTDIR}/f_src/reseq.o: f_src/reseq.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/reseq.o f_src/reseq.F90

${OBJECTDIR}/f_src/resolv.o: f_src/resolv.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/resolv.o f_src/resolv.F90

${OBJECTDIR}/f_src/resolv_I.o: f_src/resolv_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/resolv_I.o f_src/resolv_I.F90

${OBJECTDIR}/f_src/rijkl_I.o: f_src/rijkl_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/rijkl_I.o f_src/rijkl_I.F90

${OBJECTDIR}/f_src/rotat_I.o: f_src/rotat_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/rotat_I.o f_src/rotat_I.F90

${OBJECTDIR}/f_src/rotatd_I.o: f_src/rotatd_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/rotatd_I.o f_src/rotatd_I.F90

${OBJECTDIR}/f_src/rotate.o: f_src/rotate.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/rotate.o f_src/rotate.F90

${OBJECTDIR}/f_src/rotate_C.o: f_src/rotate_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/rotate_C.o f_src/rotate_C.F90

${OBJECTDIR}/f_src/rotate_I.o: f_src/rotate_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/rotate_I.o f_src/rotate_I.F90

${OBJECTDIR}/f_src/rotlmo.o: f_src/rotlmo.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/rotlmo.o f_src/rotlmo.F90

${OBJECTDIR}/f_src/rotmat_I.o: f_src/rotmat_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/rotmat_I.o f_src/rotmat_I.F90

${OBJECTDIR}/f_src/rotmol.o: f_src/rotmol.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/rotmol.o f_src/rotmol.F90

${OBJECTDIR}/f_src/rotmol_I.o: f_src/rotmol_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/rotmol_I.o f_src/rotmol_I.F90

${OBJECTDIR}/f_src/rsc_I.o: f_src/rsc_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/rsc_I.o f_src/rsc_I.F90

${OBJECTDIR}/f_src/rsp.o: f_src/rsp.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/rsp.o f_src/rsp.F90

${OBJECTDIR}/f_src/scfcri.o: f_src/scfcri.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/scfcri.o f_src/scfcri.F90

${OBJECTDIR}/f_src/schmib.o: f_src/schmib.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/schmib.o f_src/schmib.F90

${OBJECTDIR}/f_src/schmib_I.o: f_src/schmib_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/schmib_I.o f_src/schmib_I.F90

${OBJECTDIR}/f_src/schmit.o: f_src/schmit.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/schmit.o f_src/schmit.F90

${OBJECTDIR}/f_src/schmit_I.o: f_src/schmit_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/schmit_I.o f_src/schmit_I.F90

${OBJECTDIR}/f_src/scprm_I.o: f_src/scprm_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/scprm_I.o f_src/scprm_I.F90

${OBJECTDIR}/f_src/search_I.o: f_src/search_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/search_I.o f_src/search_I.F90

${OBJECTDIR}/f_src/second.o: f_src/second.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/second.o f_src/second.F90

${OBJECTDIR}/f_src/second_I.o: f_src/second_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/second_I.o f_src/second_I.F90

${OBJECTDIR}/f_src/selmos.o: f_src/selmos.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/selmos.o f_src/selmos.F90

${OBJECTDIR}/f_src/set.o: f_src/set.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/set.o f_src/set.F90

${OBJECTDIR}/f_src/set_I.o: f_src/set_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/set_I.o f_src/set_I.F90

${OBJECTDIR}/f_src/set_up_MOZYME_arrays.o: f_src/set_up_MOZYME_arrays.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/set_up_MOZYME_arrays.o f_src/set_up_MOZYME_arrays.F90

${OBJECTDIR}/f_src/set_up_RAPID.o: f_src/set_up_RAPID.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/set_up_RAPID.o f_src/set_up_RAPID.F90

${OBJECTDIR}/f_src/set_up_dentate.o: f_src/set_up_dentate.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/set_up_dentate.o f_src/set_up_dentate.F90

${OBJECTDIR}/f_src/setup3_I.o: f_src/setup3_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/setup3_I.o f_src/setup3_I.F90

${OBJECTDIR}/f_src/setup_mopac_arrays.o: f_src/setup_mopac_arrays.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/setup_mopac_arrays.o f_src/setup_mopac_arrays.F90

${OBJECTDIR}/f_src/setupg.o: f_src/setupg.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/setupg.o f_src/setupg.F90

${OBJECTDIR}/f_src/setupk.o: f_src/setupk.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/setupk.o f_src/setupk.F90

${OBJECTDIR}/f_src/solrot.o: f_src/solrot.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/solrot.o f_src/solrot.F90

${OBJECTDIR}/f_src/solrot_I.o: f_src/solrot_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/solrot_I.o f_src/solrot_I.F90

${OBJECTDIR}/f_src/sort.o: f_src/sort.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/sort.o f_src/sort.F90

${OBJECTDIR}/f_src/sort_I.o: f_src/sort_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/sort_I.o f_src/sort_I.F90

${OBJECTDIR}/f_src/sp_two_electron.o: f_src/sp_two_electron.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/sp_two_electron.o f_src/sp_two_electron.F90

${OBJECTDIR}/f_src/spcore_I.o: f_src/spcore_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/spcore_I.o f_src/spcore_I.F90

${OBJECTDIR}/f_src/spline_I.o: f_src/spline_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/spline_I.o f_src/spline_I.F90

${OBJECTDIR}/f_src/ss_I.o: f_src/ss_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/ss_I.o f_src/ss_I.F90

${OBJECTDIR}/f_src/suma2_I.o: f_src/suma2_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/suma2_I.o f_src/suma2_I.F90

${OBJECTDIR}/f_src/supdot.o: f_src/supdot.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/supdot.o f_src/supdot.F90

${OBJECTDIR}/f_src/supdot_I.o: f_src/supdot_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/supdot_I.o f_src/supdot_I.F90

${OBJECTDIR}/f_src/superd.o: f_src/superd.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/superd.o f_src/superd.F90

${OBJECTDIR}/f_src/superd_I.o: f_src/superd_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/superd_I.o f_src/superd_I.F90

${OBJECTDIR}/f_src/surfa_I.o: f_src/surfa_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/surfa_I.o f_src/surfa_I.F90

${OBJECTDIR}/f_src/surfac_I.o: f_src/surfac_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/surfac_I.o f_src/surfac_I.F90

${OBJECTDIR}/f_src/swap.o: f_src/swap.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/swap.o f_src/swap.F90

${OBJECTDIR}/f_src/swap_I.o: f_src/swap_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/swap_I.o f_src/swap_I.F90

${OBJECTDIR}/f_src/switch.o: f_src/switch.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/switch.o f_src/switch.F90

${OBJECTDIR}/f_src/symdec_I.o: f_src/symdec_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/symdec_I.o f_src/symdec_I.F90

${OBJECTDIR}/f_src/symh.o: f_src/symh.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/symh.o f_src/symh.F90

${OBJECTDIR}/f_src/symh_I.o: f_src/symh_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/symh_I.o f_src/symh_I.F90

${OBJECTDIR}/f_src/symmetry_C.o: f_src/symmetry_C.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/symmetry_C.o f_src/symmetry_C.F90

${OBJECTDIR}/f_src/symopr.o: f_src/symopr.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/symopr.o f_src/symopr.F90

${OBJECTDIR}/f_src/symopr_I.o: f_src/symopr_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/symopr_I.o f_src/symopr_I.F90

${OBJECTDIR}/f_src/symp_I.o: f_src/symp_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/symp_I.o f_src/symp_I.F90

${OBJECTDIR}/f_src/sympop.o: f_src/sympop.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/sympop.o f_src/sympop.F90

${OBJECTDIR}/f_src/sympop_I.o: f_src/sympop_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/sympop_I.o f_src/sympop_I.F90

${OBJECTDIR}/f_src/symr.o: f_src/symr.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/symr.o f_src/symr.F90

${OBJECTDIR}/f_src/symr_I.o: f_src/symr_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/symr_I.o f_src/symr_I.F90

${OBJECTDIR}/f_src/symt.o: f_src/symt.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/symt.o f_src/symt.F90

${OBJECTDIR}/f_src/symt_I.o: f_src/symt_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/symt_I.o f_src/symt_I.F90

${OBJECTDIR}/f_src/symtry.o: f_src/symtry.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/symtry.o f_src/symtry.F90

${OBJECTDIR}/f_src/symtry_I.o: f_src/symtry_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/symtry_I.o f_src/symtry_I.F90

${OBJECTDIR}/f_src/tf_I.o: f_src/tf_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/tf_I.o f_src/tf_I.F90

${OBJECTDIR}/f_src/tidy.o: f_src/tidy.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/tidy.o f_src/tidy.F90

${OBJECTDIR}/f_src/time_I.o: f_src/time_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/time_I.o f_src/time_I.F90

${OBJECTDIR}/f_src/timer.o: f_src/timer.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/timer.o f_src/timer.F90

${OBJECTDIR}/f_src/timer_I.o: f_src/timer_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/timer_I.o f_src/timer_I.F90

${OBJECTDIR}/f_src/timout.o: f_src/timout.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/timout.o f_src/timout.F90

${OBJECTDIR}/f_src/transf_I.o: f_src/transf_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/transf_I.o f_src/transf_I.F90

${OBJECTDIR}/f_src/trsub_I.o: f_src/trsub_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/trsub_I.o f_src/trsub_I.F90

${OBJECTDIR}/f_src/trudgu_I.o: f_src/trudgu_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/trudgu_I.o f_src/trudgu_I.F90

${OBJECTDIR}/f_src/trugdu_I.o: f_src/trugdu_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/trugdu_I.o f_src/trugdu_I.F90

${OBJECTDIR}/f_src/trugud_I.o: f_src/trugud_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/trugud_I.o f_src/trugud_I.F90

${OBJECTDIR}/f_src/tx_I.o: f_src/tx_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/tx_I.o f_src/tx_I.F90

${OBJECTDIR}/f_src/txtype.o: f_src/txtype.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/txtype.o f_src/txtype.F90

${OBJECTDIR}/f_src/upcase.o: f_src/upcase.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/upcase.o f_src/upcase.F90

${OBJECTDIR}/f_src/upcase_I.o: f_src/upcase_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/upcase_I.o f_src/upcase_I.F90

${OBJECTDIR}/f_src/update.o: f_src/update.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/update.o f_src/update.F90

${OBJECTDIR}/f_src/update_I.o: f_src/update_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/update_I.o f_src/update_I.F90

${OBJECTDIR}/f_src/values.o: f_src/values.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/values.o f_src/values.F90

${OBJECTDIR}/f_src/vastkind.o: f_src/vastkind.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/vastkind.o f_src/vastkind.F90

${OBJECTDIR}/f_src/vecprt.o: f_src/vecprt.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/vecprt.o f_src/vecprt.F90

${OBJECTDIR}/f_src/vecprt_I.o: f_src/vecprt_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/vecprt_I.o f_src/vecprt_I.F90

${OBJECTDIR}/f_src/vecprt_for_MOZYME.o: f_src/vecprt_for_MOZYME.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/vecprt_for_MOZYME.o f_src/vecprt_for_MOZYME.F90

${OBJECTDIR}/f_src/volume.o: f_src/volume.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/volume.o f_src/volume.F90

${OBJECTDIR}/f_src/volume_I.o: f_src/volume_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/volume_I.o f_src/volume_I.F90

${OBJECTDIR}/f_src/w2mat_I.o: f_src/w2mat_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/w2mat_I.o f_src/w2mat_I.F90

${OBJECTDIR}/f_src/worder_I.o: f_src/worder_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/worder_I.o f_src/worder_I.F90

${OBJECTDIR}/f_src/wrdkey_I.o: f_src/wrdkey_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/wrdkey_I.o f_src/wrdkey_I.F90

${OBJECTDIR}/f_src/wrtkey.o: f_src/wrtkey.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/wrtkey.o f_src/wrtkey.F90

${OBJECTDIR}/f_src/wrtkey_I.o: f_src/wrtkey_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/wrtkey_I.o f_src/wrtkey_I.F90

${OBJECTDIR}/f_src/wrttxt.o: f_src/wrttxt.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/wrttxt.o f_src/wrttxt.F90

${OBJECTDIR}/f_src/wrttxt_I.o: f_src/wrttxt_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/wrttxt_I.o f_src/wrttxt_I.F90

${OBJECTDIR}/f_src/wstore_I.o: f_src/wstore_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/wstore_I.o f_src/wstore_I.F90

${OBJECTDIR}/f_src/xyzcry.o: f_src/xyzcry.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/xyzcry.o f_src/xyzcry.F90

${OBJECTDIR}/f_src/xyzcry_I.o: f_src/xyzcry_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/xyzcry_I.o f_src/xyzcry_I.F90

${OBJECTDIR}/f_src/xyzint.o: f_src/xyzint.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/xyzint.o f_src/xyzint.F90

${OBJECTDIR}/f_src/xyzint_I.o: f_src/xyzint_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/xyzint_I.o f_src/xyzint_I.F90

${OBJECTDIR}/f_src/zerom_I.o: f_src/zerom_I.F90
	${MKDIR} -p ${OBJECTDIR}/f_src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/f_src/zerom_I.o f_src/zerom_I.F90

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} *.mod

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
