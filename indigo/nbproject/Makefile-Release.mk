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
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
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
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_abbreviations_core.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_abbreviations_expand.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_array.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_basic.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_calc.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_composition.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_deconvolution.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_fingerprints.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_io.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_layout.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_loaders.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_macros.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_mapping.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_match.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_misc.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_molecule.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_object.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_options.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_product_enumerator.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_properties.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_reaction.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_savers.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_scaffold.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_stereo.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_structure_checker.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_tautomer_enumerator.o \
	${OBJECTDIR}/Indigo/api/c/indigo/src/option_manager.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_c/bitarray.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_c/nano_posix.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_c/os_dir_posix.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/bitinworker.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/bitoutworker.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/cancellation_handler.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/chunk_storage.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/common_exceptions_impl.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/crc32.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/csv_reader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/d_bitset.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/exception.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/gray_codes.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/io_base.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/locale_guard.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/os_sync_wrapper.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/os_thread_wrapper.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/output.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/profiling.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/properties_map.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/scanner.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/shmem_posix.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/shmem_win32.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/smart_output.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/string_pool.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/tlscont.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/gzip/gzip_output.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/gzip/gzip_scanner.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/lzw/lzw_decoder.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/lzw/lzw_dictionary.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/lzw/lzw_encoder.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/math/best_fit.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/math/line3f.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/math/lseg3f.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/math/matr3x3d.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/math/plane3f.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/math/random.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/math/statistics.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/math/transform3f.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/math/vec2f.o \
	${OBJECTDIR}/Indigo/core/indigo-core/common/math/vec3f.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/automorphism_search.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/aux_path_finder.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/biconnected_decomposer.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/cycle_basis.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/cycle_enumerator.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/dfs_walk.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/edge_rotation_matcher.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/edge_subgraph_enumerator.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/embedding_enumerator.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/embeddings_storage.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/filter.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_affine_matcher.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_constrained_bmatching_finder.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_decomposer.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_fast_access.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_iterators.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_perfect_matching.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_subchain_enumerator.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_subtree_enumerator.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/max_common_subgraph.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/morgan_code.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/path_enumerator.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/scaffold_detection.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/shortest_path_finder.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/simple_cycle_basis.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/skew_symmetric_flow_finder.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/skew_symmetric_network.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/spanning_tree.o \
	${OBJECTDIR}/Indigo/core/indigo-core/graph/src/subgraph_hash.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/attachment_layout.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/layout_pattern.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/layout_pattern_holder.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/layout_pattern_smart.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/metalayout.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_cleaner_2d.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_assign.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_assign_smart.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_attach.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_border.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_cycle.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_geom.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_refine.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_smart.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_macrocycle_lattice.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_macrocycles.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/reaction_layout.o \
	${OBJECTDIR}/Indigo/core/indigo-core/layout/src/refinement_state.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/base_molecule.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/canonical_smiles_saver.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/cmf_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/cmf_saver.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/cml_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/cml_saver.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/crippen.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/elements.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/haworth_projection_finder.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/hybridization.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/icm_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/icm_saver.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/inchi_parser.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/inchi_wrapper.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/lipinski.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/max_common_submolecule.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/metadata_storage.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_3d_constraints.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_allene_stereo.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_arom.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_arom_match.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_auto_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_automorphism_search.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_cdxml_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_cdxml_saver.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_chain_fingerprints.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_cip_calculator.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_cis_trans.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_dearom.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_electrons_localizer.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_exact_matcher.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_exact_substructure_matcher.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_fingerprint.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_gross_formula.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_gross_formula_options.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_hash.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_inchi.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_inchi_component.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_inchi_layers.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_inchi_utils.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_ionize.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_json_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_json_saver.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_layered_molecules.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_mass.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_mass_options.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_morgan_fingerprint_builder.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_name_parser.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_neighbourhood_counters.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_pi_systems_matcher.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_rgroups.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_rgroups_composition.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_savers.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_scaffold_detection.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_sgroups.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_standardize.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_standardize_options.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_stereocenter_options.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_stereocenters.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_substructure_matcher.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_chain.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_enumerator.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_match.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_matcher.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_substructure_matcher.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_superstructure.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_utils.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tgroups.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molfile_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molfile_saver.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/monomer_commons.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/monomers_lib.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/multiple_cdx_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/multiple_cml_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/parse_utils.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/query_molecule.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/rdf_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/sdf_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/smiles_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/smiles_saver.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/structure_checker.o \
	${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/tpsa.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/base_reaction.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/base_reaction_substructure_matcher.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/canonical_rsmiles_saver.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/crf_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/crf_saver.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/icr_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/icr_saver.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/query_reaction.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_auto_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_automapper.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_cdxml_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_cdxml_saver.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_cml_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_cml_saver.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_enumerator_state.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_exact_matcher.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_fingerprint.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_gross_formula.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_hash.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_json_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_json_saver.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_neighborhood_counters.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_product_enumerator.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_substructure_matcher.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_transformation.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/rsmiles_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/rsmiles_saver.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/rxnfile_loader.o \
	${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/rxnfile_saver.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src/ichilnct.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll_a.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll_a2.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll_b.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichi_bns.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichi_io.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichican2.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichicano.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichicans.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichierr.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiisot.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichimak2.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichimake.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichimap1.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichimap2.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichimap4.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichinorm.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiparm.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiprt1.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiprt2.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiprt3.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiqueu.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiread.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiring.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr1.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr2.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr3.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr4.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr5.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr6.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr7.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichisort.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichister.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichitaut.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ikey_base26.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ikey_dll.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/inchi_gui.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/mol2atom.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt1.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt2.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt3.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt4.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/readinch.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/runichi.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/runichi2.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/runichi3.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/runichi4.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/sha2.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/strutil.o \
	${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/util.o \
	${OBJECTDIR}/Indigo/third_party/tinyxml2/tinyxml2.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-std=c++17
CXXFLAGS=-std=c++17

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libindigo.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libindigo.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libindigo.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libindigo.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libindigo.a

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo.o: Indigo/api/c/indigo/src/indigo.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo.o Indigo/api/c/indigo/src/indigo.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_abbreviations_core.o: Indigo/api/c/indigo/src/indigo_abbreviations_core.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_abbreviations_core.o Indigo/api/c/indigo/src/indigo_abbreviations_core.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_abbreviations_expand.o: Indigo/api/c/indigo/src/indigo_abbreviations_expand.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_abbreviations_expand.o Indigo/api/c/indigo/src/indigo_abbreviations_expand.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_array.o: Indigo/api/c/indigo/src/indigo_array.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_array.o Indigo/api/c/indigo/src/indigo_array.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_basic.o: Indigo/api/c/indigo/src/indigo_basic.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_basic.o Indigo/api/c/indigo/src/indigo_basic.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_calc.o: Indigo/api/c/indigo/src/indigo_calc.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_calc.o Indigo/api/c/indigo/src/indigo_calc.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_composition.o: Indigo/api/c/indigo/src/indigo_composition.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_composition.o Indigo/api/c/indigo/src/indigo_composition.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_deconvolution.o: Indigo/api/c/indigo/src/indigo_deconvolution.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_deconvolution.o Indigo/api/c/indigo/src/indigo_deconvolution.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_fingerprints.o: Indigo/api/c/indigo/src/indigo_fingerprints.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_fingerprints.o Indigo/api/c/indigo/src/indigo_fingerprints.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_io.o: Indigo/api/c/indigo/src/indigo_io.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_io.o Indigo/api/c/indigo/src/indigo_io.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_layout.o: Indigo/api/c/indigo/src/indigo_layout.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_layout.o Indigo/api/c/indigo/src/indigo_layout.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_loaders.o: Indigo/api/c/indigo/src/indigo_loaders.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_loaders.o Indigo/api/c/indigo/src/indigo_loaders.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_macros.o: Indigo/api/c/indigo/src/indigo_macros.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_macros.o Indigo/api/c/indigo/src/indigo_macros.c

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_mapping.o: Indigo/api/c/indigo/src/indigo_mapping.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_mapping.o Indigo/api/c/indigo/src/indigo_mapping.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_match.o: Indigo/api/c/indigo/src/indigo_match.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_match.o Indigo/api/c/indigo/src/indigo_match.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_misc.o: Indigo/api/c/indigo/src/indigo_misc.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_misc.o Indigo/api/c/indigo/src/indigo_misc.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_molecule.o: Indigo/api/c/indigo/src/indigo_molecule.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_molecule.o Indigo/api/c/indigo/src/indigo_molecule.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_object.o: Indigo/api/c/indigo/src/indigo_object.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_object.o Indigo/api/c/indigo/src/indigo_object.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_options.o: Indigo/api/c/indigo/src/indigo_options.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_options.o Indigo/api/c/indigo/src/indigo_options.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_product_enumerator.o: Indigo/api/c/indigo/src/indigo_product_enumerator.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_product_enumerator.o Indigo/api/c/indigo/src/indigo_product_enumerator.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_properties.o: Indigo/api/c/indigo/src/indigo_properties.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_properties.o Indigo/api/c/indigo/src/indigo_properties.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_reaction.o: Indigo/api/c/indigo/src/indigo_reaction.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_reaction.o Indigo/api/c/indigo/src/indigo_reaction.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_savers.o: Indigo/api/c/indigo/src/indigo_savers.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_savers.o Indigo/api/c/indigo/src/indigo_savers.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_scaffold.o: Indigo/api/c/indigo/src/indigo_scaffold.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_scaffold.o Indigo/api/c/indigo/src/indigo_scaffold.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_stereo.o: Indigo/api/c/indigo/src/indigo_stereo.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_stereo.o Indigo/api/c/indigo/src/indigo_stereo.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_structure_checker.o: Indigo/api/c/indigo/src/indigo_structure_checker.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_structure_checker.o Indigo/api/c/indigo/src/indigo_structure_checker.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_tautomer_enumerator.o: Indigo/api/c/indigo/src/indigo_tautomer_enumerator.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/indigo_tautomer_enumerator.o Indigo/api/c/indigo/src/indigo_tautomer_enumerator.cpp

${OBJECTDIR}/Indigo/api/c/indigo/src/option_manager.o: Indigo/api/c/indigo/src/option_manager.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/api/c/indigo/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/api/c/indigo/src/option_manager.o Indigo/api/c/indigo/src/option_manager.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_c/bitarray.o: Indigo/core/indigo-core/common/base_c/bitarray.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_c
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_c/bitarray.o Indigo/core/indigo-core/common/base_c/bitarray.c

${OBJECTDIR}/Indigo/core/indigo-core/common/base_c/nano_posix.o: Indigo/core/indigo-core/common/base_c/nano_posix.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_c
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_c/nano_posix.o Indigo/core/indigo-core/common/base_c/nano_posix.c

${OBJECTDIR}/Indigo/core/indigo-core/common/base_c/os_dir_posix.o: Indigo/core/indigo-core/common/base_c/os_dir_posix.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_c
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_c/os_dir_posix.o Indigo/core/indigo-core/common/base_c/os_dir_posix.c

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/bitinworker.o: Indigo/core/indigo-core/common/base_cpp/bitinworker.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/bitinworker.o Indigo/core/indigo-core/common/base_cpp/bitinworker.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/bitoutworker.o: Indigo/core/indigo-core/common/base_cpp/bitoutworker.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/bitoutworker.o Indigo/core/indigo-core/common/base_cpp/bitoutworker.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/cancellation_handler.o: Indigo/core/indigo-core/common/base_cpp/cancellation_handler.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/cancellation_handler.o Indigo/core/indigo-core/common/base_cpp/cancellation_handler.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/chunk_storage.o: Indigo/core/indigo-core/common/base_cpp/chunk_storage.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/chunk_storage.o Indigo/core/indigo-core/common/base_cpp/chunk_storage.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/common_exceptions_impl.o: Indigo/core/indigo-core/common/base_cpp/common_exceptions_impl.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/common_exceptions_impl.o Indigo/core/indigo-core/common/base_cpp/common_exceptions_impl.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/crc32.o: Indigo/core/indigo-core/common/base_cpp/crc32.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/crc32.o Indigo/core/indigo-core/common/base_cpp/crc32.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/csv_reader.o: Indigo/core/indigo-core/common/base_cpp/csv_reader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/csv_reader.o Indigo/core/indigo-core/common/base_cpp/csv_reader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/d_bitset.o: Indigo/core/indigo-core/common/base_cpp/d_bitset.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/d_bitset.o Indigo/core/indigo-core/common/base_cpp/d_bitset.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/exception.o: Indigo/core/indigo-core/common/base_cpp/exception.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/exception.o Indigo/core/indigo-core/common/base_cpp/exception.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/gray_codes.o: Indigo/core/indigo-core/common/base_cpp/gray_codes.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/gray_codes.o Indigo/core/indigo-core/common/base_cpp/gray_codes.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/io_base.o: Indigo/core/indigo-core/common/base_cpp/io_base.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/io_base.o Indigo/core/indigo-core/common/base_cpp/io_base.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/locale_guard.o: Indigo/core/indigo-core/common/base_cpp/locale_guard.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/locale_guard.o Indigo/core/indigo-core/common/base_cpp/locale_guard.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/os_sync_wrapper.o: Indigo/core/indigo-core/common/base_cpp/os_sync_wrapper.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/os_sync_wrapper.o Indigo/core/indigo-core/common/base_cpp/os_sync_wrapper.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/os_thread_wrapper.o: Indigo/core/indigo-core/common/base_cpp/os_thread_wrapper.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/os_thread_wrapper.o Indigo/core/indigo-core/common/base_cpp/os_thread_wrapper.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/output.o: Indigo/core/indigo-core/common/base_cpp/output.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/output.o Indigo/core/indigo-core/common/base_cpp/output.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/profiling.o: Indigo/core/indigo-core/common/base_cpp/profiling.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/profiling.o Indigo/core/indigo-core/common/base_cpp/profiling.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/properties_map.o: Indigo/core/indigo-core/common/base_cpp/properties_map.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/properties_map.o Indigo/core/indigo-core/common/base_cpp/properties_map.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/scanner.o: Indigo/core/indigo-core/common/base_cpp/scanner.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/scanner.o Indigo/core/indigo-core/common/base_cpp/scanner.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/shmem_posix.o: Indigo/core/indigo-core/common/base_cpp/shmem_posix.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/shmem_posix.o Indigo/core/indigo-core/common/base_cpp/shmem_posix.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/shmem_win32.o: Indigo/core/indigo-core/common/base_cpp/shmem_win32.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/shmem_win32.o Indigo/core/indigo-core/common/base_cpp/shmem_win32.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/smart_output.o: Indigo/core/indigo-core/common/base_cpp/smart_output.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/smart_output.o Indigo/core/indigo-core/common/base_cpp/smart_output.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/string_pool.o: Indigo/core/indigo-core/common/base_cpp/string_pool.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/string_pool.o Indigo/core/indigo-core/common/base_cpp/string_pool.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/tlscont.o: Indigo/core/indigo-core/common/base_cpp/tlscont.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/base_cpp/tlscont.o Indigo/core/indigo-core/common/base_cpp/tlscont.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/gzip/gzip_output.o: Indigo/core/indigo-core/common/gzip/gzip_output.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/gzip
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/gzip/gzip_output.o Indigo/core/indigo-core/common/gzip/gzip_output.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/gzip/gzip_scanner.o: Indigo/core/indigo-core/common/gzip/gzip_scanner.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/gzip
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/gzip/gzip_scanner.o Indigo/core/indigo-core/common/gzip/gzip_scanner.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/lzw/lzw_decoder.o: Indigo/core/indigo-core/common/lzw/lzw_decoder.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/lzw
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/lzw/lzw_decoder.o Indigo/core/indigo-core/common/lzw/lzw_decoder.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/lzw/lzw_dictionary.o: Indigo/core/indigo-core/common/lzw/lzw_dictionary.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/lzw
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/lzw/lzw_dictionary.o Indigo/core/indigo-core/common/lzw/lzw_dictionary.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/lzw/lzw_encoder.o: Indigo/core/indigo-core/common/lzw/lzw_encoder.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/lzw
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/lzw/lzw_encoder.o Indigo/core/indigo-core/common/lzw/lzw_encoder.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/math/best_fit.o: Indigo/core/indigo-core/common/math/best_fit.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/math
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/math/best_fit.o Indigo/core/indigo-core/common/math/best_fit.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/math/line3f.o: Indigo/core/indigo-core/common/math/line3f.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/math
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/math/line3f.o Indigo/core/indigo-core/common/math/line3f.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/math/lseg3f.o: Indigo/core/indigo-core/common/math/lseg3f.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/math
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/math/lseg3f.o Indigo/core/indigo-core/common/math/lseg3f.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/math/matr3x3d.o: Indigo/core/indigo-core/common/math/matr3x3d.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/math
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/math/matr3x3d.o Indigo/core/indigo-core/common/math/matr3x3d.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/math/plane3f.o: Indigo/core/indigo-core/common/math/plane3f.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/math
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/math/plane3f.o Indigo/core/indigo-core/common/math/plane3f.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/math/random.o: Indigo/core/indigo-core/common/math/random.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/math
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/math/random.o Indigo/core/indigo-core/common/math/random.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/math/statistics.o: Indigo/core/indigo-core/common/math/statistics.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/math
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/math/statistics.o Indigo/core/indigo-core/common/math/statistics.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/math/transform3f.o: Indigo/core/indigo-core/common/math/transform3f.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/math
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/math/transform3f.o Indigo/core/indigo-core/common/math/transform3f.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/math/vec2f.o: Indigo/core/indigo-core/common/math/vec2f.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/math
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/math/vec2f.o Indigo/core/indigo-core/common/math/vec2f.cpp

${OBJECTDIR}/Indigo/core/indigo-core/common/math/vec3f.o: Indigo/core/indigo-core/common/math/vec3f.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/common/math
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/common/math/vec3f.o Indigo/core/indigo-core/common/math/vec3f.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/automorphism_search.o: Indigo/core/indigo-core/graph/src/automorphism_search.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/automorphism_search.o Indigo/core/indigo-core/graph/src/automorphism_search.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/aux_path_finder.o: Indigo/core/indigo-core/graph/src/aux_path_finder.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/aux_path_finder.o Indigo/core/indigo-core/graph/src/aux_path_finder.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/biconnected_decomposer.o: Indigo/core/indigo-core/graph/src/biconnected_decomposer.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/biconnected_decomposer.o Indigo/core/indigo-core/graph/src/biconnected_decomposer.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/cycle_basis.o: Indigo/core/indigo-core/graph/src/cycle_basis.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/cycle_basis.o Indigo/core/indigo-core/graph/src/cycle_basis.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/cycle_enumerator.o: Indigo/core/indigo-core/graph/src/cycle_enumerator.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/cycle_enumerator.o Indigo/core/indigo-core/graph/src/cycle_enumerator.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/dfs_walk.o: Indigo/core/indigo-core/graph/src/dfs_walk.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/dfs_walk.o Indigo/core/indigo-core/graph/src/dfs_walk.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/edge_rotation_matcher.o: Indigo/core/indigo-core/graph/src/edge_rotation_matcher.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/edge_rotation_matcher.o Indigo/core/indigo-core/graph/src/edge_rotation_matcher.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/edge_subgraph_enumerator.o: Indigo/core/indigo-core/graph/src/edge_subgraph_enumerator.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/edge_subgraph_enumerator.o Indigo/core/indigo-core/graph/src/edge_subgraph_enumerator.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/embedding_enumerator.o: Indigo/core/indigo-core/graph/src/embedding_enumerator.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/embedding_enumerator.o Indigo/core/indigo-core/graph/src/embedding_enumerator.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/embeddings_storage.o: Indigo/core/indigo-core/graph/src/embeddings_storage.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/embeddings_storage.o Indigo/core/indigo-core/graph/src/embeddings_storage.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/filter.o: Indigo/core/indigo-core/graph/src/filter.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/filter.o Indigo/core/indigo-core/graph/src/filter.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph.o: Indigo/core/indigo-core/graph/src/graph.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph.o Indigo/core/indigo-core/graph/src/graph.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_affine_matcher.o: Indigo/core/indigo-core/graph/src/graph_affine_matcher.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_affine_matcher.o Indigo/core/indigo-core/graph/src/graph_affine_matcher.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_constrained_bmatching_finder.o: Indigo/core/indigo-core/graph/src/graph_constrained_bmatching_finder.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_constrained_bmatching_finder.o Indigo/core/indigo-core/graph/src/graph_constrained_bmatching_finder.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_decomposer.o: Indigo/core/indigo-core/graph/src/graph_decomposer.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_decomposer.o Indigo/core/indigo-core/graph/src/graph_decomposer.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_fast_access.o: Indigo/core/indigo-core/graph/src/graph_fast_access.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_fast_access.o Indigo/core/indigo-core/graph/src/graph_fast_access.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_iterators.o: Indigo/core/indigo-core/graph/src/graph_iterators.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_iterators.o Indigo/core/indigo-core/graph/src/graph_iterators.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_perfect_matching.o: Indigo/core/indigo-core/graph/src/graph_perfect_matching.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_perfect_matching.o Indigo/core/indigo-core/graph/src/graph_perfect_matching.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_subchain_enumerator.o: Indigo/core/indigo-core/graph/src/graph_subchain_enumerator.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_subchain_enumerator.o Indigo/core/indigo-core/graph/src/graph_subchain_enumerator.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_subtree_enumerator.o: Indigo/core/indigo-core/graph/src/graph_subtree_enumerator.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/graph_subtree_enumerator.o Indigo/core/indigo-core/graph/src/graph_subtree_enumerator.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/max_common_subgraph.o: Indigo/core/indigo-core/graph/src/max_common_subgraph.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/max_common_subgraph.o Indigo/core/indigo-core/graph/src/max_common_subgraph.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/morgan_code.o: Indigo/core/indigo-core/graph/src/morgan_code.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/morgan_code.o Indigo/core/indigo-core/graph/src/morgan_code.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/path_enumerator.o: Indigo/core/indigo-core/graph/src/path_enumerator.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/path_enumerator.o Indigo/core/indigo-core/graph/src/path_enumerator.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/scaffold_detection.o: Indigo/core/indigo-core/graph/src/scaffold_detection.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/scaffold_detection.o Indigo/core/indigo-core/graph/src/scaffold_detection.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/shortest_path_finder.o: Indigo/core/indigo-core/graph/src/shortest_path_finder.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/shortest_path_finder.o Indigo/core/indigo-core/graph/src/shortest_path_finder.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/simple_cycle_basis.o: Indigo/core/indigo-core/graph/src/simple_cycle_basis.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/simple_cycle_basis.o Indigo/core/indigo-core/graph/src/simple_cycle_basis.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/skew_symmetric_flow_finder.o: Indigo/core/indigo-core/graph/src/skew_symmetric_flow_finder.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/skew_symmetric_flow_finder.o Indigo/core/indigo-core/graph/src/skew_symmetric_flow_finder.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/skew_symmetric_network.o: Indigo/core/indigo-core/graph/src/skew_symmetric_network.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/skew_symmetric_network.o Indigo/core/indigo-core/graph/src/skew_symmetric_network.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/spanning_tree.o: Indigo/core/indigo-core/graph/src/spanning_tree.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/spanning_tree.o Indigo/core/indigo-core/graph/src/spanning_tree.cpp

${OBJECTDIR}/Indigo/core/indigo-core/graph/src/subgraph_hash.o: Indigo/core/indigo-core/graph/src/subgraph_hash.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/graph/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/graph/src/subgraph_hash.o Indigo/core/indigo-core/graph/src/subgraph_hash.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/attachment_layout.o: Indigo/core/indigo-core/layout/src/attachment_layout.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/attachment_layout.o Indigo/core/indigo-core/layout/src/attachment_layout.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/layout_pattern.o: Indigo/core/indigo-core/layout/src/layout_pattern.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/layout_pattern.o Indigo/core/indigo-core/layout/src/layout_pattern.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/layout_pattern_holder.o: Indigo/core/indigo-core/layout/src/layout_pattern_holder.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/layout_pattern_holder.o Indigo/core/indigo-core/layout/src/layout_pattern_holder.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/layout_pattern_smart.o: Indigo/core/indigo-core/layout/src/layout_pattern_smart.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/layout_pattern_smart.o Indigo/core/indigo-core/layout/src/layout_pattern_smart.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/metalayout.o: Indigo/core/indigo-core/layout/src/metalayout.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/metalayout.o Indigo/core/indigo-core/layout/src/metalayout.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_cleaner_2d.o: Indigo/core/indigo-core/layout/src/molecule_cleaner_2d.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_cleaner_2d.o Indigo/core/indigo-core/layout/src/molecule_cleaner_2d.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout.o: Indigo/core/indigo-core/layout/src/molecule_layout.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout.o Indigo/core/indigo-core/layout/src/molecule_layout.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph.o: Indigo/core/indigo-core/layout/src/molecule_layout_graph.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph.o Indigo/core/indigo-core/layout/src/molecule_layout_graph.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_assign.o: Indigo/core/indigo-core/layout/src/molecule_layout_graph_assign.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_assign.o Indigo/core/indigo-core/layout/src/molecule_layout_graph_assign.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_assign_smart.o: Indigo/core/indigo-core/layout/src/molecule_layout_graph_assign_smart.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_assign_smart.o Indigo/core/indigo-core/layout/src/molecule_layout_graph_assign_smart.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_attach.o: Indigo/core/indigo-core/layout/src/molecule_layout_graph_attach.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_attach.o Indigo/core/indigo-core/layout/src/molecule_layout_graph_attach.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_border.o: Indigo/core/indigo-core/layout/src/molecule_layout_graph_border.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_border.o Indigo/core/indigo-core/layout/src/molecule_layout_graph_border.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_cycle.o: Indigo/core/indigo-core/layout/src/molecule_layout_graph_cycle.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_cycle.o Indigo/core/indigo-core/layout/src/molecule_layout_graph_cycle.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_geom.o: Indigo/core/indigo-core/layout/src/molecule_layout_graph_geom.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_geom.o Indigo/core/indigo-core/layout/src/molecule_layout_graph_geom.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_refine.o: Indigo/core/indigo-core/layout/src/molecule_layout_graph_refine.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_refine.o Indigo/core/indigo-core/layout/src/molecule_layout_graph_refine.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_smart.o: Indigo/core/indigo-core/layout/src/molecule_layout_graph_smart.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_graph_smart.o Indigo/core/indigo-core/layout/src/molecule_layout_graph_smart.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_macrocycle_lattice.o: Indigo/core/indigo-core/layout/src/molecule_layout_macrocycle_lattice.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_macrocycle_lattice.o Indigo/core/indigo-core/layout/src/molecule_layout_macrocycle_lattice.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_macrocycles.o: Indigo/core/indigo-core/layout/src/molecule_layout_macrocycles.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/molecule_layout_macrocycles.o Indigo/core/indigo-core/layout/src/molecule_layout_macrocycles.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/reaction_layout.o: Indigo/core/indigo-core/layout/src/reaction_layout.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/reaction_layout.o Indigo/core/indigo-core/layout/src/reaction_layout.cpp

${OBJECTDIR}/Indigo/core/indigo-core/layout/src/refinement_state.o: Indigo/core/indigo-core/layout/src/refinement_state.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/layout/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/layout/src/refinement_state.o Indigo/core/indigo-core/layout/src/refinement_state.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/base_molecule.o: Indigo/core/indigo-core/molecule/src/base_molecule.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/base_molecule.o Indigo/core/indigo-core/molecule/src/base_molecule.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/canonical_smiles_saver.o: Indigo/core/indigo-core/molecule/src/canonical_smiles_saver.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/canonical_smiles_saver.o Indigo/core/indigo-core/molecule/src/canonical_smiles_saver.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/cmf_loader.o: Indigo/core/indigo-core/molecule/src/cmf_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/cmf_loader.o Indigo/core/indigo-core/molecule/src/cmf_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/cmf_saver.o: Indigo/core/indigo-core/molecule/src/cmf_saver.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/cmf_saver.o Indigo/core/indigo-core/molecule/src/cmf_saver.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/cml_loader.o: Indigo/core/indigo-core/molecule/src/cml_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/cml_loader.o Indigo/core/indigo-core/molecule/src/cml_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/cml_saver.o: Indigo/core/indigo-core/molecule/src/cml_saver.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/cml_saver.o Indigo/core/indigo-core/molecule/src/cml_saver.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/crippen.o: Indigo/core/indigo-core/molecule/src/crippen.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/crippen.o Indigo/core/indigo-core/molecule/src/crippen.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/elements.o: Indigo/core/indigo-core/molecule/src/elements.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/elements.o Indigo/core/indigo-core/molecule/src/elements.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/haworth_projection_finder.o: Indigo/core/indigo-core/molecule/src/haworth_projection_finder.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/haworth_projection_finder.o Indigo/core/indigo-core/molecule/src/haworth_projection_finder.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/hybridization.o: Indigo/core/indigo-core/molecule/src/hybridization.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/hybridization.o Indigo/core/indigo-core/molecule/src/hybridization.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/icm_loader.o: Indigo/core/indigo-core/molecule/src/icm_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/icm_loader.o Indigo/core/indigo-core/molecule/src/icm_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/icm_saver.o: Indigo/core/indigo-core/molecule/src/icm_saver.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/icm_saver.o Indigo/core/indigo-core/molecule/src/icm_saver.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/inchi_parser.o: Indigo/core/indigo-core/molecule/src/inchi_parser.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/inchi_parser.o Indigo/core/indigo-core/molecule/src/inchi_parser.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/inchi_wrapper.o: Indigo/core/indigo-core/molecule/src/inchi_wrapper.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/inchi_wrapper.o Indigo/core/indigo-core/molecule/src/inchi_wrapper.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/lipinski.o: Indigo/core/indigo-core/molecule/src/lipinski.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/lipinski.o Indigo/core/indigo-core/molecule/src/lipinski.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/max_common_submolecule.o: Indigo/core/indigo-core/molecule/src/max_common_submolecule.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/max_common_submolecule.o Indigo/core/indigo-core/molecule/src/max_common_submolecule.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/metadata_storage.o: Indigo/core/indigo-core/molecule/src/metadata_storage.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/metadata_storage.o Indigo/core/indigo-core/molecule/src/metadata_storage.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule.o: Indigo/core/indigo-core/molecule/src/molecule.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule.o Indigo/core/indigo-core/molecule/src/molecule.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_3d_constraints.o: Indigo/core/indigo-core/molecule/src/molecule_3d_constraints.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_3d_constraints.o Indigo/core/indigo-core/molecule/src/molecule_3d_constraints.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_allene_stereo.o: Indigo/core/indigo-core/molecule/src/molecule_allene_stereo.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_allene_stereo.o Indigo/core/indigo-core/molecule/src/molecule_allene_stereo.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_arom.o: Indigo/core/indigo-core/molecule/src/molecule_arom.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_arom.o Indigo/core/indigo-core/molecule/src/molecule_arom.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_arom_match.o: Indigo/core/indigo-core/molecule/src/molecule_arom_match.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_arom_match.o Indigo/core/indigo-core/molecule/src/molecule_arom_match.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_auto_loader.o: Indigo/core/indigo-core/molecule/src/molecule_auto_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_auto_loader.o Indigo/core/indigo-core/molecule/src/molecule_auto_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_automorphism_search.o: Indigo/core/indigo-core/molecule/src/molecule_automorphism_search.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_automorphism_search.o Indigo/core/indigo-core/molecule/src/molecule_automorphism_search.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_cdxml_loader.o: Indigo/core/indigo-core/molecule/src/molecule_cdxml_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_cdxml_loader.o Indigo/core/indigo-core/molecule/src/molecule_cdxml_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_cdxml_saver.o: Indigo/core/indigo-core/molecule/src/molecule_cdxml_saver.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_cdxml_saver.o Indigo/core/indigo-core/molecule/src/molecule_cdxml_saver.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_chain_fingerprints.o: Indigo/core/indigo-core/molecule/src/molecule_chain_fingerprints.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_chain_fingerprints.o Indigo/core/indigo-core/molecule/src/molecule_chain_fingerprints.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_cip_calculator.o: Indigo/core/indigo-core/molecule/src/molecule_cip_calculator.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_cip_calculator.o Indigo/core/indigo-core/molecule/src/molecule_cip_calculator.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_cis_trans.o: Indigo/core/indigo-core/molecule/src/molecule_cis_trans.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_cis_trans.o Indigo/core/indigo-core/molecule/src/molecule_cis_trans.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_dearom.o: Indigo/core/indigo-core/molecule/src/molecule_dearom.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_dearom.o Indigo/core/indigo-core/molecule/src/molecule_dearom.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_electrons_localizer.o: Indigo/core/indigo-core/molecule/src/molecule_electrons_localizer.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_electrons_localizer.o Indigo/core/indigo-core/molecule/src/molecule_electrons_localizer.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_exact_matcher.o: Indigo/core/indigo-core/molecule/src/molecule_exact_matcher.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_exact_matcher.o Indigo/core/indigo-core/molecule/src/molecule_exact_matcher.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_exact_substructure_matcher.o: Indigo/core/indigo-core/molecule/src/molecule_exact_substructure_matcher.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_exact_substructure_matcher.o Indigo/core/indigo-core/molecule/src/molecule_exact_substructure_matcher.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_fingerprint.o: Indigo/core/indigo-core/molecule/src/molecule_fingerprint.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_fingerprint.o Indigo/core/indigo-core/molecule/src/molecule_fingerprint.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_gross_formula.o: Indigo/core/indigo-core/molecule/src/molecule_gross_formula.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_gross_formula.o Indigo/core/indigo-core/molecule/src/molecule_gross_formula.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_gross_formula_options.o: Indigo/core/indigo-core/molecule/src/molecule_gross_formula_options.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_gross_formula_options.o Indigo/core/indigo-core/molecule/src/molecule_gross_formula_options.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_hash.o: Indigo/core/indigo-core/molecule/src/molecule_hash.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_hash.o Indigo/core/indigo-core/molecule/src/molecule_hash.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_inchi.o: Indigo/core/indigo-core/molecule/src/molecule_inchi.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_inchi.o Indigo/core/indigo-core/molecule/src/molecule_inchi.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_inchi_component.o: Indigo/core/indigo-core/molecule/src/molecule_inchi_component.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_inchi_component.o Indigo/core/indigo-core/molecule/src/molecule_inchi_component.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_inchi_layers.o: Indigo/core/indigo-core/molecule/src/molecule_inchi_layers.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_inchi_layers.o Indigo/core/indigo-core/molecule/src/molecule_inchi_layers.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_inchi_utils.o: Indigo/core/indigo-core/molecule/src/molecule_inchi_utils.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_inchi_utils.o Indigo/core/indigo-core/molecule/src/molecule_inchi_utils.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_ionize.o: Indigo/core/indigo-core/molecule/src/molecule_ionize.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_ionize.o Indigo/core/indigo-core/molecule/src/molecule_ionize.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_json_loader.o: Indigo/core/indigo-core/molecule/src/molecule_json_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_json_loader.o Indigo/core/indigo-core/molecule/src/molecule_json_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_json_saver.o: Indigo/core/indigo-core/molecule/src/molecule_json_saver.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_json_saver.o Indigo/core/indigo-core/molecule/src/molecule_json_saver.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_layered_molecules.o: Indigo/core/indigo-core/molecule/src/molecule_layered_molecules.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_layered_molecules.o Indigo/core/indigo-core/molecule/src/molecule_layered_molecules.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_mass.o: Indigo/core/indigo-core/molecule/src/molecule_mass.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_mass.o Indigo/core/indigo-core/molecule/src/molecule_mass.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_mass_options.o: Indigo/core/indigo-core/molecule/src/molecule_mass_options.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_mass_options.o Indigo/core/indigo-core/molecule/src/molecule_mass_options.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_morgan_fingerprint_builder.o: Indigo/core/indigo-core/molecule/src/molecule_morgan_fingerprint_builder.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_morgan_fingerprint_builder.o Indigo/core/indigo-core/molecule/src/molecule_morgan_fingerprint_builder.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_name_parser.o: Indigo/core/indigo-core/molecule/src/molecule_name_parser.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_name_parser.o Indigo/core/indigo-core/molecule/src/molecule_name_parser.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_neighbourhood_counters.o: Indigo/core/indigo-core/molecule/src/molecule_neighbourhood_counters.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_neighbourhood_counters.o Indigo/core/indigo-core/molecule/src/molecule_neighbourhood_counters.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_pi_systems_matcher.o: Indigo/core/indigo-core/molecule/src/molecule_pi_systems_matcher.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_pi_systems_matcher.o Indigo/core/indigo-core/molecule/src/molecule_pi_systems_matcher.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_rgroups.o: Indigo/core/indigo-core/molecule/src/molecule_rgroups.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_rgroups.o Indigo/core/indigo-core/molecule/src/molecule_rgroups.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_rgroups_composition.o: Indigo/core/indigo-core/molecule/src/molecule_rgroups_composition.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_rgroups_composition.o Indigo/core/indigo-core/molecule/src/molecule_rgroups_composition.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_savers.o: Indigo/core/indigo-core/molecule/src/molecule_savers.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_savers.o Indigo/core/indigo-core/molecule/src/molecule_savers.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_scaffold_detection.o: Indigo/core/indigo-core/molecule/src/molecule_scaffold_detection.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_scaffold_detection.o Indigo/core/indigo-core/molecule/src/molecule_scaffold_detection.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_sgroups.o: Indigo/core/indigo-core/molecule/src/molecule_sgroups.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_sgroups.o Indigo/core/indigo-core/molecule/src/molecule_sgroups.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_standardize.o: Indigo/core/indigo-core/molecule/src/molecule_standardize.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_standardize.o Indigo/core/indigo-core/molecule/src/molecule_standardize.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_standardize_options.o: Indigo/core/indigo-core/molecule/src/molecule_standardize_options.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_standardize_options.o Indigo/core/indigo-core/molecule/src/molecule_standardize_options.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_stereocenter_options.o: Indigo/core/indigo-core/molecule/src/molecule_stereocenter_options.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_stereocenter_options.o Indigo/core/indigo-core/molecule/src/molecule_stereocenter_options.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_stereocenters.o: Indigo/core/indigo-core/molecule/src/molecule_stereocenters.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_stereocenters.o Indigo/core/indigo-core/molecule/src/molecule_stereocenters.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_substructure_matcher.o: Indigo/core/indigo-core/molecule/src/molecule_substructure_matcher.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_substructure_matcher.o Indigo/core/indigo-core/molecule/src/molecule_substructure_matcher.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_chain.o: Indigo/core/indigo-core/molecule/src/molecule_tautomer_chain.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_chain.o Indigo/core/indigo-core/molecule/src/molecule_tautomer_chain.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_enumerator.o: Indigo/core/indigo-core/molecule/src/molecule_tautomer_enumerator.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_enumerator.o Indigo/core/indigo-core/molecule/src/molecule_tautomer_enumerator.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_match.o: Indigo/core/indigo-core/molecule/src/molecule_tautomer_match.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_match.o Indigo/core/indigo-core/molecule/src/molecule_tautomer_match.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_matcher.o: Indigo/core/indigo-core/molecule/src/molecule_tautomer_matcher.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_matcher.o Indigo/core/indigo-core/molecule/src/molecule_tautomer_matcher.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_substructure_matcher.o: Indigo/core/indigo-core/molecule/src/molecule_tautomer_substructure_matcher.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_substructure_matcher.o Indigo/core/indigo-core/molecule/src/molecule_tautomer_substructure_matcher.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_superstructure.o: Indigo/core/indigo-core/molecule/src/molecule_tautomer_superstructure.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_superstructure.o Indigo/core/indigo-core/molecule/src/molecule_tautomer_superstructure.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_utils.o: Indigo/core/indigo-core/molecule/src/molecule_tautomer_utils.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tautomer_utils.o Indigo/core/indigo-core/molecule/src/molecule_tautomer_utils.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tgroups.o: Indigo/core/indigo-core/molecule/src/molecule_tgroups.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molecule_tgroups.o Indigo/core/indigo-core/molecule/src/molecule_tgroups.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molfile_loader.o: Indigo/core/indigo-core/molecule/src/molfile_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molfile_loader.o Indigo/core/indigo-core/molecule/src/molfile_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molfile_saver.o: Indigo/core/indigo-core/molecule/src/molfile_saver.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/molfile_saver.o Indigo/core/indigo-core/molecule/src/molfile_saver.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/monomer_commons.o: Indigo/core/indigo-core/molecule/src/monomer_commons.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/monomer_commons.o Indigo/core/indigo-core/molecule/src/monomer_commons.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/monomers_lib.o: Indigo/core/indigo-core/molecule/src/monomers_lib.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/monomers_lib.o Indigo/core/indigo-core/molecule/src/monomers_lib.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/multiple_cdx_loader.o: Indigo/core/indigo-core/molecule/src/multiple_cdx_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/multiple_cdx_loader.o Indigo/core/indigo-core/molecule/src/multiple_cdx_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/multiple_cml_loader.o: Indigo/core/indigo-core/molecule/src/multiple_cml_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/multiple_cml_loader.o Indigo/core/indigo-core/molecule/src/multiple_cml_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/parse_utils.o: Indigo/core/indigo-core/molecule/src/parse_utils.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/parse_utils.o Indigo/core/indigo-core/molecule/src/parse_utils.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/query_molecule.o: Indigo/core/indigo-core/molecule/src/query_molecule.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/query_molecule.o Indigo/core/indigo-core/molecule/src/query_molecule.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/rdf_loader.o: Indigo/core/indigo-core/molecule/src/rdf_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/rdf_loader.o Indigo/core/indigo-core/molecule/src/rdf_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/sdf_loader.o: Indigo/core/indigo-core/molecule/src/sdf_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/sdf_loader.o Indigo/core/indigo-core/molecule/src/sdf_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/smiles_loader.o: Indigo/core/indigo-core/molecule/src/smiles_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/smiles_loader.o Indigo/core/indigo-core/molecule/src/smiles_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/smiles_saver.o: Indigo/core/indigo-core/molecule/src/smiles_saver.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/smiles_saver.o Indigo/core/indigo-core/molecule/src/smiles_saver.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/structure_checker.o: Indigo/core/indigo-core/molecule/src/structure_checker.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/structure_checker.o Indigo/core/indigo-core/molecule/src/structure_checker.cpp

${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/tpsa.o: Indigo/core/indigo-core/molecule/src/tpsa.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/molecule/src/tpsa.o Indigo/core/indigo-core/molecule/src/tpsa.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/base_reaction.o: Indigo/core/indigo-core/reaction/src/base_reaction.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/base_reaction.o Indigo/core/indigo-core/reaction/src/base_reaction.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/base_reaction_substructure_matcher.o: Indigo/core/indigo-core/reaction/src/base_reaction_substructure_matcher.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/base_reaction_substructure_matcher.o Indigo/core/indigo-core/reaction/src/base_reaction_substructure_matcher.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/canonical_rsmiles_saver.o: Indigo/core/indigo-core/reaction/src/canonical_rsmiles_saver.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/canonical_rsmiles_saver.o Indigo/core/indigo-core/reaction/src/canonical_rsmiles_saver.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/crf_loader.o: Indigo/core/indigo-core/reaction/src/crf_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/crf_loader.o Indigo/core/indigo-core/reaction/src/crf_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/crf_saver.o: Indigo/core/indigo-core/reaction/src/crf_saver.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/crf_saver.o Indigo/core/indigo-core/reaction/src/crf_saver.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/icr_loader.o: Indigo/core/indigo-core/reaction/src/icr_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/icr_loader.o Indigo/core/indigo-core/reaction/src/icr_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/icr_saver.o: Indigo/core/indigo-core/reaction/src/icr_saver.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/icr_saver.o Indigo/core/indigo-core/reaction/src/icr_saver.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/query_reaction.o: Indigo/core/indigo-core/reaction/src/query_reaction.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/query_reaction.o Indigo/core/indigo-core/reaction/src/query_reaction.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction.o: Indigo/core/indigo-core/reaction/src/reaction.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction.o Indigo/core/indigo-core/reaction/src/reaction.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_auto_loader.o: Indigo/core/indigo-core/reaction/src/reaction_auto_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_auto_loader.o Indigo/core/indigo-core/reaction/src/reaction_auto_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_automapper.o: Indigo/core/indigo-core/reaction/src/reaction_automapper.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_automapper.o Indigo/core/indigo-core/reaction/src/reaction_automapper.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_cdxml_loader.o: Indigo/core/indigo-core/reaction/src/reaction_cdxml_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_cdxml_loader.o Indigo/core/indigo-core/reaction/src/reaction_cdxml_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_cdxml_saver.o: Indigo/core/indigo-core/reaction/src/reaction_cdxml_saver.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_cdxml_saver.o Indigo/core/indigo-core/reaction/src/reaction_cdxml_saver.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_cml_loader.o: Indigo/core/indigo-core/reaction/src/reaction_cml_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_cml_loader.o Indigo/core/indigo-core/reaction/src/reaction_cml_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_cml_saver.o: Indigo/core/indigo-core/reaction/src/reaction_cml_saver.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_cml_saver.o Indigo/core/indigo-core/reaction/src/reaction_cml_saver.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_enumerator_state.o: Indigo/core/indigo-core/reaction/src/reaction_enumerator_state.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_enumerator_state.o Indigo/core/indigo-core/reaction/src/reaction_enumerator_state.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_exact_matcher.o: Indigo/core/indigo-core/reaction/src/reaction_exact_matcher.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_exact_matcher.o Indigo/core/indigo-core/reaction/src/reaction_exact_matcher.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_fingerprint.o: Indigo/core/indigo-core/reaction/src/reaction_fingerprint.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_fingerprint.o Indigo/core/indigo-core/reaction/src/reaction_fingerprint.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_gross_formula.o: Indigo/core/indigo-core/reaction/src/reaction_gross_formula.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_gross_formula.o Indigo/core/indigo-core/reaction/src/reaction_gross_formula.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_hash.o: Indigo/core/indigo-core/reaction/src/reaction_hash.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_hash.o Indigo/core/indigo-core/reaction/src/reaction_hash.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_json_loader.o: Indigo/core/indigo-core/reaction/src/reaction_json_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_json_loader.o Indigo/core/indigo-core/reaction/src/reaction_json_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_json_saver.o: Indigo/core/indigo-core/reaction/src/reaction_json_saver.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_json_saver.o Indigo/core/indigo-core/reaction/src/reaction_json_saver.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_neighborhood_counters.o: Indigo/core/indigo-core/reaction/src/reaction_neighborhood_counters.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_neighborhood_counters.o Indigo/core/indigo-core/reaction/src/reaction_neighborhood_counters.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_product_enumerator.o: Indigo/core/indigo-core/reaction/src/reaction_product_enumerator.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_product_enumerator.o Indigo/core/indigo-core/reaction/src/reaction_product_enumerator.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_substructure_matcher.o: Indigo/core/indigo-core/reaction/src/reaction_substructure_matcher.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_substructure_matcher.o Indigo/core/indigo-core/reaction/src/reaction_substructure_matcher.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_transformation.o: Indigo/core/indigo-core/reaction/src/reaction_transformation.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/reaction_transformation.o Indigo/core/indigo-core/reaction/src/reaction_transformation.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/rsmiles_loader.o: Indigo/core/indigo-core/reaction/src/rsmiles_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/rsmiles_loader.o Indigo/core/indigo-core/reaction/src/rsmiles_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/rsmiles_saver.o: Indigo/core/indigo-core/reaction/src/rsmiles_saver.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/rsmiles_saver.o Indigo/core/indigo-core/reaction/src/rsmiles_saver.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/rxnfile_loader.o: Indigo/core/indigo-core/reaction/src/rxnfile_loader.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/rxnfile_loader.o Indigo/core/indigo-core/reaction/src/rxnfile_loader.cpp

${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/rxnfile_saver.o: Indigo/core/indigo-core/reaction/src/rxnfile_saver.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/core/indigo-core/reaction/src/rxnfile_saver.o Indigo/core/indigo-core/reaction/src/rxnfile_saver.cpp

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src/ichilnct.o: Indigo/third_party/inchi/INCHI_API/libinchi/src/ichilnct.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src/ichilnct.o Indigo/third_party/inchi/INCHI_API/libinchi/src/ichilnct.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll.o: Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll.o Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll_a.o: Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll_a.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll_a.o Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll_a.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll_a2.o: Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll_a2.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll_a2.o Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll_a2.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll_b.o: Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll_b.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll_b.o Indigo/third_party/inchi/INCHI_API/libinchi/src/inchi_dll_b.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichi_bns.o: Indigo/third_party/inchi/INCHI_BASE/src/ichi_bns.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichi_bns.o Indigo/third_party/inchi/INCHI_BASE/src/ichi_bns.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichi_io.o: Indigo/third_party/inchi/INCHI_BASE/src/ichi_io.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichi_io.o Indigo/third_party/inchi/INCHI_BASE/src/ichi_io.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichican2.o: Indigo/third_party/inchi/INCHI_BASE/src/ichican2.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichican2.o Indigo/third_party/inchi/INCHI_BASE/src/ichican2.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichicano.o: Indigo/third_party/inchi/INCHI_BASE/src/ichicano.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichicano.o Indigo/third_party/inchi/INCHI_BASE/src/ichicano.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichicans.o: Indigo/third_party/inchi/INCHI_BASE/src/ichicans.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichicans.o Indigo/third_party/inchi/INCHI_BASE/src/ichicans.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichierr.o: Indigo/third_party/inchi/INCHI_BASE/src/ichierr.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichierr.o Indigo/third_party/inchi/INCHI_BASE/src/ichierr.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiisot.o: Indigo/third_party/inchi/INCHI_BASE/src/ichiisot.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiisot.o Indigo/third_party/inchi/INCHI_BASE/src/ichiisot.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichimak2.o: Indigo/third_party/inchi/INCHI_BASE/src/ichimak2.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichimak2.o Indigo/third_party/inchi/INCHI_BASE/src/ichimak2.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichimake.o: Indigo/third_party/inchi/INCHI_BASE/src/ichimake.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichimake.o Indigo/third_party/inchi/INCHI_BASE/src/ichimake.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichimap1.o: Indigo/third_party/inchi/INCHI_BASE/src/ichimap1.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichimap1.o Indigo/third_party/inchi/INCHI_BASE/src/ichimap1.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichimap2.o: Indigo/third_party/inchi/INCHI_BASE/src/ichimap2.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichimap2.o Indigo/third_party/inchi/INCHI_BASE/src/ichimap2.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichimap4.o: Indigo/third_party/inchi/INCHI_BASE/src/ichimap4.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichimap4.o Indigo/third_party/inchi/INCHI_BASE/src/ichimap4.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichinorm.o: Indigo/third_party/inchi/INCHI_BASE/src/ichinorm.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichinorm.o Indigo/third_party/inchi/INCHI_BASE/src/ichinorm.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiparm.o: Indigo/third_party/inchi/INCHI_BASE/src/ichiparm.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiparm.o Indigo/third_party/inchi/INCHI_BASE/src/ichiparm.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiprt1.o: Indigo/third_party/inchi/INCHI_BASE/src/ichiprt1.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiprt1.o Indigo/third_party/inchi/INCHI_BASE/src/ichiprt1.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiprt2.o: Indigo/third_party/inchi/INCHI_BASE/src/ichiprt2.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiprt2.o Indigo/third_party/inchi/INCHI_BASE/src/ichiprt2.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiprt3.o: Indigo/third_party/inchi/INCHI_BASE/src/ichiprt3.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiprt3.o Indigo/third_party/inchi/INCHI_BASE/src/ichiprt3.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiqueu.o: Indigo/third_party/inchi/INCHI_BASE/src/ichiqueu.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiqueu.o Indigo/third_party/inchi/INCHI_BASE/src/ichiqueu.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiread.o: Indigo/third_party/inchi/INCHI_BASE/src/ichiread.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiread.o Indigo/third_party/inchi/INCHI_BASE/src/ichiread.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiring.o: Indigo/third_party/inchi/INCHI_BASE/src/ichiring.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichiring.o Indigo/third_party/inchi/INCHI_BASE/src/ichiring.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr1.o: Indigo/third_party/inchi/INCHI_BASE/src/ichirvr1.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr1.o Indigo/third_party/inchi/INCHI_BASE/src/ichirvr1.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr2.o: Indigo/third_party/inchi/INCHI_BASE/src/ichirvr2.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr2.o Indigo/third_party/inchi/INCHI_BASE/src/ichirvr2.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr3.o: Indigo/third_party/inchi/INCHI_BASE/src/ichirvr3.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr3.o Indigo/third_party/inchi/INCHI_BASE/src/ichirvr3.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr4.o: Indigo/third_party/inchi/INCHI_BASE/src/ichirvr4.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr4.o Indigo/third_party/inchi/INCHI_BASE/src/ichirvr4.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr5.o: Indigo/third_party/inchi/INCHI_BASE/src/ichirvr5.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr5.o Indigo/third_party/inchi/INCHI_BASE/src/ichirvr5.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr6.o: Indigo/third_party/inchi/INCHI_BASE/src/ichirvr6.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr6.o Indigo/third_party/inchi/INCHI_BASE/src/ichirvr6.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr7.o: Indigo/third_party/inchi/INCHI_BASE/src/ichirvr7.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichirvr7.o Indigo/third_party/inchi/INCHI_BASE/src/ichirvr7.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichisort.o: Indigo/third_party/inchi/INCHI_BASE/src/ichisort.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichisort.o Indigo/third_party/inchi/INCHI_BASE/src/ichisort.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichister.o: Indigo/third_party/inchi/INCHI_BASE/src/ichister.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichister.o Indigo/third_party/inchi/INCHI_BASE/src/ichister.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichitaut.o: Indigo/third_party/inchi/INCHI_BASE/src/ichitaut.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ichitaut.o Indigo/third_party/inchi/INCHI_BASE/src/ichitaut.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ikey_base26.o: Indigo/third_party/inchi/INCHI_BASE/src/ikey_base26.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ikey_base26.o Indigo/third_party/inchi/INCHI_BASE/src/ikey_base26.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ikey_dll.o: Indigo/third_party/inchi/INCHI_BASE/src/ikey_dll.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/ikey_dll.o Indigo/third_party/inchi/INCHI_BASE/src/ikey_dll.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/inchi_gui.o: Indigo/third_party/inchi/INCHI_BASE/src/inchi_gui.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/inchi_gui.o Indigo/third_party/inchi/INCHI_BASE/src/inchi_gui.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/mol2atom.o: Indigo/third_party/inchi/INCHI_BASE/src/mol2atom.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/mol2atom.o Indigo/third_party/inchi/INCHI_BASE/src/mol2atom.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt1.o: Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt1.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt1.o Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt1.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt2.o: Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt2.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt2.o Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt2.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt3.o: Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt3.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt3.o Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt3.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt4.o: Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt4.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt4.o Indigo/third_party/inchi/INCHI_BASE/src/mol_fmt4.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/readinch.o: Indigo/third_party/inchi/INCHI_BASE/src/readinch.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/readinch.o Indigo/third_party/inchi/INCHI_BASE/src/readinch.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/runichi.o: Indigo/third_party/inchi/INCHI_BASE/src/runichi.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/runichi.o Indigo/third_party/inchi/INCHI_BASE/src/runichi.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/runichi2.o: Indigo/third_party/inchi/INCHI_BASE/src/runichi2.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/runichi2.o Indigo/third_party/inchi/INCHI_BASE/src/runichi2.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/runichi3.o: Indigo/third_party/inchi/INCHI_BASE/src/runichi3.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/runichi3.o Indigo/third_party/inchi/INCHI_BASE/src/runichi3.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/runichi4.o: Indigo/third_party/inchi/INCHI_BASE/src/runichi4.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/runichi4.o Indigo/third_party/inchi/INCHI_BASE/src/runichi4.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/sha2.o: Indigo/third_party/inchi/INCHI_BASE/src/sha2.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/sha2.o Indigo/third_party/inchi/INCHI_BASE/src/sha2.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/strutil.o: Indigo/third_party/inchi/INCHI_BASE/src/strutil.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/strutil.o Indigo/third_party/inchi/INCHI_BASE/src/strutil.c

${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/util.o: Indigo/third_party/inchi/INCHI_BASE/src/util.c
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src
	${RM} "$@.d"
	$(COMPILE.c) -O3 -DTARGET_API_LIB -IIndigo/core/indigo-core/common/base_c -IIndigo/api/c/indigo -IIndigo/core/indigo-core/common -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/inchi/INCHI_BASE/src/util.o Indigo/third_party/inchi/INCHI_BASE/src/util.c

${OBJECTDIR}/Indigo/third_party/tinyxml2/tinyxml2.o: Indigo/third_party/tinyxml2/tinyxml2.cpp
	${MKDIR} -p ${OBJECTDIR}/Indigo/third_party/tinyxml2
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -DTARGET_API_LIB -IIndigo/third_party/object_threadsafe -IIndigo/third_party/tinyxml2 -IIndigo/third_party/rapidjson -IIndigo/core/indigo-core/common -IIndigo/core/indigo-core/molecule -IIndigo/core/indigo-core/common/base_c -IIndigo/core/indigo-core/common/base_cpp -IIndigo/core/indigo-core -IIndigo/api/c/indigo -IIndigo/third_party/cppcodec/cppcodec -IIndigo/third_party/cppcodec -IIndigo/api/cpp -IIndigo/api/cpp/src -IIndigo/third_party/inchi/INCHI_BASE/src -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Indigo/third_party/tinyxml2/tinyxml2.o Indigo/third_party/tinyxml2/tinyxml2.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
