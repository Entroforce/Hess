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
	${OBJECTDIR}/_ext/58cc6651/elements.o \
	${OBJECTDIR}/src/bond/connectivity.o \
	${OBJECTDIR}/src/bond/functional_groups.o \
	${OBJECTDIR}/src/bond/hybridization.o \
	${OBJECTDIR}/src/bond/ring_operations.o \
	${OBJECTDIR}/src/empirical-docking/fast_gradient_upd.o \
	${OBJECTDIR}/src/empirical-docking/find_ligand_pairs.o \
	${OBJECTDIR}/src/empirical-docking/global_optimizer.o \
	${OBJECTDIR}/src/empirical-docking/ils_primitives.o \
	${OBJECTDIR}/src/empirical-docking/ils_random.o \
	${OBJECTDIR}/src/empirical-docking/matrix4.o \
	${OBJECTDIR}/src/empirical-docking/monte_carlo_entropy.o \
	${OBJECTDIR}/src/empirical-docking/optimizable_molecule.o \
	${OBJECTDIR}/src/empirical-docking/optimization_methods.o \
	${OBJECTDIR}/src/empirical-docking/transform_ph.o \
	${OBJECTDIR}/src/interface/interface.o \
	${OBJECTDIR}/src/model/logger.o \
	${OBJECTDIR}/src/model/molecule.o \
	${OBJECTDIR}/src/model/ph_transformers.o \
	${OBJECTDIR}/src/model/vector3.o \
	${OBJECTDIR}/src/parser/parser.o \
	${OBJECTDIR}/src/scoring_function/features.o \
	${OBJECTDIR}/src/scoring_function/scoring.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libhess-base.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libhess-base.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libhess-base.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libhess-base.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libhess-base.a

${OBJECTDIR}/_ext/58cc6651/elements.o: ../openbabel/elements.cpp
	${MKDIR} -p ${OBJECTDIR}/_ext/58cc6651
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/58cc6651/elements.o ../openbabel/elements.cpp

${OBJECTDIR}/src/bond/connectivity.o: src/bond/connectivity.cpp
	${MKDIR} -p ${OBJECTDIR}/src/bond
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bond/connectivity.o src/bond/connectivity.cpp

${OBJECTDIR}/src/bond/functional_groups.o: src/bond/functional_groups.cpp
	${MKDIR} -p ${OBJECTDIR}/src/bond
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bond/functional_groups.o src/bond/functional_groups.cpp

${OBJECTDIR}/src/bond/hybridization.o: src/bond/hybridization.cpp
	${MKDIR} -p ${OBJECTDIR}/src/bond
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bond/hybridization.o src/bond/hybridization.cpp

${OBJECTDIR}/src/bond/ring_operations.o: src/bond/ring_operations.cpp
	${MKDIR} -p ${OBJECTDIR}/src/bond
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bond/ring_operations.o src/bond/ring_operations.cpp

${OBJECTDIR}/src/empirical-docking/fast_gradient_upd.o: src/empirical-docking/fast_gradient_upd.cpp
	${MKDIR} -p ${OBJECTDIR}/src/empirical-docking
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/empirical-docking/fast_gradient_upd.o src/empirical-docking/fast_gradient_upd.cpp

${OBJECTDIR}/src/empirical-docking/find_ligand_pairs.o: src/empirical-docking/find_ligand_pairs.cpp
	${MKDIR} -p ${OBJECTDIR}/src/empirical-docking
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/empirical-docking/find_ligand_pairs.o src/empirical-docking/find_ligand_pairs.cpp

${OBJECTDIR}/src/empirical-docking/global_optimizer.o: src/empirical-docking/global_optimizer.cpp
	${MKDIR} -p ${OBJECTDIR}/src/empirical-docking
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/empirical-docking/global_optimizer.o src/empirical-docking/global_optimizer.cpp

${OBJECTDIR}/src/empirical-docking/ils_primitives.o: src/empirical-docking/ils_primitives.cpp
	${MKDIR} -p ${OBJECTDIR}/src/empirical-docking
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/empirical-docking/ils_primitives.o src/empirical-docking/ils_primitives.cpp

${OBJECTDIR}/src/empirical-docking/ils_random.o: src/empirical-docking/ils_random.cpp
	${MKDIR} -p ${OBJECTDIR}/src/empirical-docking
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/empirical-docking/ils_random.o src/empirical-docking/ils_random.cpp

${OBJECTDIR}/src/empirical-docking/matrix4.o: src/empirical-docking/matrix4.cpp
	${MKDIR} -p ${OBJECTDIR}/src/empirical-docking
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/empirical-docking/matrix4.o src/empirical-docking/matrix4.cpp

${OBJECTDIR}/src/empirical-docking/monte_carlo_entropy.o: src/empirical-docking/monte_carlo_entropy.cpp
	${MKDIR} -p ${OBJECTDIR}/src/empirical-docking
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/empirical-docking/monte_carlo_entropy.o src/empirical-docking/monte_carlo_entropy.cpp

${OBJECTDIR}/src/empirical-docking/optimizable_molecule.o: src/empirical-docking/optimizable_molecule.cpp
	${MKDIR} -p ${OBJECTDIR}/src/empirical-docking
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/empirical-docking/optimizable_molecule.o src/empirical-docking/optimizable_molecule.cpp

${OBJECTDIR}/src/empirical-docking/optimization_methods.o: src/empirical-docking/optimization_methods.cpp
	${MKDIR} -p ${OBJECTDIR}/src/empirical-docking
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/empirical-docking/optimization_methods.o src/empirical-docking/optimization_methods.cpp

${OBJECTDIR}/src/empirical-docking/transform_ph.o: src/empirical-docking/transform_ph.cpp
	${MKDIR} -p ${OBJECTDIR}/src/empirical-docking
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/empirical-docking/transform_ph.o src/empirical-docking/transform_ph.cpp

${OBJECTDIR}/src/interface/interface.o: src/interface/interface.cpp
	${MKDIR} -p ${OBJECTDIR}/src/interface
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/interface/interface.o src/interface/interface.cpp

${OBJECTDIR}/src/model/logger.o: src/model/logger.cpp
	${MKDIR} -p ${OBJECTDIR}/src/model
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/model/logger.o src/model/logger.cpp

${OBJECTDIR}/src/model/molecule.o: src/model/molecule.cpp
	${MKDIR} -p ${OBJECTDIR}/src/model
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/model/molecule.o src/model/molecule.cpp

${OBJECTDIR}/src/model/ph_transformers.o: src/model/ph_transformers.cpp
	${MKDIR} -p ${OBJECTDIR}/src/model
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/model/ph_transformers.o src/model/ph_transformers.cpp

${OBJECTDIR}/src/model/vector3.o: src/model/vector3.cpp
	${MKDIR} -p ${OBJECTDIR}/src/model
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/model/vector3.o src/model/vector3.cpp

${OBJECTDIR}/src/parser/parser.o: src/parser/parser.cpp
	${MKDIR} -p ${OBJECTDIR}/src/parser
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/parser/parser.o src/parser/parser.cpp

${OBJECTDIR}/src/scoring_function/features.o: src/scoring_function/features.cpp
	${MKDIR} -p ${OBJECTDIR}/src/scoring_function
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/scoring_function/features.o src/scoring_function/features.cpp

${OBJECTDIR}/src/scoring_function/scoring.o: src/scoring_function/scoring.cpp
	${MKDIR} -p ${OBJECTDIR}/src/scoring_function
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -Isrc -Isrc/interface -Isrc/bond -Isrc/model -Isrc/empirical-docking -I../openbabel -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../lbfgspp -I../indigo/Indigo/api/cpp/src -I../indigo/Indigo/api/c/indigo-inchi -I../indigo/Indigo/api/c/indigo-renderer -I../indigo/Indigo/api/c/bingo-nosql -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/scoring_function/scoring.o src/scoring_function/scoring.cpp

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
