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
CCC=icpx
CXX=icpx
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
	${OBJECTDIR}/main.o


# C Compiler Flags
CFLAGS=-lifcore -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread -lm

# CC Compiler Flags
CCFLAGS=-lifcore -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread -lm
CXXFLAGS=-lifcore -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread -lm

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=../hess-se/dist/Release/GNU-Linux/libhess-se.a ../hess-base/dist/Release/GNU-Linux/libhess-base.a ../indigo/dist/Release/GNU-Linux/libindigo.a ../cmopac/dist/Release/GNU-Linux/libcmopac.a -lpthread -lz

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hess-se-docking

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hess-se-docking: ../hess-se/dist/Release/GNU-Linux/libhess-se.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hess-se-docking: ../hess-base/dist/Release/GNU-Linux/libhess-base.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hess-se-docking: ../indigo/dist/Release/GNU-Linux/libindigo.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hess-se-docking: ../cmopac/dist/Release/GNU-Linux/libcmopac.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hess-se-docking: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	icpx -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hess-se-docking ${OBJECTFILES} ${LDLIBSOPTIONS} -lifcore -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread -lm

${OBJECTDIR}/main.o: main.c
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -O2 -I../hess-base/src/parser -I../hess-base/src -I../hess-base/src/interface -I../hess-base/src/model -I../hess-base/src/bond -I../hess-base/src/empirical-docking -I../indigo/Indigo/api/c/indigo -I../indigo/Indigo/third_party/object_threadsafe -I../indigo/Indigo/third_party/rapidjson -I../indigo/Indigo/third_party/tinyxml2 -I../indigo/Indigo/third_party/inchi/INCHI_BASE/src -I../indigo/Indigo/third_party/cppcodec/cppcodec -I../indigo/Indigo/third_party/cppcodec -I../indigo/Indigo/core/indigo-core -I../indigo/Indigo/core/indigo-core/common/base_c -I../indigo/Indigo/core/indigo-core/common/base_cpp -I../indigo/Indigo/core/indigo-core/molecule -I../indigo/Indigo/core/indigo-core/common -I../indigo/Indigo -I../hess-se/src/interface -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.c

# Subprojects
.build-subprojects:
	cd ../hess-se && ${MAKE}  -f Makefile CONF=Release
	cd ../hess-base && ${MAKE}  -f Makefile CONF=Release
	cd ../indigo && ${MAKE}  -f Makefile CONF=Release
	cd ../cmopac && ${MAKE}  -f Makefile CONF=Release

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:
	cd ../hess-se && ${MAKE}  -f Makefile CONF=Release clean
	cd ../hess-base && ${MAKE}  -f Makefile CONF=Release clean
	cd ../indigo && ${MAKE}  -f Makefile CONF=Release clean
	cd ../cmopac && ${MAKE}  -f Makefile CONF=Release clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
