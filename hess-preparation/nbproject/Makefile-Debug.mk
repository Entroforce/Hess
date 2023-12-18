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
CC=g++
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Debug
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
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=../hess-base/dist/Debug/GNU-Linux/libhess-base.a ../indigo/dist/Debug/GNU-Linux/libindigo.a -lpthread -lz

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hess-preparation

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hess-preparation: ../hess-base/dist/Debug/GNU-Linux/libhess-base.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hess-preparation: ../indigo/dist/Debug/GNU-Linux/libindigo.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hess-preparation: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	g++ -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hess-preparation ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/main.o: main.c
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.c) -g -DTARGET_API_LIB -I../hess-base/src -I../hess-base/src/model -I../hess-base/src/interface -I../Indigo/third_party/rapidjson -I../Indigo/third_party/tinyxml2 -I../Indigo/third_party/object_threadsafe -I../Indigo/third_party/inchi/INCHI_BASE/src -I../Indigo/core/indigo-core/molecule -I../Indigo/core/indigo-core/common/base_cpp -I../Indigo/core/indigo-core/common/base_c -I../Indigo/core/indigo-core/common -I../Indigo/core/indigo-core -I../Indigo/api/c/indigo -I../Indigo -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.c

# Subprojects
.build-subprojects:
	cd ../hess-base && ${MAKE}  -f Makefile CONF=Debug
	cd ../indigo && ${MAKE}  -f Makefile CONF=Debug

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:
	cd ../hess-base && ${MAKE}  -f Makefile CONF=Debug clean
	cd ../indigo && ${MAKE}  -f Makefile CONF=Debug clean

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
