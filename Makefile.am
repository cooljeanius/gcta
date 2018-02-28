ACLOCAL_AMFLAGS = -I m4 --install

# -----------------------------------------------------------------
#   Makefile for GCTA 
#   
#   Supported platforms
#       Unix / Linux                UNIX
#       Windows                     WIN
#       MacOS                       MAC
# ---------------------------------------------------------------------

# Set this variable to either UNIX, MAC or WIN
SYS = UNIX
OUTPUT = gcta64

# Use sinlge precision to store matrix
SINGLE_PRECISION = 

# Put C++ compiler here; Windows has its own specific version
CXX_UNIX = @CXX@
CXX_WIN = C:\CodeBlocks\MinGW\bin\mingw32-g++.exe
CXX_MAC = @CXX@

# Put the path of library EIGEN here
# Library EIGEN (version 3.0.3) is freely available at http://eigen.tuxfamily.org/
EigenLib = ./eigen

# Any other compiler flags here ( -Wall, -g, etc)
CPPFLAGS = @CPPFLAGS@
CXXFLAGS = @CXXFLAGS@
CXXFLAGS += -w -O3 -I $(EigenLib) -DEIGEN_NO_DEBUG -msse2 -I.
LDFLAGS = @LDFLAGS@

@ifGNUmake@ifdef SINGLE_PRECISION
@ifGNUmake@ CXXFLAGS += -DSINGLE_PRECISION=1
@ifGNUmake@endif

LIB = @LIB@

# Some system specific flags

@ifGNUmake@ifeq ($(SYS),WIN)
@ifGNUmake@ CXXFLAGS += -DWIN -static -I ../../Lib/zlib
@ifGNUmake@ LIB += ../../Lib/zlib/zlib.lib
@ifGNUmake@ CXX = $(CXX_WIN)
@ifGNUmake@endif

@ifGNUmake@ifeq ($(SYS),UNIX)
@ifGNUmake@ CXXFLAGS += -DUNIX -static
@ifGNUmake@ LIB += -lz
@ifGNUmake@ CXX = $(CXX_UNIX)
@ifGNUmake@endif

@ifGNUmake@ifeq ($(SYS),MAC)
@ifGNUmake@ CXXFLAGS += -DUNIX -Dfopen64=fopen @MAINT@ -isysroot /Developer/SDKs/MacOSX10.5.sdk -mmacosx-version-min=10.5
@ifGNUmake@ LIB += -lz
@ifGNUmake@ CXX = $(CXX_MAC)
@ifGNUmake@endif

@ifGNUmake@ifeq ($(SYS),SOLARIS)
@ifGNUmake@ CXX = $(CXX_UNIX)
@ifGNUmake@endif

HDR = @HDR@

HDR += CommFunc.h \
	   cdflib.h \
	   dcdflib.h \
           gcta.h \
	   ipmpar.h \
           StatFunc.h \
           StrFunc.h \
           zfstream.h
SRC = CommFunc.cpp \
           data.cpp \
	   dcdflib.cpp \
           est_hsq.cpp \
           grm.cpp \
           gwas_simu.cpp \
           ld.cpp \
           option.cpp \
           raw_geno.cpp \
           StatFunc.cpp \
           StrFunc.cpp \
           zfstream.cpp

OBJ = $(SRC:.cpp=.o)

all : $(OUTPUT) 

$(OUTPUT) :
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $(OUTPUT) $(OBJ) $(LDFLAGS) $(LIB) 

$(OBJ) : $(HDR)

.cpp.o : 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $*.cpp

SUFFIXES : .cpp .c .o $(SUFFIXES)

$(OUTPUT) : $(OBJ)

FORCE:

clean-local:
	@RM@ -f *.o *~
