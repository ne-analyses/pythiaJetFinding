os = $(shell uname -s)

#INCFLAGS      = -I$(ROOTSYS)/include -I$(FASTJETDIR)/include
INCFLAGS      = -I$(ROOTSYS)/include -I$(FASTJETDIR)/include -I/opt/local/include -I$(PYTHIA8DIR)/include

ifeq ($(os),Linux)
CXXFLAGS      = -std=c++11
else
CXXFLAGS      = -O -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init -std=c++11
## for debugging:
# CXXFLAGS      = -g -O0 -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
endif

ifeq ($(os),Linux)
LDFLAGS       = -g
LDFLAGSS      = -g --shared 
else
LDFLAGS       = -O -Xlinker -bind_at_load -flat_namespace
LDFLAGSS      = -flat_namespace -undefined suppress
LDFLAGSSS     = -bundle
endif

ifeq ($(os),Linux)
CXX          = g++
else
CXX          = clang
endif


ROOTLIBS      = $(shell root-config --libs)

LIBPATH       = $(ROOTLIBS) -L$(FASTJETDIR)/lib -L$(PYTHIA8DIR)/lib
LIBS          = -lfastjet -lfastjettools -lpythia8


# for cleanup
SDIR          = src
ODIR          = src/obj
BDIR          = bin


###############################################################################
################### Remake when these headers are touched #####################
###############################################################################
jetfinder : g++ -std=c++11 src/jetFindAnalysis.cxx /wsu/home/dx/dx54/dx5412/software/pythia8219/lib/libpythia8.a -o src/obj/jetFindAnalysis.o -I/wsu/home/dx/dx54/dx5412/software/pythia8219/include -I/nfs/rhi/STAR/software/root5/include -I/opt/local/include -I/nfs/rhi/STAR/software/fastjet3/include -O2 -ansi  -fPIC -Wl,-rpath,../lib -ldl

###############################################################################
# standard rules
$(ODIR)/%.o : $(SDIR)/%.cxx $(INCS)
	@echo 
	@echo COMPILING
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

$(BDIR)/%  : $(ODIR)/%.o 
	@echo 
	@echo LINKING
	$(CXX) $(LDFLAGS) $(LIBPATH) $(LIBS) $^ -o $@

###############################################################################

###############################################################################
############################# Main Targets ####################################
###############################################################################
all : $(BDIR)/jetFindAnalysis

#$(ODIR)/qa_v1.o 		: $(SDIR)/qa_v1.cxx
$(ODIR)/jetFindAnalysis.o      : $(SDIR)/jetFindAnalysis.cxx

#data analysis
#$(BDIR)/qa_v1		: $(ODIR)/qa_v1.o
$(BDIR)/jetFindAnalysis     : $(ODIR)/jetFindAnalysis.o

###############################################################################
##################################### MISC ####################################
###############################################################################

clean :
	@echo 
	@echo CLEANING
	rm -vf $(ODIR)/*.o
	rm -vf $(BDIR)/*
	rm -vf lib/*


