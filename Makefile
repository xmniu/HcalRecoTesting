# Root variables
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs) -lMinuit2 
ROOTGLIBS    := $(shell root-config --glibs)
#Note - Minuit and Minuit2 libraries kept at the moment, can remove Minuit after testing

ODIR=objs

# Programs
CXX          = g++
CXXFLAGS     = -Wall
LD	     = g++
LDFLAGS      = -g
SOFLAGS      = -shared

# Assign or Add variables
CXXFLAGS    += $(ROOTCFLAGS) 
LIBS        += $(ROOTLIBS)

CXXRCS       = $(patsubst %.cxx,src/%.cxx,$(CXRCS))

CXXOBJS      = $(patsubst %.cxx,$(ODIR)/%.o,$(CXRCS))

OBJECTS = $(ODIR)/Analysis.o $(ODIR)/readparameters.o $(ODIR)/HybridMinimizer.o $(ODIR)/HcalTimeSlew.o $(ODIR)/HcalPulseShape.o $(ODIR)/PulseShapeFitOOTPileupCorrection.o $(ODIR)/HcalPulseShapes.o $(ODIR)/HLTv2.o $(ODIR)/inverseGaussCDF.o $(ODIR)/PedestalSub.o $(ODIR)/TimeSlewPar.o

all: $(OBJECTS)
	$(CXX) $(LDFLAGS) $(LIBS) $(CXXFLAGS) -o Analysis -g $(OBJECTS) 

#Analyzer/Analysis.h: Analyzer/readparameters/readparameters.h Analyzer/HybridMinimizer.h 

$(ODIR)/Analysis.o: Analyzer/Analysis.cpp Analyzer/Analysis.h Analyzer/readparameters/readparameters.h Analyzer/readparameters/readparameters.cxx Analyzer/HybridMinimizer.h Analyzer/HybridMinimizer.cc Analyzer/HcalPulseShapes.h Analyzer/HcalPulseShapes.cc Analyzer/HcalTimeSlew.h Analyzer/HcalTimeSlew.cc Analyzer/PulseShapeFitOOTPileupCorrection.h Analyzer/PulseShapeFitOOTPileupCorrection.cc Analyzer/isFinite.h Analyzer/HcalPulseShape.h Analyzer/HcalPulseShape.cc Analyzer/PedestalSub.h Analyzer/PedestalSub.cc Analyzer/HLTv2.h Analyzer/HLTv2.cc Analyzer/inverseGaussCDF.hh Analyzer/inverseGaussCDF.cc Analyzer/TimeSlewPar.h Analyzer/TimeSlewPar.cc
	$(CXX) $(CXXFLAGS) -c Analyzer/Analysis.cpp -o $@

$(ODIR)/PulseShapeFitOOTPileupCorrection.o: Analyzer/PulseShapeFitOOTPileupCorrection.cc Analyzer/PulseShapeFitOOTPileupCorrection.h Analyzer/HcalPulseShapes.h Analyzer/HcalPulseShapes.cc Analyzer/HcalPulseShape.h Analyzer/HcalPulseShape.cc
	$(CXX) $(CXXFLAGS) -c Analyzer/PulseShapeFitOOTPileupCorrection.cc -o $@

$(ODIR)/readparameters.o: Analyzer/readparameters/readparameters.h Analyzer/readparameters/readparameters.cxx 
	$(CXX) $(CXXFLAGS) -c Analyzer/readparameters/readparameters.cxx -o $@

$(ODIR)/HybridMinimizer.o: Analyzer/HybridMinimizer.h Analyzer/HybridMinimizer.cc
	$(CXX) $(CXXFLAGS) -c Analyzer/HybridMinimizer.cc -o $@

$(ODIR)/HcalTimeSlew.o: Analyzer/HcalTimeSlew.h Analyzer/HcalTimeSlew.cc
	$(CXX) $(CXXFLAGS) -c Analyzer/HcalTimeSlew.cc -o $@

$(ODIR)/HcalPulseShapes.o: Analyzer/HcalPulseShapes.h Analyzer/HcalPulseShapes.cc Analyzer/HcalPulseShape.h Analyzer/HcalPulseShape.cc
	$(CXX) $(CXXFLAGS) -c Analyzer/HcalPulseShapes.cc -o $@

$(ODIR)/HcalPulseShape.o: Analyzer/HcalPulseShape.h Analyzer/HcalPulseShape.cc
	$(CXX) $(CXXFLAGS) -c Analyzer/HcalPulseShape.cc -o $@

$(ODIR)/HLTv2.o: Analyzer/HLTv2.h Analyzer/HLTv2.cc
	$(CXX) $(CXXFLAGS) -c Analyzer/HLTv2.cc -o $@

$(ODIR)/inverseGaussCDF.o: Analyzer/inverseGaussCDF.hh Analyzer/inverseGaussCDF.cc
	$(CXX) $(CXXFLAGS) -c Analyzer/inverseGaussCDF.cc -o $@

$(ODIR)/PedestalSub.o: Analyzer/PedestalSub.h Analyzer/PedestalSub.cc
	$(CXX) $(CXXFLAGS) -c Analyzer/PedestalSub.cc -o $@

$(ODIR)/TimeSlewPar.o: Analyzer/TimeSlewPar.h Analyzer/TimeSlewPar.cc
	$(CXX) $(CXXFLAGS) -c Analyzer/TimeSlewPar.cc -o $@

clean :
	rm -f $(ODIR)/*.o
	rm -f Analysis 
	rm -f *~;
