# Root variables
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs) -lMinuit2 
ROOTGLIBS    := $(shell root-config --glibs)
#Note - Minuit and Minuit2 libraries kept at the moment, can remove Minuit after testing

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

CXXOBJS      = $(patsubst %.cxx,%.o,$(CXRCS))

OBJECTS = Analysis.o readparameters.o HybridMinimizer.o HcalTimeSlew.o HcalPulseShape.o PulseShapeFitOOTPileupCorrection.o HcalPulseShapes.o HLTAnalyzer.o

all: $(OBJECTS)
	$(CXX) $(LDFLAGS) $(LIBS) $(CXXFLAGS) -o Analysis -g $(OBJECTS) 

#Analyzer/Analysis.h: Analyzer/readparameters/readparameters.h Analyzer/HybridMinimizer.h 

Analysis.o: Analyzer/Analysis.cpp Analyzer/Analysis.h Analyzer/readparameters/readparameters.h Analyzer/readparameters/readparameters.cxx Analyzer/HybridMinimizer.h Analyzer/HybridMinimizer.cc Analyzer/HcalPulseShapes.h Analyzer/HcalPulseShapes.cc Analyzer/HcalTimeSlew.h Analyzer/HcalTimeSlew.cc Analyzer/PulseShapeFitOOTPileupCorrection.h Analyzer/PulseShapeFitOOTPileupCorrection.cc Analyzer/isFinite.h Analyzer/HcalPulseShape.h Analyzer/HcalPulseShape.cc Analyzer/HLTAnalyzer.h Analyzer/HLTAnalyzer.cc
	$(CXX) $(CXXFLAGS) -c Analyzer/Analysis.cpp
	
PulseShapeFitOOTPileupCorrection.o: Analyzer/PulseShapeFitOOTPileupCorrection.cc Analyzer/PulseShapeFitOOTPileupCorrection.h Analyzer/HcalPulseShapes.h Analyzer/HcalPulseShapes.cc Analyzer/HcalPulseShape.h Analyzer/HcalPulseShape.cc
	$(CXX) $(CXXFLAGS) -c Analyzer/PulseShapeFitOOTPileupCorrection.cc

readparameters.o: Analyzer/readparameters/readparameters.h Analyzer/readparameters/readparameters.cxx 
	$(CXX) $(CXXFLAGS) -c Analyzer/readparameters/readparameters.cxx
	
HybridMinimizer.o: Analyzer/HybridMinimizer.h Analyzer/HybridMinimizer.cc
	$(CXX) $(CXXFLAGS) -c Analyzer/HybridMinimizer.cc
	
HcalTimeSlew.o: Analyzer/HcalTimeSlew.h Analyzer/HcalTimeSlew.cc
	$(CXX) $(CXXFLAGS) -c Analyzer/HcalTimeSlew.cc
	
HcalPulseShapes.o: Analyzer/HcalPulseShapes.h Analyzer/HcalPulseShapes.cc Analyzer/HcalPulseShape.h Analyzer/HcalPulseShape.cc
	$(CXX) $(CXXFLAGS) -c Analyzer/HcalPulseShapes.cc

HcalPulseShape.o: Analyzer/HcalPulseShape.h Analyzer/HcalPulseShape.cc
	$(CXX) $(CXXFLAGS) -c Analyzer/HcalPulseShape.cc
	
HLTAnalyzer.o: Analyzer/HLTAnalyzer.h Analyzer/HLTAnalyzer.cc
	$(CXX) $(CXXFLAGS) -c Analyzer/HLTAnalyzer.cc

clean :
	rm -f *.o
	rm -f Analysis 
	rm -f *~;