CXX = g++
CXXFLAGS = -g -Wall -fPIC -I$(MY_CCQENU)
ROOTFLAGS = `root-config --cflags --glibs` -lTreePlayer
CFLAGS = -Wall -fPIC
#OBJS =  CCQENuBinning.o CCQENuCutsNSF.o CCQENuPlotUtils.o CCQENuUtils.o CCQENuEvent.o CCQENuTruth.o GlobalParameters.o NeutronBlob.o NeutronBlobBinning.o NeutronBlobCuts.o ComputeUtils.o
OBJS = CVUniverse.o NukeCC_Cuts.o NukeCCUtilsNSF.o NukeCC_Binning.o
#SRC = $(PWD)/*.cxx
SRC = $(NUKECCSRC_ANA)/src/*.cxx
INCLUDE += -I$(NUKECCSRC_ANA) -I$(PWD)/*.h  -I$(PLOTUTILSROOT)/../include -I$(MINERVAKERNELROOT) -I$(NUKECCSRC_INCLUDE_PATH) -I$(TEST)/ 
LDFLAGS += -shared -Wl,-soname
LDLIBS += -L$(NUKECCSRC_ANA)/src -lNukeCC -lUnfoldUtils

TARGET = libNukeCC.so

$(TARGET): $(OBJS) $(ROOTDICTCXX)
	$(CXX) $(LDFLAGS),$(TARGET) -o $(TARGET) $(OBJS) $(ROOTFLAGS)
	if [ ! -d ${PWD}/lib ]; then mkdir $(PWD)/lib; fi
	ln -sf $(NUKECCSRC_ANA)/src/$(TARGET) ${NUKECCSRC_ANA}/lib/$(TARGET)

$(OBJS): %.o: %.cxx 
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@ 

.PHONY: clean

clean:
	rm -f $(NUKECCSRC_ANA)/lib/$(TARGET)
	rm -f $(TARGET)
	rm -f *.o
