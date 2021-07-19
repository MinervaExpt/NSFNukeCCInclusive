CXX := g++
ROOTFLAGS := `root-config --cflags --glibs --libs`
CXXFLAGS := -g -Wall -fPIC -I$(PLOTUTILSROOT) 
# Make sure your ROOTFLAGS includes -lTreePlayer, which is where the proxy magic lives
INCLUDE += -I $(CCQENU_ANA)/ -I$(CCQENU_INCLUDE_PATH)
LDLIBS += -L$(PLOTUTILSROOT)/$(CMTCONFIG) -lplotutils -lCintex -L. 


BINARIES = $(shell ls *.cxx | sed 's/\GridCanvas.cxx// ' | sed 's/\GridCanvas_Dict.cxx// ' | sed 's/\.cxx// ')
all : $(BINARIES) libGridCanvas.so


libGridCanvas.so: GridCanvas.o GridCanvas_Dict.o
	$(CXX) $(ROOTFLAGS) -shared -o $@ $^

%.o: %.cxx plot.h
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS)  -c $< -o $@ 

%.o: %.cxx plot_ME.h
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS)  -c $< -o $@ 

GridCanvas_Dict.cxx: GridCanvas.h LinkDef.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

# All of the binaries have the same format, so use a "static pattern
# rule". Each binary "foo" depends on "foo.o" and we build it with the
# recipe given ($@ will be the name of the binary)
$(BINARIES) : %: %.o libGridCanvas.so
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) $(LDLIBS) -lGridCanvas -o $@ $^

clean:
	rm -f $(BINARIES) *.o *.so GridCanvas_Dict.*

.PHONY: clean
