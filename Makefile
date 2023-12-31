ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)
CXXFLAGS  += $(ROOTCFLAGS)
LIBS       = $(ROOTLIBS) -lASImage
GLIBS      = $(ROOTGLIBS)
GXX	   = g++ -Wall -O3

#simpix_start:  simpix_start.cpp
#	$(GXX) -o simpix_start simpix_start.cpp $(ROOTCFLAGS) $(LIBS) $(ROOTGLIBS)

simpix: simpix.cpp
	$(GXX) -o simpix simpix.cpp $(ROOTCFLAGS) $(LIBS) $(ROOTGLIBS)


clean:
#	rm -f simpix_start out.png
	rm -f simpix out2.png
