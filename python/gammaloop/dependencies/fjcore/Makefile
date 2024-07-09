CXX ?= g++

all: libfjcore

libfjcore: fjcore.o wrapper.o
	ar cru libfjcore.a  fjcore.o wrapper.o

clean:
	rm -vf libfjcore.a
 
distclean: clean
	rm -f fjcore.o wrapper.o libfjcore.a

.cc.o:         $<
	$(CXX) $(CXXFLAGS) -O2 -fPIC -c $< -o $@

