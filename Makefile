#include Makefile.def
include Makefile_lc1.def

programObjs = functions.o main.o

main: $(programObjs) 
	$(CXX) $(CXXFLAGS) $(programObjs) $(LIBS) -o main

.cpp.o: $<
	$(CXX) $(CXXFLAGS) -c $<
	
clean:
	rm -f main *.exe *.o *~ *.oa*

main.o: main.cpp main.h functions.h
functions.o: functions.cpp functions.h
