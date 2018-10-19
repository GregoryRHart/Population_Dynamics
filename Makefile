#include Makefile.def
include Makefile_lc1.def

programObjs = main.o functions.o 

main: $(programObjs)

	$(CXX) $(programObjs) -o main $(LIBS)

.cpp.o: $<

	$(CXX) $(CXXFLAGS) -c $<

clean:

	rm -f main *.exe *.o *~

main.o: main.cpp main.h functions.h
functions.o: functions.cpp functions.h 
