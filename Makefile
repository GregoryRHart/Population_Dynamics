#include Makefile.def
include Makefile_lc1.def

programObjs = functions.o main.o
deviceObjs = d_functions.oa #d_loop.o

main: $(programObjs) $(deviceObjs)

	$(CXX) $(CXXFLAGS) $(programObjs) $(LIBS) $(deviceObjs) -o main

.cpp.o: $<

	$(CXX) $(CXXFLAGS) -c $<
	
d_functions.oa: d_functions.cpp

	$(CXX) -std=c++11 -c d_functions.cpp -o d_functions.oa

d_loop.o: d_loop.cu
	nvcc  $(NVFLAGS) -o d_loop.o -c d_loop.cu -Xptxas -v 

clean:

	rm -f main *.exe *.o *~ *.oa*

main.o: main.cpp main.h functions.h
functions.o: functions.cpp functions.h
d_functions.oa: d_functions.cpp d_functions.h
d_loop.o: d_loop.cu d_functions.h
