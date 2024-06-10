PCC = g++
PCFLAGS = -std=c++14 -g -fopenmp -march=native -O3

ALL= undirected_unweighted_kBFS undirected_weighted_kBFS undirected_weighted_kBFS_analysis

all: $(ALL)

% : %.cpp
	$(PCC) $(PCFLAGS) -o $@ $< 

.PHONY : clean

clean :
	rm -f *.o $(ALL)

cleansrc :
	rm -f *.o $(ALL)
