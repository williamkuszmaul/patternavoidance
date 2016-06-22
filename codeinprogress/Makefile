include ../make.include

# If you invoke make in this directory, make sure the things we depend on are built.
all:
	$(MAKE) -C ../hashtables all-local 
	$(MAKE) all-local
all-local: wordavoidance
wordavoidance: wordavoidance.o ../permutationalgs/perm.o ../permutationalgs/permutilities.o ../hashtables/hashmap.o ../hashtables/hashdb.o
	$(CXX) $(CXXFLAGS) $^ -o $@
clean:
	rm -f wordavoidance *.o
