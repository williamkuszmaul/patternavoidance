include ../make.include

# If you invoke make in this directory, make sure the things we depend on are built.
all:
	$(MAKE) -C ../hashtables all-local 
	$(MAKE) -C ../permutationalgs all-local 
	$(MAKE) all-local
all-local: buildpatternsets.o oeislocal.o masstest tabledrawer
masstest: masstest.o oeislocal.o buildpatternsets.o ../permutationalgs/perm.o ../permutationalgs/countavoiders.o ../permutationalgs/permutilities.o ../hashtables/hashmap.o ../hashtables/hashdb.o
	$(CXX) $(CXXFLAGS) $^ -o $@
tabledrawer: tabledrawer.o ../permutationalgs/perm.o ../permutationalgs/countavoiders.o ../permutationalgs/permutilities.o ../hashtables/hashmap.o ../hashtables/hashdb.o
	$(CXX) $(CXXFLAGS) $^ -o $@
clean:
	rm -f example masstest tabledrawer *.o
