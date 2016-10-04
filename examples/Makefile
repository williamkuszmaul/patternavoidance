include ../make.include
all:
	$(MAKE) -C ../permutationalgs
	$(MAKE) -C ../hashtables
	$(MAKE) all-local
all-local: example runtests
example:  example.o  ../permutationalgs/perm.o ../permutationalgs/countavoiders.o ../permutationalgs/countpatterns.o ../permutationalgs/permutilities.o ../hashtables/hashdb.o ../hashtables/hashmap.o
	$(CXX) $(CXXFLAGS) $^ -o $@ 
runtests: runtests.o ../permutationalgs/perm.o ../permutationalgs/countavoiders.o ../permutationalgs/countpatterns.o ../permutationalgs/permutilities.o ../hashtables/hashdb.o ../hashtables/hashmap.o
	$(CXX) $(CXXFLAGS) $^ -o $@ 
example.o:  ../permutationalgs/countavoiders.h ../permutationalgs/countpatterns.h                                    ../permutationalgs/perm.h ../hashtables/hashdb.h
runtests.o: ../permutationalgs/countavoiders.h ../permutationalgs/permutilities.h ../permutationalgs/countpatterns.h ../permutationalgs/perm.h ../hashtables/hashdb.h
clean:
	rm -f example runtests *.o
check: example runtests
	./example
	echo this is broken: ./runtests 

