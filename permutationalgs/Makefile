include ../make.include
default all all-local: countavoiders.o countpatterns.o perm.o permutilities.o
clean:
	rm -rf *.o
countavoiders.o: countavoiders.h ../hashtables/hashdb.h ../hashtables/hashmap.h perm.h
countpatterns.o: permutilities.h ../hashtables/hashdb.h ../hashtables/hashmap.h perm.h
perm.o: perm.h
permutilities.o: ../hashtables/hashdb.h perm.h

# Nothing for check in this directory
check:



