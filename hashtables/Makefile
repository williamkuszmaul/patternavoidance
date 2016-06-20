include ../make.include
default all all-local: hashdb.o hashmap.o hashmaptest
clean:
	rm -f *.o  hashmaptest
hashmaptest: hashmaptest.o hashmap.o
	$(CXX) $(CXXFLAGS) $^ -o $@
check: hashmaptest
	./hashmaptest

hashdb.o: hashdb.h
hashmap.o: hashmap.h
