# hashtables depends on nothing
# datacollection and examples depend on hashtables and permutationalgs
DIRS = hashtables permutationalgs datacollection examples
all: all-local
check clean all-local:
	for x in $(DIRS); do $(MAKE) -C $$x $@; done
.PHONY: check clean all
