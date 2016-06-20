# hashtables depends on nothing
# datacollection depends on hashtables and permutationalgs
DIRS = hashtables permutationalgs datacollection examples
all: all-local
check clean all-local:
	for x in $(DIRS); do $(MAKE) -C $$x $@; done
.PHONY: check clean all
