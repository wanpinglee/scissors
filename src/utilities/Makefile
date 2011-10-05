# define our build subdirectories
SUB_DIRS = bam \
		hashTable \
		smithwaterman \
		miscellaneous

all:
	@for dir in $(SUB_DIRS); do \
		$(MAKE) --no-print-directory --directory=$$dir; \
	done

.PHONY: all
