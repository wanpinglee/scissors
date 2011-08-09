# define our object and binary directories
OBJ_DIR = ../obj
BIN_DIR = ../bin


# define our source subdirectories
SUB_DIRS = hasher

all:
	@echo "Building......"
	@for dir in $(SUB_DIRS); do \
		$(MAKE) --no-print-directory --directory=$$dir; \
	done

.PHONY: all

clean:
	@echo "Cleaning......"
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*
.PHONY: clean
