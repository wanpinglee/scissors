
FLAGS = -Wall -g

LIBS = -lbam -lscissors
LIBS_DIR = -L../outsources/samtools -L./
INCLUDE = -I../


EXECUTABLE = example

all: example.cpp
	@echo "Creating scissors.a"
	@cd ../ && $(MAKE) api --no-print-directory
	@echo "Compiling example.cpp"
	$(CXX) $(FLAGS) -o $(EXECUTABLE) example.cpp $(INCLUDE) $(LIBS_DIR) $(LIBS)

.PHONY: all

clean:
	@cd ../ && $(MAKE) clean --no-print-directory
	@rm example
.PHONY: clean
