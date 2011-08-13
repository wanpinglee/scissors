OBJ_DIR = ../../../obj

CSOURCES = md5.c

COBJECTS = $(CSOURCES:.c=.o)
BUILT_COBJECTS = $(patsubst %, $(OBJ_DIR)/%, $(COBJECTS) )


HEADERS = SR_Types.h \
		SR_Utilities.h
SOURCES = SR_Error.c
OBJECTS = $(SOURCES:.c=.o)
BUILT_OBJECTS = $(patsubst %, $(OBJ_DIR)/%, $(OBJECTS) )




all: $(BUILT_OBJECTS) $(BUILT_COBJECTS)

$(BUILT_OBJECTS): $(SOURCES) $(HEADERS)
	@echo "  * compliing" $(*F).c
	@$(CC) -c $(CFLAGS) $(*F).c -o $@

$(BUILT_COBJECTS): $(CSOURCES)
	@echo "  * compliing" $(*F).c
	@$(CC) -c $(CFLAGS) $(*F).c -o $@

.PHONY: all
