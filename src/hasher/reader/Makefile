OBJ_DIR = ../../../obj

SOURCES = SR_InHashTable.c\
		SR_Reference.c \
		SR_HashRegionTable.c

OBJECTS = $(SOURCES:.c=.o)
BUILT_OBJECTS = $(patsubst %, $(OBJ_DIR)/%, $(OBJECTS) )


all: $(BUILT_OBJECTS)

$(BUILT_OBJECTS): $(SOURCES) $(HEADERS)
	@echo "  * compliing" $(*F).c
	@$(CC) -c $(CFLAGS) $(*F).c -o $@

.PHONY: all
