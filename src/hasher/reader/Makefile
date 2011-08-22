OBJ_DIR = ../../../obj

SOURCES = SR_InHashTable.c\
		SR_Reference.c \
		SR_HashRegionTable.c

OBJECTS = $(patsubst %, $(OBJ_DIR)/%, $(SOURCES) )


all: $(OBJECTS)

$(OBJECTS): $(SOURCES) $(HEADERS)
	@echo "  * compliing" $(*F).c
	@$(CC) -c $(CFLAGS) $(*F).c -o $@

.PHONY: all
