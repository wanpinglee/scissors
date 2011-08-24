OBJ_DIR = ../../../obj

SOURCES = SR_InHashTable.c\
		SR_Reference.c \
		SR_HashRegionTable.c

OBJECTS = $(patsubst %.c, $(OBJ_DIR)/%.o, $(SOURCES) )
INCLUDE = -I../..

all: $(OBJECTS)

$(OBJECTS): $(SOURCES) $(HEADERS)
	@echo "  * compliing" $(*F).c
	@$(CC) -c $(CFLAGS) $(*F).c $(INCLUDE) -o $@

.PHONY: all
