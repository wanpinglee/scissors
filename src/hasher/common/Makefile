OBJ_DIR = ../../../obj

HEADERS = SR_Types.h \
		SR_Utilities.h

SOURCES = md5.c \
		SR_Error.c \
		SR_Reference.c

OBJECTS = $(SOURCES:.c=.o)
BUILT_OBJECTS = $(patsubst %, $(OBJ_DIR)/%, $(OBJECTS) )


all: $(BUILT_OBJECTS)

$(BUILT_OBJECTS): $(SOURCES) $(HEADERS)
	@echo "  * compliing" $(*F).c
	@$(CC) -c $(CFLAGS) $(*F).c -o $@

.PHONY: all
