OBJ_DIR = ../../../obj


HEADERS = SR_Types.h \
		SR_Utilities.h
SOURCES = SR_Error.c
OBJECTS = $(patsubst %.c, $(OBJ_DIR)/%.o, $(SOURCES) )

all: $(OBJECTS) 

$(OBJECTS): $(SOURCES) $(HEADERS)
	@echo "  * compliing" $(*F).c
	@$(CC) -c $(CFLAGS) $(*F).c -o $@


.PHONY: all
