OBJ_DIR = ../../../obj

HEADERS = SR_Types.h \
		SR_Utilities.h

SOURCES = md5.c \
		SR_Error.c \
		SR_Reference.c

OBJECTS = $(SOURCES:.c=.o)

all: $(OBJECTS)

$(OBJECTS): $(SOURCES) $(HEADERS)
	@test -d $(OBJ_DIR) || mkdir $(OBJ_DIR)
	@echo "  * compliing" $(*F).c
	@$(CC) -c -o $(OBJ_DIR)/$@ $(*F).c $(CFLAGS)

.PHONY: all
