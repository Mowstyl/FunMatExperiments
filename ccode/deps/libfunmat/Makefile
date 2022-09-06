TARGET = libfunmat
CC = gcc

IDIR = deps
ODIR = obj
LDIR = headers
SDIR = src
OUTD = bin

LIBS = -lm
CFLAGS = -I$(LDIR) -std=iso9899:1999 -Wall -Wextra -Werror -pedantic-errors -m64

ifeq ($(OS), Windows_NT)
RM = rm -rf
MKDIR_P = mkdir -p
BAR = \\
EXT = .dll
else
RM = rm -rf
MKDIR_P = mkdir
BAR = /
EXT = .so
endif

.PHONY: clean default

default: directories $(TARGET)

directories:
	$(MKDIR_P) $(ODIR)
	$(MKDIR_P) $(OUTD)

$(TARGET): $(ODIR)/funmatrix.o
	$(CC) -shared -o $(OUTD)/$@$(EXT) $(ODIR)/funmatrix.o

$(ODIR)/funmatrix.o: $(LDIR)/precision.h $(LDIR)/funmatrix.h $(SDIR)/funmatrix.c
	$(CC) -c $(CFLAGS) -fpic $(SDIR)/funmatrix.c -o $(ODIR)/funmatrix.o

clean:
	$(RM) $(ODIR)
	$(RM) $(OUTD)

run:
	$(OUTD)/$(TARGET)$(EXT) $(ARGS)
