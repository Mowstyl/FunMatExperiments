TARGET = libqgates
CC = gcc

IDIR = deps
ODIR = obj
LDIR = headers
SDIR = src
OUTD = bin

LIBS = -lm -lfunmat
CFLAGS = -fopenmp -L$(IDIR) -I$(IDIR) -I$(LDIR) -std=iso9899:1999 -Wall -Wextra -Werror -pedantic-errors -m64

ifeq ($(OS), Windows_NT)
RM = rm -rf
MKDIR_P = mkdir -p
BAR = \\
EXT = .exe
LEXT = .dll
else
RM = rm -rf
MKDIR_P = mkdir -p
BAR = /
EXT = .out
LEXT = .so
endif

.PHONY: clean default directories

default: directories $(TARGET)

directories:
	$(MKDIR_P) $(ODIR)
	$(MKDIR_P) $(OUTD)
	cp $(IDIR)$(BAR)libfunmat$(BAR)bin$(BAR)libfunmat$(LEXT) $(IDIR)$(BAR)libfunmat$(LEXT)
	cp $(IDIR)$(BAR)libfunmat$(BAR)headers$(BAR)funmatrix.h $(IDIR)$(BAR)funmatrix.h

_DEPS = $(wildcard $(LDIR)/*.h)
DEPS = $(patsubst %,%,$(_DEPS))

_LIBH = $(wildcard $(IDIR)/*.h)
LIBH =$(patsubst %,%,$(_LIBH))

_OBJ = $(wildcard $(SDIR)/*.c)
OBJ = $(patsubst $(SDIR)/%.c, $(ODIR)/%.o,$(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS) $(LIBH)
	$(CC) -c $(CFLAGS) -fpic -o $@ $<

$(TARGET): $(ODIR)/qgates.o
	$(CC) -shared -o $(OUTD)/$@$(LEXT) $^ $(CFLAGS) $(LIBS)

clean:
	$(RM) $(ODIR)
	$(RM) $(OUTD)
	$(RM) $(IDIR)$(BAR)libfunmat$(LEXT)
	$(RM) $(IDIR)$(BAR)funmatrix.h

run:
	$(OUTD)/$(TARGET)$(EXT) $(ARGS)
