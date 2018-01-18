include ../../make.inc

SRC = Faddeeva.c

OBJ = $(SRC:.c=.o) 

#
# C source
#
.SUFFIXES: .o .c
.c.o:
	$(CC) -c $(CC_OPTS) $<


# Targets

# we only need the *mod files
lib: $(OBJ)
	cp *.o ../../

# does not delete *.x files
clean:
	rm -rf *.o 




