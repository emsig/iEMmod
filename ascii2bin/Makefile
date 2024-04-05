# Makefile

include ../Make_include

ALL: ascii2bin 

SRC	= ascii2bin.c

OBJ	= $(SRC:%.c=%.o)

ascii2bin:	$(OBJ)
	$(CC) $(OPTC) $(CFLAGS) -o ascii2bin $(OBJ) $(LIBS)

install: ascii2bin
	cp ascii2bin $B

clean:
		rm -f core $(OBJ) ascii2bin

realclean: clean
		rm -f $B/ascii2bin




