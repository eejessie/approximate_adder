# This is a c++ code. Using g++ compiler.
CC=g++

# Compiler flags.
CFLAGS= -g

MAINSRC = main.cc 
#OTHSRC1 = proposed.cc enumerate.cc MC.cc previous.cc helper.cc bb_adder.cc #bb_adder_v2.cc #bb_adder_v2_ext.cc
OTHSRC1 = proposed.cc enumerate.cc MC.cc previous.cc helper.cc bb_adder.cc 
OTHSRC2 = compute_ER.cc compute_ED.cc


SRC = $(MAINSRC) $(OTHSRC1) #$(OTHSRC2) 
OBJ = $(SRC:.cc=.o)
TARGET = main

# make all runs.
all: $(TARGET)

%.o: %.cc
	$(CC) $(CFLAGS) -o $@  -c $< 

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJ)  $(LIBS)


# make clean
clean:
	rm -f $(OBJ) $(TARGET)

