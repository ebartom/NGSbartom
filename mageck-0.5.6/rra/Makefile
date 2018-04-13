# define the C compiler to use
# CC = gcc
CC = g++

# define any compile-time flags
CFLAGS = -Wall -g -O2

# define any directories containing header files other than /usr/include
#
INCLUDES = -I./include

# define the C source files
APIS = ./src/rngs.cpp ./src/words.cpp ./src/rvgs.cpp ./src/math_api.cpp ./src/fileio.cpp
MAIN1 = ./src/RRA.cpp
# MAIN2 = ./src/CrisprNorm.c

# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
API_OBJS = $(APIS:.cpp=.o)
MAIN1_OBJS = $(MAIN1:.cpp=.o)
MAIN2_OBJS = $(MAIN2:.c=.o)

# define the executable file 
MAIN1_APP = ../bin/RRA
# MAIN2_APP = ../bin/CrisprNorm

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

# all:    $(MAIN1_APP) $(MAIN2_APP)
all:    $(MAIN1_APP) 

$(MAIN1_APP): $(API_OBJS) $(MAIN1_OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN1_APP) $(API_OBJS) $(MAIN1_OBJS) -lm 

# $(MAIN2_APP): $(API_OBJS) $(MAIN2_OBJS)
# 	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN2_APP) $(API_OBJS) $(MAIN2_OBJS) -lm 

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@
.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) $(API_OBJS) $(MAIN1_OBJS) $(MAIN2_OBJS) $(MAIN1_APP) 

depend: $(SRCS)
	makedepend $(INCLUDES) $^
