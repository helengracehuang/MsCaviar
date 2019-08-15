CC=g++
DIC=$(PWD)
CFLAGS=-c -Wall -g  -I $(DIC) 
LDFLAGS= -I $(DIC)/armadillo/include/ -DARMA_DONT_USE_WRAPPER -llapack -lblas -lgslcblas  -lgsl
SOURCES1=Mcaviar.cpp MPostCal.cpp MUtil.cpp 
EXECUTABLE1=MCAVIAR

all: $(SOURCES1) $(EXECUTABLE1) 
	
$(EXECUTABLE1): $(SOURCES1) 
	$(CC) $(SOURCES1)   $(LDFLAGS) -o $@

clean:
	rm MCAVIAR
