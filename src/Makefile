
EXECUTABLE := cudaCollide

CU_FILES   := VInternal.cu 

CU_DEPS    :=

CC_FILES   := nbody.cpp VCScene.cpp display.cpp
# CC_FILES   := main.cpp VCScene.cpp display.cpp


ARCH=$(shell uname | sed -e 's/-.*//g')

CXX=g++ -m64
CXXFLAGS=-O3 -Wall
LDFLAGS=-L/usr/local/depot/cuda-10.2/lib64/ -lcudart
NVCC=nvcc
NVCCFLAGS=-O3 -m64 --gpu-architecture compute_61 -ccbin /usr/bin/gcc

all: $(EXECUTABLE)

default: $(EXECUTABLE)

OBJS = VInternal.o display.o VCScene.o nbody.o 
# OBJS = VInternal.o display.o VCScene.o main.o

TARGET = libVCollide.a


# $(TARGET) : $(OBJS) 
# 	ar ruv $(TARGET) $(OBJS)
$(EXECUTABLE): $(OBJS)
		$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)
# simple.o : simple.cpp
# 	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)


./%.o: %.cu
		$(NVCC) $< $(NVCCFLAGS) -c -o $@

./%.o: %.cpp
		$(CXX) $< $(CXXFLAGS) -c -o $@

