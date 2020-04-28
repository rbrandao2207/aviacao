CXX=g++
CXXFLAGS=-std=c++17 -Wall -g3 -I /usr/local/include/boost_1_70_0
LDFLAGS=-lpqxx -lpq -pthread -lboost_serialization

AUX_OBJS=LogitEstimation.o auxFunctions.o

MAIN_OBJS=$(AUX_OBJS) main.o

PROGS=main

all:	$(PROGS)

main:	$(MAIN_OBJS)
						$(CXX) $(CXXFLAGS) -o $@ $(MAIN_OBJS) $(LDFLAGS)
