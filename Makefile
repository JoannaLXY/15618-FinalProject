APP_NAME=barneshut
OBJS += BHTree.o
OBJS += main.o

CXX = mpic++ -std=c++11
CXXFLAGS = -g -I. -fopenmp -O3 -Wno-unknown-pragmas#-Wall -Wextra

default: $(APP_NAME)

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp
	$(CXX) $< $(CXXFLAGS) -c -o $@

clean:
	/bin/rm -rf *~ *.o $(APP_NAME) *.class
