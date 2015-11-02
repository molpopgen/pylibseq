CXX=c++
CXXFLAGS:=-std=c++11

DIR=check_deps

all: $(DIR)/test.o
	$(CXX) $(CXXFLAGS) -o  $(DIR)/test  $(DIR)/test.o $(LDFLAGS) -lsequence

clean:
	rm -f  $(DIR)/test  $(DIR)/test.o
