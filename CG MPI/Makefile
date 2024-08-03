CXX := mpicxx -std=c++17 -Wall -Wshadow -Werror -O3 -DNDEBUG  # Add all the flags here


.PHONY: all clean

TARGET:cgsolve

all: $TARGET

$TARGET: $TARGET.cpp Grid.hpp Timer.h
	$(CXX) -o $TARGET $TARGET.cpp -lm

clean:
	rm -f $TARGET
