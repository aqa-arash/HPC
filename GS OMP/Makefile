CXX = g++
CXXFLAGS += -std=c++17 -Wall -Winline -Wshadow -fopenmp -O3 -march=native -ffast-math -DNDEBUG

TARGET = gssolve
OBJS = $(TARGET).o

.PHONY: all
all: $(TARGET)

$(TARGET): $(OBJS) Makefile
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

$(TARGET).o: $(TARGET).cpp Makefile
	$(CXX) -c $(CXXFLAGS)  $(TARGET).cpp

.PHONY: clean
clean:
	@$(RM) -rf *.o $(TARGET)
