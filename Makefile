CXX := g++
source := HepCalculation.cc

# The `-lm` tag links the C++ math library and `-g` compiles the
# project with debug symbols.
flags := -lm -g

all: HepCalculation

HepCalculation: $(source)
	$(CXX) -o HepCalculation $(source) $(flags)
