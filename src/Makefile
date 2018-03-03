SRC=$(wildcard GenForm-svn/*.cpp)
OBJECTS=$(SRC:.cpp=.o)
DEST=GenForm

.PHONY: all noclean

all: GenForm

GenForm: $(OBJECTS)
	$(CXX) -o $(DEST) $(OBJECTS)

noclean:
	rm -rf $(OBJECTS) $(DEST)
