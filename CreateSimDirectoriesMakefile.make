#CreateSimDirectoriesMakefile.make
PROG=CreateSimDirectoriesMakefile
SRC=$(PROG).cpp

$(PROG):
	g++ -o $(PROG) $(SRC)

run:
	./$(PROG)