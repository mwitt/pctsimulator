#SimulatePct.make
PROG=SimulatePct
SRC=$(PROG).cpp

$(PROG):
	g++ -o $(PROG) $(SRC)

run:
	./$(PROG)
