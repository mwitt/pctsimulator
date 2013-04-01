#MainMakefile.make

RunSimulator:
	make -f CreateSimDirectoriesMakefile.make
	make -f CreateSimDirectoriesMakefile.make run 
	make -f CreateSimDirectories.make
	make -f SimulatePct.make
	make -f SimulatePct.make run

clean:
	rm -f CreateSimDirectories.make
	rm -f CreateSimDirectoriesMakefile
	rm -f SimulatePct
