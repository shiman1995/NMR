all: supercell.x
IF=ifort

supercell.x: supercell.f90
	$(IF) -CB $< -o $@

clean:
	rm -f *.x
