
include make.inc

info:
	@ echo "Author: Hou, Sian"
	@ echo "E-mail:	sianhou1987@outlook.com."
	@ echo "all   : Complier the programs."
	@ echo "clean :	Clean the programs."

all:
	(cd lib && make all)
	(cd src && make all)
	(cd mpisrc && make all)

clean:
	(cd lib && make clean)
	(cd src && make clean)
	(cd mpisrc && make clean)
	(cd bin && rm * -f)
