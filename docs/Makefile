SOURCES = Project.toml make.jl src/index.md src/install.md src/examples/intro.md src/examples/plotting.md src/examples/simulating.md src/examples/modeling.md src/examples/imaging.md 

default: build/index.html

clean:
	rm -f *~ */*~ */*/*~

build/index.html: $(SOURCES)
	julia --color=yes make.jl

.PHOONY: default clean
