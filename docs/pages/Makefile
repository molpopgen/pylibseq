.SUFFIXES: .rst

%.html: %.ipynb
	PYTHONPATH=../.. jupyter nbconvert --ExecutePreprocessor.timeout=3600 --execute $<

%.md: %.ipynb %.html
	jupyter nbconvert --to markdown $<

%.rst: %.ipynb %.html
	jupyter nbconvert --to rst $<

RST=overview.rst
MD=

all: $(RST) $(MD)

clean: 
	rm -f $(RST)
	rm -f trajectories.nbconvert_files/*

