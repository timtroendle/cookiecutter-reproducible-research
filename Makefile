PYTHON = PYTHONPATH=./src python
RESULTS = build/results.pickle
PLOT = build/plot.png

.PHONY : help
help : Makefile
		@sed -n 's/^##//p' $<

build:
	mkdir ./build

## clean : removes all results that have been computed earlier
.PHONY: clean
clean:
	rm -rf ./build/*

## test  : runs the test suite
.PHONY: test
test:
	py.test

$(RESULTS): src/model.py
	$(PYTHON) $< 4 5 $@

$(PLOT): src/vis.py $(RESULTS)
	$(PYTHON) $^ $@

## paper : runs all computational steps and creates the final paper
.PHONY: paper
paper: | build build/paper.pdf

build/paper.pdf: $(PLOT)
build/paper.pdf: report/literature.bib report/main.md report/pandoc-metadata.yml
	cd ./report && \
	pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc \
		main.md pandoc-metadata.yml -t latex -o ../build/paper.pdf
