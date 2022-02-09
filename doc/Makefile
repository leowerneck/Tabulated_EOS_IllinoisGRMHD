all: notes.pdf

notes.pdf: notes.tex notes.bib
	pdflatex notes.tex
	biber notes
	pdflatex notes.tex
	pdflatex notes.tex
	rm -rf *.out *.bbl *.bcf *.aux *.xml *.log *.blg *.toc

clean:
	rm -rf *.out *.bbl *.bcf *.aux *.xml *.log *.blg *.toc *.pdf
