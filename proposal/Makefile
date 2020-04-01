.PHONY: all clean

all: howmanybears.pdf

PDFLATEX = pdflatex -interaction=nonstopmode

howmanybears.pdf : references.bib

clean: 
	-rm *.aux *.log *.lof *.lot *.fff *.ttt *.out *.bbl *.blg

%.pdf : %.tex %.bbl
	while ( $(PDFLATEX) $<;  grep -q "Rerun to get" $*.log ) do true ; done

%.aux : %.tex
	-$(PDFLATEX) $<

%.bbl : %.aux
	-bibtex $<

%.svg : %.pdf
	inkscape $< --export-plain-svg=$@

%.png : %.pdf
	convert -density 300 $< -flatten $@

%.pdf : %.ink.svg
	inkscape $< --export-pdf=$@

%.eps : %.pdf
	inkscape --without-gui --export-eps=$@ $<

%.pdf : %.xlsx
	libreoffice --headless --convert-to pdf:calc_pdf_Export $<
	pdfcrop $@ $@
