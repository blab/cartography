TITLE=cartography

rm -f *.log
rm -f *.aux
rm -f *.out
rm -f *.bbl
rm -f *.bcf
rm -f *.blg
rm -f .DS_Store
rm -f ${TITLE}.pdf

pdflatex -draftmode ${TITLE}
bibtex ${TITLE}
pdflatex -draftmode ${TITLE}
pdflatex ${TITLE}

# pdflatex -draftmode ${TITLE}_supplement
# pdflatex -draftmode ${TITLE}_supplement
# pdflatex ${TITLE}_supplement
