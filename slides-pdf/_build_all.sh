#-----------------------------------------------------------------------
# List of Rnw files.

RNWFILES=$(ls *.Rnw)
echo $RNWFILES

rm *.aux *.log *.nav *.out *.snm *.tex *.toc *.bbl *.bcf *.blg *.run.xml *.vrb

# Runs rmarkdown::render() in each Rmd file down in the tree.
for RNW in $RNWFILES; do
    FILE="${RNW%.*}"
    echo $FILE
    Rscript -e "require(knitr); knit(\"$RNW\")"
    pdflatex $FILE
    if grep -q '\addbibresource{.*}' "$RNW"; then
        biber $FILE
        pdflatex $FILE
    fi
    pdflatex $FILE
    pdflatex $FILE
done

rm *.aux *.log *.nav *.out *.snm *.tex *.toc *.bbl *.bcf *.blg *.run.xml *.vrb


#-----------------------------------------------------------------------
