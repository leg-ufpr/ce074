#-----------------------------------------------------------------------
# List of Rnw files.

RNWFILES=$(ls *.Rnw)
echo $RNWFILES

# Runs rmarkdown::render() in each Rmd file down in the tree.
for RNW in $RNWFILES; do
    FILE="${RNW%.*}"
    echo $FILE
    Rscript -e "require(knitr); knit(\"$RNW\")"
    pdflatex $FILE
    pdflatex $FILE
    pdflatex $FILE
done

#-----------------------------------------------------------------------