if [ "$1" = "-h" ]; then
    echo 'Usage:'
    echo '  $1: filename without extension, e.g. `slides`.'
    echo '  $2: {0, 1}, remove LaTeX auxiliary files (.log, .aux, ...).'
else
    FILENAME=$1
    echo "\nConverting Rnw to tex.\n"
    knit $FILENAME.Rnw
    echo "\nConverting tex to PDF.\n"
    pdflatex $FILENAME
    if grep -q '\addbibresource{.*}' "$FILENAME.Rnw"; then
        biber $FILENAME
        pdflatex $FILENAME
    fi
    pdflatex $FILENAME
    pdflatex $FILENAME
fi

if [ "$#" -eq 2 ]; then
    if [ "$2" -eq 1 ]; then
        echo "\nDeleting LaTeX auxiliary files.\n"
        rubber --clean $FILENAME
    fi
fi
