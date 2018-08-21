# Upload.
rsync -avzp \
      ./docs/ \
      --progress \
      --rsh="ssh -p$PATAXOP" "$PATAXO:~/public_html/ensino/CPI/"

# Visit the homepage.
firefox http://leg.ufpr.br/~walmes/ensino/CPI/
