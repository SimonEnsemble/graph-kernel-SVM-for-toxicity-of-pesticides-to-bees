for f in mol_w_highlights_347_True mol_w_highlights_347_False
do
    rsvg-convert -f pdf -o $f.pdf $f.svg # to convert svg to pdf
    echo "see $f.pdf"
done
