basnam=`basename $1 .tex`
xelatex -shell-escape $1
bibtex $basnam
xelatex -shell-escape $1
xelatex -shell-escape $1
