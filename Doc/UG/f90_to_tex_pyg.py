from __future__ import print_function

from pygments import highlight
from pygments.lexers import FortranLexer
from pygments.formatters import LatexFormatter

import sys

lines = open(sys.argv[1]).readlines()

PRELINE = '''
\\documentclass[a4paper,10pt,fleqn]{article}
\\usepackage[a4paper]{geometry}
\\geometry{verbose,tmargin=2cm,bmargin=2cm,lmargin=1.5cm,rmargin=1.5cm}

\\setlength{\parskip}{\smallskipamount}
\\setlength{\parindent}{0pt}

\\usepackage{hyperref}
\\usepackage{url}
\\usepackage{fancyvrb}
\\usepackage{xcolor}\n\n

\\definecolor{mygray}{rgb}{0.95,0.95,0.95}
\\usepackage{mdframed}
\\BeforeBeginEnvironment{Verbatim}{\\begin{mdframed}[backgroundcolor=mygray]}
\\AfterEndEnvironment{Verbatim}{\\end{mdframed}}
'''

MIDLINE = '''
\\begin{document}
'''

print(PRELINE)
print( LatexFormatter().get_style_defs() )

print(MIDLINE)

#STRVERBSTART = '''
#\\begin{Verbatim}[commandchars=\\\\\{\},codes={\catcode`\$=3\catcode`\^=7\catcode`\_=8}]
#'''

STRVERBSTART = '''
\\begin{Verbatim}[commandchars=\\\\\{\}]'''

startVerb = False
for l in lines:
    if( not ('!!>' in l) ):
        if( startVerb==False ):
            #print('\\begin{Verbatim}[commandchars=\\\\\{\\},frame=single]')
            print(STRVERBSTART)
            startVerb = True
        str1 = highlight(l, FortranLexer(), LatexFormatter()).split('\n')
        Nlen = len(str1)
        for il in range(1,Nlen-2):
            if str1[il] != '':
                print(str1[il])
    else:
        if( startVerb == True ):
            print('\\end{Verbatim}')
            startVerb = False
        print(l.replace('!!>', '').strip())

if( startVerb == True ):
    print('\\end{Verbatim}')
    startVerb = False


ENDLINE = '''
\\end{document}
'''

print(ENDLINE)
