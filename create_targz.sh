#/bin/sh
filname="ffr-LFDFT-0.2.0.tar.gz"
rm -v ../$filname
tar cvzf ../$filname src/ works/ pseudopotentials/ platform/ tests/
