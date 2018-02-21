# shell script for Makefile: run all .R scripts in this directory

for filename in *.R
do
  Rscript $filename
done
