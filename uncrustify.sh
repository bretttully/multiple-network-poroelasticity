#! /bin/sh

DIRECTORIES="src"
CPP_FILES=`find $DIRECTORIES -name "*.cpp"`
H_FILES=`find $DIRECTORIES -name "*.h"`
files="${CPP_FILES} ${H_FILES}"

mkdir -p out

for item in $files ; do
    dn=$(dirname $item)
    mkdir -p out/$dn
    uncrustify -f $item -c "./uncrustify.cfg" -o indentoutput.tmp
    mv indentoutput.tmp "$item"

done

