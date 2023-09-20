#!/bin/sh
"$1" 's/#include \"\.\/.*\//#include \"/' "$2" | "$1"  's/#include \"\.\//#include \"/' > "$3"
