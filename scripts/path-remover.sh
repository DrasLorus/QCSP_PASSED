#!/bin/sh
sh -c "$1 's/#include \"\.\/.*\//#include \"/' $2 > $3"
