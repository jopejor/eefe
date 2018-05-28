#!/bin/bash

# Build the files using mcp 


regular="./$1/"
expression="*/*/*/*.csv"
regular_expresion=$regular$expression


echo "$regular_expresion"


mcp "$regular_expresion" '#1_#2_#3_#4.csv'

# generate a file of all of them

filename="$2.csv"

cat *.csv > "$filename"

Rscript "RTS_susbet.r" "$filename"




