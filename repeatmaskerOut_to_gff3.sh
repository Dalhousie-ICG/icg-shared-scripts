#!/bin/bash
#################################
# RM.out to gff3 with Colors!!! #
#################################
#Author: Clément Goubert
#Date created: 08-18-2020
#
# ***Changelog***
#
# Last update by Joran Martijn (clean up code, add coloring for rRNA and tRNA, change coloring for RC/Helitron and simple repeats)
# - v1.2   03-29-2021 #  CURRENT (adds Maverick and Crypton in separate subclass)
# - v1.1   01-18-2021 # (add Penelope and DIRS as separate colors)
# - v1.    08-18-2020
# 
#################################
#
# arguments:
# $1 = RM.out


# usage statement
if [[ $1 == "" ]]; then
    echo "no input provided (RepeatMasker ".out" file)"
    echo "USAGE: ./rm2gff3.sh RM.out"
    exit
fi 


# code body
## create temporary file
## this helps to prevent an unreadable extremely long one-liner
tmp=`mktemp`


# shuffle columns around to fit GFF3 format
## set tab as output field separator
gawk 'BEGIN{OFS="\t"} 
    NR>3 {
        print $5, "RepeatMasker-4.1.2p1", "similarity", $6, $7, $1, $9, ".", "Target="$10" "$11" "$12" "$13" "$14, "Div="$2";Del="$3";Ins="$4";SWscore="$1
    }' $1 > $tmp

# replace C in 'strand' field with -
# replace DNA with TIR
# replace TIR/MAVERICK with MAV/Maverick
# replace TIR/Crypton  with CRY/Crypton
sed -i -e "s/\tC\t/\t-\t/g" -e "s/DNA/TIR/g" -e "s/TIR\/Maverick/MAV\/Maverick/g" -e "s/TIR\/Crypton/CRY\/Crypton/g" $tmp


# omit (left) position
# if feature is in positive strand
## the (left) position field is $13
## else if feature is negative strand
## the (left) position is field $11
gawk -i inplace ' 
    BEGIN{OFS="\t"} {
    if ($7 == "+")
        {print $1,$2,$3,$4,$5,$6,$7,$8,$9" "$10" "$11" "$12";"$NF}
    else 
        {print $1,$2,$3,$4,$5,$6,$7,$8,$9" "$10" "$13" "$12";"$NF}
    }
' $tmp


# add color tag to the 'Attributes' field
# this allows IGV to give different types of repeat elements different colors
gawk -i inplace '
    BEGIN{print "##gff-version 3"} {
    if      ($10 ~ /LINE/)            {print $0";color=#3399ff"}
    else if ($10 ~ /SINE/)            {print $0";color=#800080"}
    else if ($10 ~ /TIR/)             {print $0";color=#ff6666"}
    else if ($10 ~ /LTR/)             {print $0";color=#00cc44"}
    else if ($10 ~ /RC/)              {print $0";color=#0000ff"}
    else if ($10 ~ /Low_complexity/)  {print $0";color=#ffff00"}
    else if ($10 ~ /Satellite/)       {print $0";color=#ff99ff"}
    else if ($10 ~ /Simple_repeat/)   {print $0";color=#8686ac"}
    else if ($10 ~ /Penelope/)        {print $0";color=#b2edba"}
    else if ($10 ~ /DIRS/)            {print $0";color=#fce7bd"}
    else if ($10 ~ /CRY/)             {print $0";color=#8f1800"}
    else if ($10 ~ /MAV/)             {print $0";color=#669999"}
    else if ($10 ~ /rRNA/)            {print $0";color=#ff0000"}
    else if ($10 ~ /tRNA/)            {print $0";color=#ffa500"}
    else                              {print $0";color=#c9c9c9"}
    }
' $tmp


# The repeat type to color mapping:
# /LINE/            {color=#3399ff"} # lightblue
# /SINE/            {color=#800080"} # purple
# /TIR/             {color=#ff6666"} # salmon
# /LTR/             {color=#00cc44"} # green
# /RC/              {color=#0000ff"} # blue
# /Low_complexity/  {color=#ffff00"} # yellow
# /Satellite/       {color=#ff99ff"} # pink
# /Simple_repeat/   {color=#8686ac"} # peri winkle purple
# /Penelope/        {color=#b2edba"} # lightgreen
# /DIRS/            {color=#fce7bd"} # sand
# /CRY/             {color=#8f1800"} # rust
# /MAV/             {color=#669999"} # maroon
# /rRNA/            {color=#ff0000"} # red
# /tRNA/            {color=#ffa500"} # orange
#                   {color=#c9c9c9"} # grey


cat $tmp

