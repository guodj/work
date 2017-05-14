#!/bin/bash
# Download SSUSI data from HTTPS
for Sate in  "f16 f17 f18" #"f16"
do
    mkdir ./FileName_${Sate}
    PathDir="-P ./FileName_${Sate}"

    # Get files' name & Save it to "./FileName"
    for YY in "2011 2012 2014"
    do
        for DD in {001..366}$
        do
            wget -nc -c "http://ssusi.jhuapl.edu/data_retriver?spc=${Sate}&type=edr-aur&year=${YY}&Doy=${DD}" ${PathDir}
        done
    done

    # Download data files Using './FileName_{Sate}'
    cat ./FileName_${Sate}/* |grep -i "PS.APL" |grep -i "NC" |cut -d '"' -f2 > filenames_${Sate}.log

    cat filenames_${Sate}.log | while read oneline
    do
        echo " "
        echo "Now is downloading:  ${oneline}"
        echo "-------------------------------"
        wget -r -nc -c "http://ssusi.jhuapl.edu/${oneline}"
    done

    echo "Download data is OK ???"

done
