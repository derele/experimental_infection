#!/bin/bash

# #get one header
while read L
do [[ "${L}" == "<Hit>" ]] && break 
    echo $L
done < $1

#loop over each file. Take the fragments  between "<Hit>" and "</Hit>".
for B in "$@" 
do 
    while read L
    do if [[  "${L}" == *"<Hit>"* ]] 
        then echo $L
            while read L
            do echo $L
                if [[  "${L}" == *"</Hit>"* ]] 
                then break
                fi
            done
        fi 
    done < $B 
done 

#get one footer
grep -A 10000 "</Iteration_hits>" $1 
