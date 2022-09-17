#!/usr/bin/bash

set -e

while IFS="" read -r folder || [ -n "$folder" ]
do
    newrst=$folder'/'$folder'.rst'
    # errase rst file if exist
    if [ -f $newrst ]
    then
        rm $newrst
    fi
    # create rst file
    touch $newrst
    while IFS="" read -r line || [ -n "$line" ]
    do
        # detect input file starting with '# folder'
        if [[ $line = *"# folder"* ]]
        then
            a=( $line )
            folderinput=${a[2]}
            inputname=$folder$'/files/'$folderinput'/'${a[4]}
            if [ -f $inputname ]
            then
                rm $inputname
            fi
            # write beginning input file
            touch $inputname
            echo '# pure bulk water system' >> $inputname
            echo '# Written by Simon Gravelle' >> $inputname
            echo '# My personal page : https://simongravelle.github.io/' >> $inputname
            echo '# My Github account: https://github.com/simongravelle/' >> $inputname
            echo '# LAMMPS tutorials for beginners: https://lammpstutorials.github.io/' >> $inputname
            echo "" >> $inputname
        else
            if [[ $line = "#space"* ]]; 
            then
                #echo $line
                i=0
                newline=''
                for word in $line
                do
                    #echo $word
                    let i=$i+1
                    if (( $i > 1 )); then
                        newline=$newline' '$word
                    fi
                done
                # add to lammps input file
                echo $newline >> $inputname
                # write rst file
                echo '    '$newline >> $newrst
            elif [[ $line = "# jump" ]]; 
            then
                # add a blank line to lammps input file
                echo "" >> $inputname
            elif [[ $line = *":width:"* ]]; 
            then
                echo '    '$line >> $newrst
            elif [[ $line = *":align:"* ]]; 
            then
                echo '    '$line >> $newrst
            elif [[ $line = *":alt:"* ]]; 
            then
                echo '    '$line >> $newrst
            else
                echo $line >> $newrst
            fi
        fi
    done < $folder'/raw-rst/'$folder'.raw'
done < folderlist.txt

# for coloring some words
#newline=''
#for val in $line
#do   
#    if grep -Fxq $val lammps.constants
#    then 
#        newline=$newline' '$val 
#    else 
#        newline=$newline' '$val 
#    fi
#done
#echo $newline
#echo $newline >> $folder'/'$folder'.rst'