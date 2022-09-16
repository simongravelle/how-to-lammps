#!/usr/bin/bash

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
            if [[ $line = "    "* ]]; 
            then
                # add to lammps input file
                echo $line >> $inputname
                # write rst file
                echo '    '$line >> $newrst
            elif [[ $line = "# jump" ]]; 
            then
                # add a blank line to lammps input file
                echo "" >> $inputname
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