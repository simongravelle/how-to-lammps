#!/usr/bin/bash

while IFS="" read -r folder || [ -n "$folder" ]
do

    touch $folder'/'$folder'.rst'

    while IFS="" read -r line || [ -n "$line" ]
    do
        #echo $line

        if [[ $line = " "* ]]; then
            newline=''
            for val in $line
            do   
                if grep -Fxq $val lammps.constants
                then 
                    newline=$newline' '$val 
                else 
                    newline=$newline' '$val 
                fi
            done
            echo $newline
        else
            echo $line >> $folder'/'$folder'.rst'
        fi


    done < $folder'/raw-rst/'$folder'.rst'
done < folderlist.txt