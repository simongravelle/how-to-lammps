#!/usr/bin/bash

while IFS="" read -r folder || [ -n "$folder" ]
do
    if [ -f $folder'/'$folder'.rst' ] ; then
        rm $folder'/'$folder'.rst'
    fi
    touch $folder'/'$folder'.rst'
    while IFS="" read -r line || [ -n "$line" ]
    do
        if [[ $line = " "* ]]; then
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
            echo '    '$line >> $folder'/'$folder'.rst'
        else
            echo $line >> $folder'/'$folder'.rst'
        fi
    done < $folder'/raw-rst/'$folder'.rst'
done < folderlist.txt