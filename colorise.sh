#!/usr/bin/bash

while IFS="" read -r folder || [ -n "$folder" ]
do

    touch $folder'/'$folder'.rst'

    while IFS="" read -r line || [ -n "$line" ]
    do
        #echo $line

        if [[ $line = " "* ]]; then
            echo $line
        else
            echo $line >> $folder'/'$folder'.rst'
        fi


    done < $folder'/raw-rst/'$folder'.rst'
done < folderlist.txt