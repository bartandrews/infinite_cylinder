#!/bin/bash

read -r -p "Are you sure that you want to delete all files in data, logs, and pickles? [y/N] " response
if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
then
    for directory in data logs pickles
    do
        find ../${directory}/* ! -name '.gitignore' -delete
    done
else
    exit
fi

