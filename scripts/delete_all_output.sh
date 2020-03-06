#!/bin/bash

read -r -p "Are you sure that you want to delete all files in data, logs, and pickles? [y/N] " response
if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
then
  cd ../data
  rm -rf '!'(".gitignore")

  cd ../logs
  rm -rf !(".gitignore")

  cd ../pickles
  rm -rf !(".gitignore")
else
    exit
fi

