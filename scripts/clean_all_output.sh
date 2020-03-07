#!/bin/bash

for directory in data logs pickles
do
    find ../${directory}/* -type f -size 0 -delete
done

