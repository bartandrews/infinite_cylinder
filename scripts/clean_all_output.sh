#!/bin/bash

cd ../data
find * -type f -size 0 -delete

cd ../logs
find * -type f -size 0 -delete

cd ../pickles
find * -type f -size 0 -delete

