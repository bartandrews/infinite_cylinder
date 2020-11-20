#!/bin/bash

SLOGS="/home/cluster/baandr/infinite_cylinder/slogs"
SJOBS_FILES="/home/cluster/baandr/infinite_cylinder/sjobs/fixed_V/nu_1_3/*"

for FILE_PATH in $(ls $SJOBS_FILES | sort -h -t '_' -k10,10n -k6,6n)  # first Vrange, then chi
do
	FILE="$(basename "$FILE_PATH")"  # strip the file from the path
	sbatch --output=$SLOGS/"$FILE".out --error=$SLOGS/"$FILE".err --job-name="$FILE" --mail-type=ALL \
	--mail-user=bandrews@physik.uzh.ch --qos=long --time=168:00:00 --mem=185G --ntasks=1 \
	--cpus-per-task=12 --threads-per-core=1 "$FILE_PATH"
done

SLOGS="/home/cluster/baandr/infinite_cylinder/slogs"
SJOBS_FILES="/home/cluster/baandr/infinite_cylinder/sjobs/fixed_V/nu_2_5/*"

for FILE_PATH in $(ls $SJOBS_FILES | sort -h -t '_' -k10,10n -k6,6n)  # first Vrange, then chi
do
	FILE="$(basename "$FILE_PATH")"  # strip the file from the path
	sbatch --output=$SLOGS/"$FILE".out --error=$SLOGS/"$FILE".err --job-name="$FILE" --mail-type=ALL \
	--mail-user=bandrews@physik.uzh.ch --qos=long --time=168:00:00 --mem=185G --ntasks=1 \
	--cpus-per-task=12 --threads-per-core=1 "$FILE_PATH"
done

SLOGS="/home/cluster/baandr/infinite_cylinder/slogs"
SJOBS_FILES="/home/cluster/baandr/infinite_cylinder/sjobs/fixed_V/nu_3_7/*"

for FILE_PATH in $(ls $SJOBS_FILES | sort -h -t '_' -k10,10n -k6,6n)  # first Vrange, then chi
do
	FILE="$(basename "$FILE_PATH")"  # strip the file from the path
	sbatch --output=$SLOGS/"$FILE".out --error=$SLOGS/"$FILE".err --job-name="$FILE" --mail-type=ALL \
	--mail-user=bandrews@physik.uzh.ch --qos=long --time=168:00:00 --mem=185G --ntasks=1 \
	--cpus-per-task=12 --threads-per-core=1 "$FILE_PATH"
done