#!/bin/bash

#SLOGS="/home/cluster/baandr/infinite_cylinder/slogs"
#SJOBS_FILES="/home/cluster/baandr/infinite_cylinder/sjobs/Hof_V_flow/r_-1/*"
#
#for FILE_PATH in $(ls $SJOBS_FILES | sort -h -t '_' -k13,13n)  # sort by q
#do
#	FILE="$(basename "$FILE_PATH")"  # strip the file from the path
#	sbatch --output=$SLOGS/"$FILE".out --error=$SLOGS/"$FILE".err --job-name="$FILE" --mail-type=ALL \
#	--mail-user=bandrews@physik.uzh.ch --qos=long --time=168:00:00 --mem=185G --ntasks=1 \
#	--cpus-per-task=12 --threads-per-core=1 "$FILE_PATH"
#done

SLOGS="/home/cluster/baandr/infinite_cylinder/slogs"
SJOBS_FILES="/home/cluster/baandr/infinite_cylinder/sjobs/Hof_V_flow/r_1/*"

for FILE_PATH in $(ls $SJOBS_FILES | sort -h -t '_' -k12,12n)  # sort by q
do
	FILE="$(basename "$FILE_PATH")"  # strip the file from the path
	sbatch --output=$SLOGS/"$FILE".out --error=$SLOGS/"$FILE".err --job-name="$FILE" --mail-type=ALL \
	--mail-user=bandrews@physik.uzh.ch --qos=long --time=168:00:00 --mem=64G --ntasks=1 \
	--cpus-per-task=16 --threads-per-core=1 "$FILE_PATH"
done

SLOGS="/home/cluster/baandr/infinite_cylinder/slogs"
SJOBS_FILES="/home/cluster/baandr/infinite_cylinder/sjobs/Hof_V_flow/r_-2/*"

for FILE_PATH in $(ls $SJOBS_FILES | sort -h -t '_' -k12,12n)  # sort by q
do
	FILE="$(basename "$FILE_PATH")"  # strip the file from the path
	sbatch --output=$SLOGS/"$FILE".out --error=$SLOGS/"$FILE".err --job-name="$FILE" --mail-type=ALL \
	--mail-user=bandrews@physik.uzh.ch --qos=long --time=168:00:00 --mem=64G --ntasks=1 \
	--cpus-per-task=16 --threads-per-core=1 "$FILE_PATH"
done

SLOGS="/home/cluster/baandr/infinite_cylinder/slogs"
SJOBS_FILES="/home/cluster/baandr/infinite_cylinder/sjobs/Hof_V_flow/r_2/*"

for FILE_PATH in $(ls $SJOBS_FILES | sort -h -t '_' -k12,12n)  # sort by q
do
	FILE="$(basename "$FILE_PATH")"  # strip the file from the path
	sbatch --output=$SLOGS/"$FILE".out --error=$SLOGS/"$FILE".err --job-name="$FILE" --mail-type=ALL \
	--mail-user=bandrews@physik.uzh.ch --qos=long --time=168:00:00 --mem=64G --ntasks=1 \
	--cpus-per-task=16 --threads-per-core=1 "$FILE_PATH"
done

SLOGS="/home/cluster/baandr/infinite_cylinder/slogs"
SJOBS_FILES="/home/cluster/baandr/infinite_cylinder/sjobs/Hof_V_flow/r_-3/*"

for FILE_PATH in $(ls $SJOBS_FILES | sort -h -t '_' -k12,12n)  # sort by q
do
	FILE="$(basename "$FILE_PATH")"  # strip the file from the path
	sbatch --output=$SLOGS/"$FILE".out --error=$SLOGS/"$FILE".err --job-name="$FILE" --mail-type=ALL \
	--mail-user=bandrews@physik.uzh.ch --qos=long --time=168:00:00 --mem=64G --ntasks=1 \
	--cpus-per-task=16 --threads-per-core=1 "$FILE_PATH"
done

SLOGS="/home/cluster/baandr/infinite_cylinder/slogs"
SJOBS_FILES="/home/cluster/baandr/infinite_cylinder/sjobs/Hof_V_flow/r_3/*"

for FILE_PATH in $(ls $SJOBS_FILES | sort -h -t '_' -k12,12n)  # sort by q
do
	FILE="$(basename "$FILE_PATH")"  # strip the file from the path
	sbatch --output=$SLOGS/"$FILE".out --error=$SLOGS/"$FILE".err --job-name="$FILE" --mail-type=ALL \
	--mail-user=bandrews@physik.uzh.ch --qos=long --time=168:00:00 --mem=64G --ntasks=1 \
	--cpus-per-task=16 --threads-per-core=1 "$FILE_PATH"
done

SLOGS="/home/cluster/baandr/infinite_cylinder/slogs"
SJOBS_FILES="/home/cluster/baandr/infinite_cylinder/sjobs/Hof_V_flow/r_-4/*"

for FILE_PATH in $(ls $SJOBS_FILES | sort -h -t '_' -k12,12n)  # sort by q
do
	FILE="$(basename "$FILE_PATH")"  # strip the file from the path
	sbatch --output=$SLOGS/"$FILE".out --error=$SLOGS/"$FILE".err --job-name="$FILE" --mail-type=ALL \
	--mail-user=bandrews@physik.uzh.ch --qos=long --time=168:00:00 --mem=64G --ntasks=1 \
	--cpus-per-task=16 --threads-per-core=1 "$FILE_PATH"
done
