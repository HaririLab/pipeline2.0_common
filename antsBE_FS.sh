#!/bin/bash

T1=$1   ##nifti image for T1
prefix=$2 ###Output Directory +prefix....Prefix appended to end of path like #/mnt/BIAC/munin2.dhe.duke.edu/Hariri/DNS.01/Analysis/All_Imaging/DNS0832/highRes_
threads=1
threads=$3
wd=$(echo $pre | rev | cut -d "/" -f2- | rev)
name=$(echo $pre | rev | cut -d "/" -f1 | rev)
cd $wd
freeDir=$wd

export SUBJECTS_DIR=$wd
antsBrainExtraction.sh -d 3 -a $T1 -e /mnt/BIAC/munin2.dhe.duke.edu/Hariri/DNS.01/Analysis/Max/templates/DNS500/DNS500template_MNI.nii.gz -m /mnt/BIAC/munin2.dhe.duke.edu/Hariri/DNS.01/Analysis/Max/templates/DNS500/DNS500template_MNI_BrainCerebellumProbabilityMask.nii.gz -o $prefix

mksubjdirs FreeSurfer
mri_convert ${prefix}BrainExtractionBrain.nii.gz FreeSurfer/mri/001.mgz
recon-all -autorecon1 -noskullstrip -s FreeSurfer -openmp $threads
cp ${freeDir}/mri/T1.mgz ${freeDir}/mri/brainmask.auto.mgz
cp ${freeDir}/mri/brainmask.auto.mgz ${freeDir}/mri/brainmask.mgz
recon-all -autorecon2 -autorecon3 -s FreeSurfer -openmp $threads
