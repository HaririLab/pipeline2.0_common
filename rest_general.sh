#!/bin/bash
#
# Script: rest_DNS.sh
# Purpose: Take a minimally preprocessed Resting State Scan and finish preprocessing so that subject is ready for Group Analyses
# Author: Maxwell Elliott

################Steps to include#######################
#1)despike
#2)motion Regress 12 params
#3)censor
#4) bandpass
#5) compcorr
#6) 

###Eventually
#surface
#graph Analysis construction


###############################################################################
#
# Environment set up
#
###############################################################################

# --- BEGIN GLOBAL DIRECTIVE -- 
#$ -o $HOME/$JOB_NAME.$JOB_ID.out
#$ -e $HOME/$JOB_NAME.$JOB_ID.out
#$ -l h_vmem=12G 
# -- END GLOBAL DIRECTIVE -- 

minProcDir=$1
sub=$(echo $minProcDir | rev | cut -d "/" -f2 | rev)
subDir=$(echo $minProcDir | rev | cut -d "/" -f2- | rev)
outDir=${minProcDir}
tmpOutDir=$TMPDIR
tmpDir=${tmpOutDir}/tmp
minProcEpi=${outDir}/epiWarped.nii.gz
templateDir=/mnt/BIAC/munin4.dhe.duke.edu/Hariri/DBIS.01/Analysis/Templates #pipenotes= update/Change away from HardCoding later
templatePre=dunedin115template_MNI_ #pipenotes= update/Change away from HardCoding later
antDir=${subDir}/antCT
antPre="highRes_" #pipenotes= Change away from HardCoding later
FDthresh=.25 #pipenotes= Change away from HardCoding later, also find citations for what you decide likely power 2014, minimun of .5 fd 20DVARS suggested
DVARSthresh=1.55 #pipenotes= Change away from HardCoding later, also find citations for what you decide

echo "----JOB [$JOB_NAME.$JOB_ID] SUBJ $sub START [`date`] on HOST [$HOSTNAME]----"

mkdir -p $tmpDir
##Nest minProc within overarching rest directory
if [[ ! -f ${minProcEpi} ]];then
	echo ""
	echo "!!!!!!!!!!!!!!!!!!!!!!No minimally processed Rest Scan Found!!!!!!!!!!!!!!!!!!!!!!!"
	echo "!!!!!!!!!!!!!!need to run epi_minProc_DBIS.sh first before this script!!!!!!!!!!!!!"
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!EXITING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo ""
	exit
fi

###Extract CompCor Components
voxSize=$(@GetAfniRes ${minProcEpi})
numTR=$(3dinfo -nv ${minProcEpi})
##Check to make sure rest scans are the same size

3dresample -input ${templateDir}/${templatePre}Brain.nii.gz -dxyz $voxSize -prefix ${tmpDir}/refTemplate4epi.nii.gz
antsApplyTransforms -d 3 -t ${antDir}/${antPre}SubjectToTemplate1Warp.nii.gz -t ${antDir}/${antPre}SubjectToTemplate0GenericAffine.mat -o ${tmpDir}/BrainSegmentationPosteriors1Warped2Template.nii.gz -r ${tmpDir}/refTemplate4epi.nii.gz -i ${antDir}/${antPre}BrainSegmentationPosteriors1.nii.gz
antsApplyTransforms -d 3 -t ${antDir}/${antPre}SubjectToTemplate1Warp.nii.gz -t ${antDir}/${antPre}SubjectToTemplate0GenericAffine.mat -o ${tmpDir}/BrainSegmentationPosteriors3Warped2Template.nii.gz -r ${tmpDir}/refTemplate4epi.nii.gz -i ${antDir}/${antPre}BrainSegmentationPosteriors3.nii.gz
3dcalc -a ${tmpDir}/BrainSegmentationPosteriors3Warped2Template.nii.gz -b ${tmpDir}/BrainSegmentationPosteriors1Warped2Template.nii.gz -expr 'step(a-0.95)+step(b-0.95)' -prefix ${tmpDir}/seg.wm.csf.nii.gz
3dmerge -1clust_depth 5 5 -prefix ${tmpDir}/seg.wm.csf.depth.nii.gz ${tmpDir}/seg.wm.csf.nii.gz
3dcalc -a ${tmpDir}/seg.wm.csf.depth.nii.gz -expr 'step(a-1)' -prefix ${tmpDir}/seg.wm.csf.erode.nii.gz ##pipenotes:for DBIS may want to edit this to move further away from WM because of smaller voxels

3dcalc -a ${tmpDir}/seg.wm.csf.erode.nii.gz -b ${outDir}/epiWarped.nii.gz -expr 'a*b' -prefix ${tmpDir}/rest.wm.csf.nii.gz
3dpc -pcsave 5 -prefix ${tmpDir}/pcRest.wm.csf ${tmpDir}/rest.wm.csf.nii.gz
mv ${tmpDir}/pcRest.wm.csf_vec.1D ${tmpOutDir}/
####Setup Censoring
awk -v thresh=$FDthresh '{if($1 > thresh) print NR}' ${outDir}/FD_FSL.1D | awk '{print ($1 - 1) " " $2}' > ${tmpOutDir}/FDcensorTRs.1D #find TRs above threshold and subtract 1 from list to 0 index for afni's liking
awk -v thresh=$DVARSthresh '{if($1 > thresh) print NR}' ${outDir}/DVARS.1D | awk '{print ($1) " " $2}' > ${tmpOutDir}/DVARScensorTRs.1D #find TRs above threshold and Don't subtract 1 from list because DVARS is based on change from first TR and has one less value, value 1 will therefore be for afni 1 index (TR number 2)

cat ${tmpOutDir}/FDcensorTRs.1D ${tmpOutDir}/DVARScensorTRs.1D | sort -g | uniq > ${tmpOutDir}/censorTRs.1D #combine DVARS and FD TRs above threshold 
###cat ${outDir}/pcRest*.wm.csf_vec.1D > ${outDir}/allCompCorr.1D

####Project everything out
####################### replaced allmotion.1D with motion_spm_deg.1D and allmotion_deriv.1D with motion_deriv.1D
clist=$(cat ${tmpOutDir}/censorTRs.1D)
lenC=$(echo $clist | wc -w )
if [[ $lenC == 0 ]];then
	3dTproject -input ${outDir}/epiWarped.nii.gz -mask ${templateDir}/${templatePre}BrainExtractionMask_2mmDil1.nii.gz  -prefix ${tmpOutDir}/epiPrepped_blur6mm.nii.gz -ort ${outDir}/motion.1D -ort ${outDir}/motion_deriv.1D -ort ${tmpOutDir}/pcRest.wm.csf_vec.1D -polort 1 -bandpass 0.008 0.10 -blur 6
##comments: Decided again a more restricted blur in mask with different compartments for cerebellum etc, because that approach seemed to be slighly harming tSNR actually and did not help with peak voxel or extent analyses when applied to Faces contrast. Decided to use a dilated Brain Extraction mask because this at least gets rid of crap that is way outside of brain. This saves space (slightly) and aids with cleaner visualizations. A GM mask can still later be applied for group analyses, this way we at least leave that up to the user.
else
	3dTproject -input ${outDir}/epiWarped.nii.gz -mask ${templateDir}/${templatePre}BrainExtractionMask_2mmDil1.nii.gz -prefix ${tmpOutDir}/epiPrepped_blur6mm.nii.gz -CENSORTR $clist -ort ${outDir}/motion.1D -ort ${outDir}/motion_deriv.1D -ort ${tmpOutDir}/pcRest.wm.csf_vec.1D -polort 1 -bandpass 0.008 0.10 -blur 6
##comments: Decided against a more restricted blur in mask with different compartments for cerebellum etc, because that approach seemed to be slighly harming tSNR actually and did not help with peak voxel or extent analyses when applied to Faces contrast. Decided to use a dilated Brain Extraction mask because this at least gets rid of crap that is way outside of brain. This saves space (slightly) and aids with cleaner visualizations. A GM mask can still later be applied for group analyses, this way we at least leave that up to the user.
fi

rm -r $tmpDir
#extract power 264 time series
mkdir -p ${outDir}/parcellations
/mnt/BIAC/munin4.dhe.duke.edu/Hariri/DNS.01/Analysis/Max/scripts/Pipelines/utils/roi2ts.R -r /mnt/BIAC/munin4.dhe.duke.edu/Hariri/DNS.01/Analysis/Max/templates/Power2011_264/power264_gm10_2mm_myConnectome.nii.gz -i ${tmpOutDir}/epiPrepped_blur6mm.nii.gz > ${outDir}/parcellations/power264.txt
##Now copy all files to the server. 
#Because we run into issues when too many processes are trying to do this in parallel, use a lock dir system to make sure that only a small number of processes are doing this step simultaneously
N_allowed=2
lockDir=/mnt/BIAC/munin4.dhe.duke.edu/Hariri/DBIS.01/Data/ALL_DATA_TO_USE/Imaging/x_x.KEEP.OUT.x_x/locks
if [ ! -e $lockDir ]; then mkdir $lockDir; fi
break=0
while true; do
	for i in `seq 1 $N_allowed`; do
		if mkdir $lockDir/writingBigFiles$i; then
			while true; do
				if mkdir -p $outDir; then
					rsync -v --stats --progress $tmpOutDir/* $outDir/ # check out -W option, --timeout
					echo rsync return code: $?
					rm -r $lockDir/writingBigFiles$i
					break=1
					break
				else
					echo "mkdir $outDir/fslFD25 failed! Sleeping 5"
					sleep 5
				fi
			done
		fi
		if [[ $break -eq 1 ]]; then
			break
		fi
	done
	if [[ $break -eq 1 ]]; then
		break
	else
		sleep 2
	fi
done

# -- BEGIN POST-USER -- 
echo "----JOB [$JOB_NAME.$JOB_ID] STOP [`date`]----" 
mv $HOME/$JOB_NAME.$JOB_ID.out $outDir/fslFD35/$JOB_NAME.$JOB_ID.out	 
# -- END POST-USER -- 


