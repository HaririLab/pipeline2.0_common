#!/bin/bash
#
# Script: anat_DBIS.sh
# Purpose: Pipeline for processing T1 anatomical images datasets that are not in DBIS or DNS. Made originally for processing myConnectome and MSC
# Author: Maxwell Elliott
# Date: 2/24/17
#		10/12/17: ARK added "GLOBAL DIRECTIVE"s to move output to subject antCT dir
#				  ARK added unzip VBM image for use in SPM


###########!!!!!!!!!Pipeline to do!!!!!!!!!!!!!#############
#1)make citations #citations
#2)Follow up on #pipeNotes using ctrl f pipeNotes.... Made these when I knew a trick or something I needed to do later
#3) 3drefit all files in MNI space with -space MNI -view tlrc
#4) maybe add a cut of brain stem, seems to be some issues with this and could add robustness
#5) add optimization to template Brain Mask, 3dzcut then inflate to help with slight cutting of top gyri

###########!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###########################

###############################################################################
#
# Environment set up
#
###############################################################################

# --- BEGIN GLOBAL DIRECTIVE -- 
#$ -o $HOME/$JOB_NAME.$JOB_ID.out
#$ -e $HOME/$JOB_NAME.$JOB_ID.out
#$ -l h_vmem=24G 
# -- END GLOBAL DIRECTIVE -- 

#sub=$1 #$1 or flag -s  #20161103_21449 #pipenotes= Change away from HardCoding later 
TOPDIR=$(findexp DBIS.01)
subDir=$1 #pipenotes= Change away from HardCoding later
sub=$(echo $subDir | rev | cut -d "/" -f1 | rev)
QADir=${subDir}/QA
antDir=${subDir}/antCT
freeDir=$subdir/FreeSurfer/${sub}
tmpDir=${antDir}/tmp
antPre="highRes_" #pipenotes= Change away from HardCoding laterF
templateDir=$TOPDIR/Analysis/Templates #pipenotes= update/Change away from HardCoding later
templatePre=dunedin115template_MNI #pipenotes= update/Change away from HardCoding later
T1=$2
#FLAIR=$3 #don't think I will be doing FS group analyses with this data for now
#T1=$2 #/mnt/BIAC/munin4.dhe.duke.edu/Hariri/DNS.01/Data/Anat/20161103_21449/bia5_21449_006.nii.gz #pipenotes= update/Change away from HardCoding later
threads=$3
if [ ${#threads} -eq 0 ]; then threads=1; fi # antsRegistrationSyN won't work properly if $threads is empty
# baseDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
baseDir=$TOPDIR/Scripts/pipeline2.0_DBIS # using BASH_SOURCE doesn't work for cluster jobs bc they are saved as local copies to nodes
export PATH=$PATH:${baseDir}/scripts/ #add dependent scripts to path #pipenotes= update/Change to DNS scripts
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$threads
export OMP_NUM_THREADS=$threads
export SUBJECTS_DIR=$TOPDIR/Analysis/All_Imaging/FreeSurfer_AllSubs/
export FREESURFER_HOME=${TOPDIR/DBIS/DNS}/Scripts/Tools/FreeSurfer/freesurfer
export ANTSPATH=${TOPDIR/DBIS/DNS}/Scripts/Tools/ants-2.2.0/bin/
export PATH=$PATH:${baseDir}/scripts/:${baseDir}/utils/:$ANTSPATH
echo "----JOB [$JOB_NAME.$JOB_ID] SUBJ $sub START [`date`] on HOST [$HOSTNAME]----" 

##Set up directory
mkdir -p $QADir
cd $subDir
mkdir -p $antDir
mkdir -p $tmpDir



#if [[ ! -f $T1 ]];then
#	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#	echo "!!!!!!!!!!!!!!!!!!!!!NO T1, skipping Anat Processing and Epi processing will also be unavailable!!!!!!!!!!!!!!!"
#	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!EXITING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#	exit
#fi

if [[ ! -f ${antDir}/${antPre}CorticalThicknessNormalizedToTemplate.nii.gz ]];then
	###Rigidly align, to avoid future processing issues
	antsRegistrationSyN.sh -d 3 -t r -f ${templateDir}/${templatePre}.nii.gz -m $T1 -n $threads -o ${antDir}/${antPre}r
	#Make Montage of sub T1 brain extraction to check quality
	echo ""
	echo "#########################################################################################################"
	echo "########################################ANTs Cortical Thickness##########################################"
	echo "#########################################################################################################"
	echo ""
	###Run antCT
	which antsCorticalThickness.sh
	antsCorticalThickness.sh -d 3 -a ${antDir}/${antPre}rWarped.nii.gz -e ${templateDir}/${templatePre}.nii.gz -m ${templateDir}/${templatePre}_BrainCerebellumProbabilityMask.nii.gz -p ${templateDir}/${templatePre}_BrainSegmentationPosteriors%d.nii.gz -t ${templateDir}/${templatePre}_Brain.nii.gz -o ${antDir}/${antPre}
else
	echo ""
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Skipping antCT, Completed Previously!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo ""
fi
##Smooth Cortical Thickness for 2nd level
# ARK updated this from FWHM 8 to 6 on 2/15/18 since that's what had been run in all DNS subjects (DNS used GMSegWarped30 but that shouldnt make much dif)
if [[ ! -f ${antDir}/${antPre}CorticalThicknessNormalizedToTemplate_blur6mm.nii ]];then
	3dBlurInMask -input ${antDir}/${antPre}CorticalThicknessNormalizedToTemplate.nii.gz -mask ${templateDir}/${templatePre}_AvgGMSegWarped25connected.nii.gz -FWHM 6 -prefix ${antDir}/${antPre}CorticalThicknessNormalizedToTemplate_blur6mm.nii.gz
	gunzip ${antDir}/${antPre}CorticalThicknessNormalizedToTemplate_blur6mm.nii.gz # need to unzip so folks can use it in SPM
fi
###Make VBM and smooth
if [[ ! -f ${antDir}/${antPre}JacModVBM_blur8mm.nii.gz ]] && [[ ! -f ${antDir}/${antPre}JacModVBM_blur8mm.nii ]];then
	antsApplyTransforms -d 3 -r ${templateDir}/${templatePre}.nii.gz -i ${antDir}/${antPre}BrainSegmentationPosteriors2.nii.gz -t ${antDir}/${antPre}SubjectToTemplate1Warp.nii.gz -t ${antDir}/${antPre}SubjectToTemplate0GenericAffine.mat -o ${antDir}/${antPre}GMwarped.nii.gz
	antsApplyTransforms -d 3 -r ${templateDir}/${templatePre}.nii.gz -i ${antDir}/${antPre}BrainSegmentationPosteriors4.nii.gz -t ${antDir}/${antPre}SubjectToTemplate1Warp.nii.gz -t ${antDir}/${antPre}SubjectToTemplate0GenericAffine.mat -o ${antDir}/${antPre}SCwarped.nii.gz
	antsApplyTransforms -d 3 -r ${templateDir}/${templatePre}.nii.gz -i ${antDir}/${antPre}BrainSegmentationPosteriors5.nii.gz -t ${antDir}/${antPre}SubjectToTemplate1Warp.nii.gz -t ${antDir}/${antPre}SubjectToTemplate0GenericAffine.mat -o ${antDir}/${antPre}BSwarped.nii.gz
	antsApplyTransforms -d 3 -r ${templateDir}/${templatePre}.nii.gz -i ${antDir}/${antPre}BrainSegmentationPosteriors6.nii.gz -t ${antDir}/${antPre}SubjectToTemplate1Warp.nii.gz -t ${antDir}/${antPre}SubjectToTemplate0GenericAffine.mat -o ${antDir}/${antPre}CBwarped.nii.gz
	3dcalc -a ${antDir}/${antPre}GMwarped.nii.gz -b ${antDir}/${antPre}SCwarped.nii.gz -c ${antDir}/${antPre}CBwarped.nii.gz -d ${antDir}/${antPre}BSwarped.nii.gz -e ${templateDir}/${templatePre}_blurMask25.nii.gz -i ${antDir}/${antPre}SubjectToTemplateLogJacobian.nii.gz -expr '(a*equals(e,1)+b*equals(e,2)+c*equals(e,3)+d*equals(e,4))*i' -prefix ${antDir}/${antPre}JacModVBM.nii.gz
	3dcalc -a ${antDir}/${antPre}GMwarped.nii.gz -b ${antDir}/${antPre}SCwarped.nii.gz -c ${antDir}/${antPre}CBwarped.nii.gz -d ${antDir}/${antPre}BSwarped.nii.gz -e ${templateDir}/${templatePre}_blurMask25.nii.gz -expr '(a*equals(e,1)+b*equals(e,2)+c*equals(e,3)+d*equals(e,4))' -prefix ${antDir}/${antPre}noModVBM.nii.gz
	3dBlurInMask -input ${antDir}/${antPre}JacModVBM.nii.gz -Mmask ${templateDir}/${templatePre}_blurMask25.nii.gz -FWHM 8 -prefix ${antDir}/${antPre}JacModVBM_blur8mm.nii.gz
	3dBlurInMask -input ${antDir}/${antPre}noModVBM.nii.gz -Mmask ${templateDir}/${templatePre}_blurMask25.nii.gz -FWHM 8 -prefix ${antDir}/${antPre}noModVBM_blur8mm.nii.gz
	gunzip ${antDir}/${antPre}JacModVBM_blur8mm.nii.gz # unzip for use in SPM
fi
###Make Brain Extraction QA montages
if [[ ! -f ${QADir}/anat.BrainExtractionCheckAxial.png ]];then
	echo ""
	echo "#########################################################################################################"
	echo "####################################Make QA montages######################################"
	echo "#########################################################################################################"
	echo ""
	##Make Cortical Thickness QA montage
	ConvertScalarImageToRGB 3 ${antDir}/${antPre}CorticalThickness.nii.gz ${tmpDir}/corticalThicknessRBG.nii.gz none red none 0 1 #convert for Ants Montage
	3dcalc -a ${tmpDir}/corticalThicknessRBG.nii.gz -expr 'step(a)' -prefix ${tmpDir}/corticalThicknessRBGstep.nii.gz 
	CreateTiledMosaic -i ${antDir}/${antPre}BrainSegmentation0N4.nii.gz -r ${tmpDir}/corticalThicknessRBG.nii.gz -o ${QADir}/anat.antCTCheck.png -a 0.35 -t -1x-1 -d 2 -p mask -s [5,mask,mask] -x ${tmpDir}/corticalThicknessRBGStep.nii.gz -f 0x1  #Create Montage taking images in axial slices every 5 slices
	ConvertScalarImageToRGB 3 ${antDir}/${antPre}ExtractedBrain0N4.nii.gz ${tmpDir}/highRes_BrainRBG.nii.gz none red none 0 10
	3dcalc -a ${tmpDir}/highRes_BrainRBG.nii.gz -expr 'step(a)' -prefix ${tmpDir}/highRes_BrainRBGstep.nii.gz
	CreateTiledMosaic -i ${antDir}/${antPre}BrainSegmentation0N4.nii.gz -r ${tmpDir}/highRes_BrainRBG.nii.gz -o ${QADir}/anat.BrainExtractionCheckAxial.png -a 0.5 -t -1x-1 -d 2 -p mask -s [5,mask,mask] -x ${tmpDir}/highRes_BrainRBGstep.nii.gz -f 0x1
	CreateTiledMosaic -i ${antDir}/${antPre}BrainSegmentation0N4.nii.gz -r ${tmpDir}/highRes_BrainRBG.nii.gz -o ${QADir}/anat.BrainExtractionCheckSag.png -a 0.5 -t -1x-1 -d 0 -p mask -s [5,mask,mask] -x ${tmpDir}/highRes_BrainRBGstep.nii.gz -f 0x1
fi

### Now run freesurfer
if [[ ! -f ${freeDir}/surf/rh.pial ]];then
	###Prep for Freesurfer with PreSkull Stripped
	#Citation: followed directions from https://surfer.nmr.mgh.harvard.edu/fswiki/UserContributions/FAQ (search skull)
	echo ""
	echo "#########################################################################################################"
	echo "#####################################FreeSurfer Surface Generation#######################################"
	echo "#########################################################################################################"
	echo ""
	###Pipenotes: Currently not doing highRes processing. Can't get it to run without crashing. Also doesn't add that much to our voxels that are already near 1mm^3
	##Set up options file to allow for sub mm voxel high res run of FreeSurfer
	#echo "mris_inflate -n 15" > ${tmpDir}/expert.opts
	#Run
	rm -r ${freeDir}
	cd $TOPDIR/Analysis/All_Imaging/FreeSurfer_AllSubs/
	#mksubjdirs ${sub}
	#cp -R ${FREESURFER_HOME}/subjects/fsaverage ${subDir}/
	echo $freeDir
	#mri_convert ${antDir}/${antPre}ExtractedBrain0N4.nii.gz ${freeDir}/mri/001.mgz
	${FREESURFER_HOME}/bin/recon-all_noLink -all -s $sub -openmp $threads -i ${antDir}/${antPre}rWarped.nii.gz ##Had to edit recon-all to remove soft links in white matter step, links not allowed on BIAC
	#cp ${freeDir}/mri/T1.mgz ${freeDir}/mri/brainmask.auto.mgz
	#cp ${freeDir}/mri/brainmask.auto.mgz ${freeDir}/mri/brainmask.mgz
	#recon-all -autorecon2 -autorecon3 -s $sub -openmp $threads
	recon-all -s $sub -localGI -openmp $threads			
else
	echo ""
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!Skipping FreeSurfer, Completed Previously!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo ""
fi
#Run SUMA
if [[ ! -f ${freeDir}/SUMA/std.60.rh.thickness.niml.dset ]];then
	echo ""
	echo "#########################################################################################################"
	echo "######################################Map Surfaces to SUMA and AFNI######################################"
	echo "#########################################################################################################"
	echo ""
	cd ${freeDir}
	@SUMA_Make_Spec_FS_lgi -NIFTI -ld 60 -sid $sub
	#ConvertDset -o_gii -input ${freeDir}/SUMA/std.60.lh.area.niml.dset -prefix ${freeDir}/SUMA/std.60.lh.area
	#ConvertDset -o_gii -input ${freeDir}/SUMA/std.60.rh.area.niml.dset -prefix ${freeDir}/SUMA/std.60.rh.area
	#ConvertDset -o_gii -input ${freeDir}/SUMA/std.60.lh.thickness.niml.dset -prefix ${freeDir}/SUMA/std.60.lh.thickness
	#ConvertDset -o_gii -input ${freeDir}/SUMA/std.60.rh.thickness.niml.dset -prefix ${freeDir}/SUMA/std.60.rh.thickness
else
	echo ""
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!Skipping SUMA_Make_Spec, Completed Previously!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo ""
fi


### copy files for vis check
cp ${QADir}/anat.BrainExtractionCheckAxial.png $TOPDIR/Graphics/Data_Check/NewPipeline/anat.BrainExtractionCheckAxial/$sub.png
cp ${QADir}/anat.BrainExtractionCheckSag.png $TOPDIR/Graphics/Data_Check/NewPipeline/anat.BrainExtractionCheckSag/$sub.png
cp ${QADir}/anat.antCTCheck.png $TOPDIR/Graphics/Data_Check/NewPipeline/anat.antCTCheck/$sub.png

#cleanup
#mv highRes_* antCT/ #pipeNotes: add more deletion and clean up to minimize space, think about deleting Freesurfer and some of SUMA output
# # # # leave this out for now until completely finished with testing!
rm -r ${antDir}/${antPre}BrainNormalizedToTemplate.nii.gz ${antDir}/${antPre}TemplateToSubject* ${subDir}/dimon.files* ${subDir}/GERT_Reco* 
rm -r ${antDir}/tmp ${freeDir}/SUMA/${sub}_.*spec  ${freeDir}/SUMA/lh.* ${freeDir}/SUMA/rh.*
gzip ${freeDir}/SUMA/*.nii 

 
# -- BEGIN POST-USER -- 
echo "----JOB [$JOB_NAME.$JOB_ID] STOP [`date`]----" 
mv $HOME/$JOB_NAME.$JOB_ID.out $antDir/$JOB_NAME.$JOB_ID.out	 
# -- END POST-USER -- 
