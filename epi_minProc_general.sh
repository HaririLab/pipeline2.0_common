#!/bin/bash
#
# Script: epi_minProc_DBIS.sh
# Purpose: Pipeline for minProc of all EPI data for DBIS
# Author: Maxwell Elliott
# Date: 06/26/17


###########!!!!!!!!!Pipeline to do!!!!!!!!!!!!!#############
#1)make citations #citations
#2)Follow up on #pipeNotes using ctrl f pipeNotes.... Made these when I knew a trick or something I needed to do later


###########!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###########################

###############################################################################
#
# Environment set up
#
###############################################################################

####Changes to be made from single band pipeline
#1) Slice Time Correction: Need a more complicated system to account for slice acquisition. Afni has tool to make volume with timing information
#2) Re-test alignment with EPI, may be alternate optimized parameters
#3) Bo unwarping, consider using epi_reg approach https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;6b6c4ba1.1404

# --- BEGIN GLOBAL DIRECTIVE -- 
#$ -o $HOME/$JOB_NAME.$JOB_ID.out
#$ -e $HOME/$JOB_NAME.$JOB_ID.out
#$ -l h_vmem=16G 
# -- END GLOBAL DIRECTIVE -- 

subDir=$1 #pipenotes= Change away from HardCoding later
outDirName=$2 #just a name not a path, this should correspond to the raw data file structure, meant to be flexible so it can be coherent within a dataset
epi=$3 #full path to epi that you want processed
sub=$(echo $subDir | rev | cut -d "/" -f1 | rev)
scriptDir=/mnt/BIAC/munin4.dhe.duke.edu/Hariri/DBIS.01/Scripts/pipeline2.0_DBIS # using BASH_SOURCE doesn't work for cluster jobs bc they are saved as local copies to nodes
QADir=${subDir}/QA
outDir=${subDir}/$outDirName
tmpOutDir=$TMPDIR
tmpDir=${tmpOutDir}/tmp
antDir=${subDir}/antCT
freeDir=${subDir}/FreeSurfer
antPre="highRes_"
templateDir=/mnt/BIAC/munin4.dhe.duke.edu/Hariri/DBIS.01/Analysis/Templates #pipenotes= update/Change away from HardCoding later
templatePre=dunedin115template_MNI_ #pipenotes= update/Change away from HardCoding later
#T1=$2 #/mnt/BIAC/munin4.dhe.duke.edu/Hariri/DNS.01/Data/Anat/20161103_21449/bia5_21449_006.nii.gz #pipenotes= update/Change away from HardCoding later
threads=$4
if [ ${#threads} -eq 0 ]; then threads=1; fi 
fieldMapDir=/mnt/BIAC/munin4.dhe.duke.edu/Hariri/DBIS.01/Data/OTAGO/${sub}/DMHDS/MR_gre_field_mapping_2mm/
export PATH=$PATH:$scriptDir/scripts/ #add dependent scripts to path #pipenotes= update/Change to DNS scripts
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$threads
export OMP_NUM_THREADS=$threads

echo "----JOB [$JOB_NAME.$JOB_ID] SUBJ $sub START [`date`] on HOST [$HOSTNAME]----"

# # ##Make sure anatomical has been run
# # if [ ! -e ${antDir}/${antPre}BrainSegmentation.nii.gz ]; then
	# # echo "#########################################################################################################"
	# # echo "###################### Processed anatomical not found, running anat_DBIS.sh #############################"
	# # echo "#########################################################################################################"
	# # sh /mnt/BIAC/munin4.dhe.duke.edu/Hariri/DBIS.01/Scripts/NewPipeline/anat_DBIS.sh $sub 1
# # fi

##Grab Epi and set up directories
##Set up directory
mkdir -p $QADir
mkdir -p $outDir
mkdir $tmpDir
cd $tmpOutDir

if [[ -f ${outDir}/epiWarped.nii.gz ]];then
	echo ""
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!epi_minProc_general has been Completed Previously!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo ""
	exit
fi

echo ""
echo "#########################################################################################################"
echo "#########################################Prep epi for warping############################################"
echo "#########################################################################################################"
echo ""
rsync -v --stats --progress $epi ${tmpDir}/epi.nii.gz
###Resample structural to voxel dimensions of epi for grid when applying warps
3dDespike -prefix ${tmpDir}/epi_d.nii.gz ${tmpDir}/epi.nii.gz #citation:Jo et al., 2014 #pipeNotes: do we want to do this on tasks?? #citation: Kalcher et al., 2013. Example of despiking in task analysis...at least some people do this... #citation: also https://afni.nimh.nih.gov/afni/community/board/read.php?1,141185,143682#msg-143682, small comment "helpful, more important in rest than task" 
#3dTshift -tzero 0 -prefix ${tmpDir}/epi_dt.nii.gz ${tmpDir}/epi_d.nii.gz #perform t-shifting and shift times to begining of TR as suggested by afni "As the comment states, the goal is to resample each voxel time series so that it is as if each volume was acquired at the beginning of each TR.  That way the slice timing will accurately match the stimulus timing." #citation: https://afni.nimh.nih.gov/pub/dist/edu/data/CD.expanded/AFNI_data6/FT_analysis/tutorial/t10_tshift.txt
3dvolreg -base 0 -prefix ${tmpDir}/epi_dv.nii.gz -zpad 1 -1Dfile ${tmpOutDir}/motion.1D ${tmpDir}/epi_d.nii.gz # volume registation and extraction of motion trace


echo ""
echo "#########################################################################################################"
echo "######################################### Run field map in SPM ##########################################"
echo "#########################################################################################################"
echo ""

#gunzip ${tmpDir}/epi_dt.nii.gz
# Loop through template MATLAB script replacing keywords
#for i in ${scriptDir}'/spm_fieldmap.m'; do
#sed -e 's@SUB_SUBJECT_SUB@'$sub'@g' \
#	-e 's@SUB_NUMTRS_SUB@'$expLen'@g' \
#	-e 's@SUB_TMPOUTDIR_SUB@'$tmpOutDir'@g' \
#	-e 's@SUB_TASK_SUB@'$task'@g' <$i> spm_fieldmap.m
#done
# run script
#/usr/local/bin/matlab -nodisplay < spm_fieldmap.m
# clean up
#rm ${tmpDir}/meanuepi_dt.nii
#gzip ${tmpDir}/wfmag_epi_dt.nii
#gzip ${tmpDir}/uepi_dt.nii
#gzip ${tmpDir}/epi_dt.nii
#mv ${tmpDir}/uepi_dt.nii.gz ${tmpDir}/epi_dtv.nii.gz
#rm ${tmpDir}/epi_dt.nii
#3drefit -TR 2 -space ORIG -view orig ${tmpDir}/epi_dtv.nii.gz

3dAutomask -prefix ${tmpDir}/epi_ExtractionMask.nii.gz ${tmpDir}/epi_dv.nii.gz #Create brain mask for extraction
3dcalc -a ${tmpDir}/epi_ExtractionMask.nii.gz -b ${tmpDir}/epi_dv.nii.gz -expr 'a*b' -prefix ${tmpDir}/epi_dvb.nii.gz #extract brain
3dTstat -prefix ${tmpDir}/epi_dvbm.nii.gz ${tmpDir}/epi_dvb.nii.gz # Create mean image for more robust alignment to sub T1
N4BiasFieldCorrection -i ${tmpDir}/epi_dvbm.nii.gz #Correct image non-uniformities in mean(not analyzed just for registration) to improve coregistration

echo ""
echo "#########################################################################################################"
echo "###################Calculates Warps to align epi to template via subjects HighRes T1#####################"
echo "#########################################################################################################"
echo ""
###Only calculate for rest1 then apply to rest2 as they are volreged to same base

##Use BBR to align EPI to Subject HighRes #citation: ftp://ftp.nmr.mgh.harvard.edu/pub/articles/greve.2009.ni.63.BBR.pdf and https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT_BBR
##Comments: This method works much better in our data for a few reasons: First, a video 1 TR from each subjects aligned and warped EPI reveals much more stable and consisted alignment in BBR than rigid registration, White matter tracts are more clearer and consistently in the same location. Second tSNR maps clearly reveal this shift, old maps are shifted down with substantially worse tSNR in dorsal gray matter. Here tSNR is better in new pipeline on the order of 200%. in ventral areas old is generally better but we've decided this is because WM has higher tSNR and this signal is artifactually being pulled down into ventral GM because of poor alignment. For these reasons we are sticking to BBR and confident in combination of BBR and ants SYN. Without BBR ants SYN with rigid performed globally worse than old pipeline. This is probably because SYN just multiplies noise from poor epi alignment
3dcalc -a ${antDir}/${antPre}BrainSegmentation.nii.gz -expr 'equals(a,3)' -prefix ${tmpDir}/ExtractedBrain0N4_wmsegtmp.nii.gz
fslreorient2std ${tmpDir}/ExtractedBrain0N4_wmsegtmp.nii.gz ${tmpDir}/ExtractedBrain0N4_FSL_wmseg.nii.gz
fslreorient2std ${tmpDir}/epi_dvbm.nii.gz ${tmpDir}/epi_dvbmFSL.nii.gz
fslreorient2std ${antDir}/${antPre}BrainSegmentation0N4.nii.gz ${tmpDir}/BrainSegmentation0N4_FSL.nii.gz
robustfov -i ${tmpDir}/BrainSegmentation0N4_FSL.nii.gz
fslreorient2std ${antDir}/${antPre}ExtractedBrain0N4.nii.gz ${tmpDir}/ExtractedBrain0N4_FSL.nii.gz
robustfov -i ${tmpDir}/ExtractedBrain0N4_FSL.nii.gz

########################################### MAX START HERE ########################################################
##BBR register 
epi_reg --epi=${tmpDir}/epi_dvbmFSL.nii.gz --t1=${tmpDir}/BrainSegmentation0N4_FSL.nii.gz --t1brain=${tmpDir}/ExtractedBrain0N4_FSL.nii.gz --out=${tmpDir}/epi2highResBBR -v
c3d_affine_tool -ref ${tmpDir}/ExtractedBrain0N4_FSL.nii.gz -src ${tmpDir}/epi_dvbmFSL.nii.gz ${tmpDir}/epi2highResBBR.mat -fsl2ras -oitk ${tmpOutDir}/epi2highRes0GenericAffine.mat
#3dTstat -prefix ${tmpDir}/epi2highResBBRmean.nii.gz ${tmpDir}/epi2highResBBR.nii.gz
voxSize=$(@GetAfniRes ${tmpDir}/epi.nii.gz)
3dresample -input ${templateDir}/${templatePre}Brain.nii.gz -dxyz 2 2 2 -prefix ${tmpDir}/refTemplate4epi.nii.gz ##Citation: Decided to switch to resample to 2mm iso after testing and showing that group activation maps are more robust and significant when this step is added

##Apply Warps #citation: https://github.com/stnava/ANTs/wiki/antsCorticalThickness-and-antsLongitudinalCorticalThickness-output and https://github.com/maxwe128/restTools/blob/master/preprocessing/norm.func.spm12sa.csh 


antsApplyTransforms -d 3 -e 3 -i ${tmpDir}/epi_dvbm.nii.gz -o ${tmpDir}/epiWarpedMean.nii.gz -r ${tmpDir}/refTemplate4epi.nii.gz -t ${antDir}/${antPre}SubjectToTemplate1Warp.nii.gz -t ${antDir}/${antPre}SubjectToTemplate0GenericAffine.mat -t ${tmpOutDir}/epi2highRes0GenericAffine.mat -n Bspline #At first Bspline looked best but we switched back to Linear because of marked improvements in tSNR maps in linear compared to Bspline compared to oldSPM methods
antsApplyTransforms -d 3 -e 3 -i ${tmpDir}/epi_dvbm.nii.gz -o ${tmpDir}/epiMean2highRes.nii.gz -r ${tmpDir}/refTemplate4epi.nii.gz -t ${tmpOutDir}/epi2highRes0GenericAffine.mat -n Linear
antsApplyTransforms -d 3 -e 3 -i ${tmpDir}/epi_dv.nii.gz -o ${tmpDir}/epiWarped.nii.gz -r ${tmpDir}/refTemplate4epi.nii.gz -t ${antDir}/${antPre}SubjectToTemplate1Warp.nii.gz -t ${antDir}/${antPre}SubjectToTemplate0GenericAffine.mat -t ${tmpOutDir}/epi2highRes0GenericAffine.mat -n Bspline 
3drefit -space MNI -view tlrc ${tmpDir}/epiWarped.nii.gz #Refit space of warped epi so that it can be viewed in MNI space within AFNI

#####Smooth Data 6mm will get output to about 11-13 FWHM on average

3dTstat -prefix ${tmpDir}/mean.epiWarped.nii.gz ${tmpDir}/epiWarped.nii.gz
3dcalc -a ${tmpDir}/epiWarped.nii.gz -b ${tmpDir}/mean.epiWarped.nii.gz -expr 'min(200, a/b*100)*step(a)*step(b)' -prefix ${tmpOutDir}/epiWarped.nii.gz ##pipenotes: this is scaling all values to have comparaple beta weights across subjects. Make sure to indicate in wiki entry that the unblurred data is not scaled!!!! Also the scaling was motivated by this post #citation: https://afni.nimh.nih.gov/pub/dist/edu/data/CD.expanded/AFNI_data6/FT_analysis/tutorial/t14_scale.txt and is not being done for rest currently.  


echo ""
echo "#########################################################################################################"
echo "###############Get Motion and QA vals and make Aligment Montage for Visual Inspection####################"
echo "#########################################################################################################"
echo ""
######Calculation of motion, DVARs, alignment correlation and plots to have for censoring and QC
#Calculate FD, #citation: Power et al., 2012 
fsl_motion_outliers -i ${tmpDir}/epi_d.nii.gz -o ${tmpDir}/conf -s ${tmpOutDir}/FD_FSL.1D --fd #Use fsl tools because the automatically do it the same as Power 2012 #citation: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLMotionOutliers
fsl_motion_outliers --nomoco --dvars -m ${tmpDir}/epi_ExtractionMask.nii.gz -o ${tmpDir}/tempConf -s fslDVARS.1D -i tmp/epi_dvb.nii.gz #calc FSL style for the purpose of comparison #pipenotes: might want to remove
#Caluclate DVARS #citation: Nichols, 2013
DVARS.sh ${tmpDir}/epi_dv.nii.gz ${tmpOutDir}/DVARS.1D  ##Use Tom Nichols standardized DVARs #pipenotes: here is one paper with a standardized DVARs threshold(1.8) #citation: https://pdfs.semanticscholar.org/8b64/7293808a4903e877f93d8241428b7596a909.pdf
#Calculate derivative of motion params for confound regression, #citation: Power et al., 2014
1d_tool.py -infile ${tmpOutDir}/motion.1D -derivative -write ${tmpOutDir}/motion_deriv.1D
#calculate TRs above threshold
awk -v thresh=".25" '{if($1 > thresh) print NR}' ${tmpOutDir}/FD_FSL.1D > ${tmpOutDir}/FD.25TRs.1D #find TRs above threshold 
awk -v thresh=".5" '{if($1 > thresh) print NR}' ${tmpOutDir}/FD_FSL.1D > ${tmpOutDir}/FD.5TRs.1D #find TRs above threshold 
###Make QC vals file, spatial correlations for each warp, and anything else from http://ccpweb.wustl.edu/pdfs/2013hcp2_barch.pdf

##Make QC/QA montages
#pipenotes: consider making carpet plots and WC-RSFC plots for at least all connectivity pipelines before and after preprocessing
3dresample -input ${templateDir}/${templatePre}BrainExtractionMask.nii.gz -master ${tmpDir}/refTemplate4epi.nii.gz -prefix ${tmpDir}/refTemplateBrainMask.nii.gz
#tSNR
3dTstat -cvarinvNOD -prefix ${tmpDir}/tSNR.nii.gz ${tmpDir}/epi_dvb.nii.gz #mean/temporal sd,  #citation: same as Marcus et al., 2013 tSNR
3dTstat -cvarinvNOD -prefix ${tmpOutDir}/tSNR.EpiWarped.nii.gz ${tmpOutDir}/epiWarped.nii.gz #pipenotes: file can be used in futre to define group masks like in SPM if wanted
tSNR=$(3dBrickStat -mask ${tmpDir}/epi_ExtractionMask.nii.gz ${tmpDir}/tSNR.nii.gz | sed "s/ *//g" | sed "s/\t\t*//g" )
#smoothness before and after warping
rawFWHM=$(3dFWHMx -mask ${tmpDir}/epi_ExtractionMask.nii.gz -combine ${tmpDir}/epi_dtvb.nii.gz | cut -d " " -f11 | sed "s/ *//g" | sed "s/\t\t*//g" )
warpedFWHM=$(3dFWHMx -mask ${tmpDir}/refTemplateBrainMask.nii.gz -combine ${tmpOutDir}/epiWarped.nii.gz | cut -d " " -f11 | sed "s/ *//g" | sed "s/\t\t*//g" )
#Motion QC Measures
nVols=$(3dinfo -nv ${tmpOutDir}/epiWarped.nii.gz)
FDavg=$(1d_tool.py -show_mmms -infile ${tmpOutDir}/FD_FSL.1D | sed 's/ /,/g' | cut -d "," -f15 | sed '/^$/d')
FDsd=$(1d_tool.py -show_mmms -infile ${tmpOutDir}/FD_FSL.1D | sed 's/ /,/g' | cut -d "," -f25 | sed '/^$/d')
DVARSavg=$(1d_tool.py -show_mmms -infile ${tmpOutDir}/DVARS.1D | sed 's/ /,/g' | cut -d "," -f15 | sed '/^$/d')
DVARSsd=$(1d_tool.py -show_mmms -infile ${tmpOutDir}/DVARS.1D | sed 's/ /,/g' | cut -d "," -f25 | sed '/^$/d')
numCenFD25=$(cat ${tmpOutDir}/FD.25TRs.1D | wc -l)
numCenFD50=$(cat ${tmpOutDir}/FD.5TRs.1D | wc -l)
FD25Per=$(echo "${numCenFD25}/${nVols}" | bc -l | cut -c1-5)
FD50Per=$(echo "${numCenFD50}/${nVols}" | bc -l | cut -c1-5)
#spatial correlations between epi and highres, epi and Template, and highRes and template as a crude index of alignment quality
antsApplyTransforms -d 3 -i ${antDir}/${antPre}ExtractedBrain0N4.nii.gz -o ${tmpDir}/BrainNormalizedToTemplate.nii.gz -t ${antDir}/${antPre}SubjectToTemplate1Warp.nii.gz -t ${antDir}/${antPre}SubjectToTemplate0GenericAffine.mat -r ${antDir}/${antPre}ExtractedBrain0N4.nii.gz -n Linear
epi2TemplateCor=$(3ddot -mask ${tmpDir}/refTemplateBrainMask.nii.gz ${tmpDir}/refTemplate4epi.nii.gz ${tmpDir}/epiWarpedMean.nii.gz | sed "s/ *//g" | sed "s/\t\t*//g" )
epi2highResCor=$(3ddot -mask ${antDir}/${antPre}BrainExtractionMask.nii.gz ${antDir}/${antPre}ExtractedBrain0N4.nii.gz ${tmpDir}/epi2highResBBR.nii.gz | sed "s/ *//g"| sed "s/\t\t*//g")
highRes2TemplateCor=$(3ddot -mask ${templateDir}/${templatePre}BrainExtractionMask.nii.gz  ${tmpDir}/BrainNormalizedToTemplate.nii.gz ${templateDir}/${templatePre}Brain.nii.gz | sed "s/ *//g" | sed "s/\t\t*//g")
#set up QC table for subject
echo "tSNR,rawFWHM,warpedFWHM,FDavg,FDsd,DVARSavg,DVARSsd,FD25Per,FD50Per,epi2TemplateCor,epi2highResCor,highRes2TemplateCor" > ${QADir}/${outDirName}.QCmeasures.txt
echo "$tSNR,$rawFWHM,$warpedFWHM,$FDavg,$FDsd,$DVARSavg,$DVARSsd,$FD25Per,$FD50Per,$epi2TemplateCor,$epi2highResCor,$highRes2TemplateCor" >> ${QADir}/${outDirName}.QCmeasures.txt
######Make Edges and montage of warped edges overlain target Struct to check alignment

##epi Alignment to HighRes with BBR T1 White Matter Edges on EPI
ConvertScalarImageToRGB 3 ${tmpDir}/epi2highResBBR_fast_wmedge.nii.gz ${tmpDir}/epi2highResBBREdgesRBG.nii.gz none red none 0 10 #convert for Ants Montage
3dcalc -a ${tmpDir}/epi2highResBBREdgesRBG.nii.gz -expr 'step(a)' -prefix ${tmpDir}/epi2highResBBREdgesRBGstep.nii.gz #Make mask to make Edges stand out
CreateTiledMosaic -i ${tmpDir}/epi2highResBBR.nii.gz -r ${tmpDir}/epi2highResBBREdgesRBG.nii.gz -o ${QADir}/${outDirName}.epi2HighResBBRAlignmentCheckAxial.png -a 0.8 -t -1x-1 -d 2 -p mask -s [5,mask,mask] -x ${tmpDir}/epi2highResBBREdgesRBGstep.nii.gz -f 0x1  #Create Montage taking images in axial slices every 5 slices
CreateTiledMosaic -i ${tmpDir}/epi2highResBBR.nii.gz -r ${tmpDir}/epi2highResBBREdgesRBG.nii.gz -o ${QADir}/${outDirName}.epi2HighResBBRAlignmentCheckCoronal.png -a 0.8 -t -1x-1 -d 1 -p mask -s [5,mask,mask] -x ${tmpDir}/epi2highResBBREdgesRBGstep.nii.gz -f 0x1
##HighRes Alignment to Template
3dedge3 -input ${tmpDir}/BrainNormalizedToTemplate.nii.gz -prefix ${tmpDir}/highRes2TemplateWarpedEdges.nii.gz  #Detect edges
ConvertScalarImageToRGB 3 ${tmpDir}/highRes2TemplateWarpedEdges.nii.gz ${tmpDir}/highRes2TemplateEdgesRBG.nii.gz none red none 0 10 #convert for Ants Montage
3dcalc -a ${tmpDir}/highRes2TemplateEdgesRBG.nii.gz -expr 'step(a)' -prefix ${tmpDir}/highRes2TemplateEdgesRBGstep.nii.gz #Make mask to make Edges stand out
CreateTiledMosaic -i ${templateDir}/${templatePre}BrainSegmentation0N4.nii.gz -r ${tmpDir}/highRes2TemplateEdgesRBG.nii.gz -o ${QADir}/${outDirName}.highRes2TemplateAlignmentCheck.png -a 0.8 -t -1x-1 -d 2 -p mask -s [5,mask,mask] -x ${tmpDir}/highRes2TemplateEdgesRBGstep.nii.gz -f 0x1  #Create Montage taking images in axial slices every 5 slices
##epiWarpedMean to Template
antsApplyTransforms -d 3 -e 3 -i ${tmpDir}/epi_dvbm.nii.gz -o ${tmpDir}/epiWarpedMeanHR.nii.gz -r ${templateDir}/${templatePre}BrainSegmentation0N4.nii.gz  -t ${antDir}/${antPre}SubjectToTemplate1Warp.nii.gz -t ${antDir}/${antPre}SubjectToTemplate0GenericAffine.mat -t ${tmpOutDir}/epi2highRes0GenericAffine.mat -n Linear
3dedge3 -input ${tmpDir}/epiWarpedMeanHR.nii.gz -prefix ${tmpDir}/epi2TemplateWarpedEdges.nii.gz  #Detect edges
ConvertScalarImageToRGB 3 ${tmpDir}/epi2TemplateWarpedEdges.nii.gz ${tmpDir}/epi2TemplateEdgesRBG.nii.gz none red none 0 10 #convert for Ants Montage
3dcalc -a ${tmpDir}/epi2TemplateEdgesRBG.nii.gz -expr 'step(a)' -prefix ${tmpDir}/epi2TemplateEdgesRBGstep.nii.gz #Make mask to make Edges stand out
CreateTiledMosaic -i ${templateDir}/${templatePre}BrainSegmentation0N4.nii.gz -r ${tmpDir}/epi2TemplateEdgesRBG.nii.gz -o ${QADir}/${outDirName}.epi2TemplateAlignmentCheck.png -a 0.8 -t -1x-1 -d 2 -p mask -s [5,mask,mask] -x ${tmpDir}/epi2TemplateEdgesRBGstep.nii.gz -f 0x1  #Create Montage taking images in axial slices every 5 slices

##Montages to check B0 unwarping
#3dresample -master ${tmpDir}/epi2highResBBR.nii.gz -inset ${tmpDir}/epi_dv.nii.gz'[0]' -prefix ${tmpDir}/epiPostb0.nii.gz
#3dresample -master ${tmpDir}/epi2highResBBR.nii.gz -inset ${tmpDir}/epi_d.nii.gz'[0]' -prefix ${tmpDir}/epiPreb0.nii.gz
#CreateTiledMosaic -i ${tmpDir}/epiPreb0.nii.gz -r ${tmpDir}/epiPreb0.nii.gz -o ${QADir}/${outDirName}.preB0.png -a 0 -t -1x-1 -d 2 -p mask -s [15,0,120] -x ${tmpDir}/epiPreb0.nii.gz -f 0x1 -p 0
#CreateTiledMosaic -i ${tmpDir}/epiPostb0.nii.gz -r ${tmpDir}/epiPostb0.nii.gz -o ${QADir}/${outDirName}.postB0.png -a 0 -t -1x-1 -d 2 -p mask -s [15,0,120] -x ${tmpDir}/epiPostb0.nii.gz -f 0x1 -p 0

### copy files for vis check
#p ${QADir}/$task.epi2TemplateAlignmentCheck.png /mnt/BIAC/munin4.dhe.duke.edu/Hariri/DBIS.01/Graphics/Data_Check/NewPipeline/$task.epi2TemplateAlignmentCheck/$sub.png

##Clean up
rm -r $tmpDir

##Now copy all files to the server. 
#Because we run into issues when too many processes are trying to do this in parallel, use a lock dir system to make sure that only a small number of processes are doing this step simultaneously
N_allowed=4
lockDir=/mnt/BIAC/munin4.dhe.duke.edu/Hariri/DBIS.01/Data/ALL_DATA_TO_USE/Imaging/x_x.KEEP.OUT.x_x/locks
if [ ! -e $lockDir ]; then mkdir $lockDir; fi
break=0
while true; do
	for i in `seq 1 $N_allowed`; do
		if mkdir $lockDir/writingBigFiles$i; then
			while true; do
				if mkdir -p $outDir; then
					rsync -v --stats --progress $tmpOutDir/* $outDir # check out -W option, --timeout
					echo rsync return code: $?
					rm -r $lockDir/writingBigFiles$i
					break=1
					break
				else
					echo "mkdir $outDir failed! Sleeping 5"
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

##run first level model
#if [ $task != "rest" ]; then
#	qsub $scriptDir/glm_$task.bash $sub
#else
#	qsub $scriptDir/rest_DBIS.sh $sub
#fi

# -- BEGIN POST-USER -- 
echo "----JOB [$JOB_NAME.$JOB_ID] STOP [`date`]----" 
mv $HOME/$JOB_NAME.$JOB_ID.out $outDir/$JOB_NAME.$JOB_ID.out	 
# -- END POST-USER -- 


##########Citations
#Jo, H. J., Gotts, S. J., Reynolds, R. C., Bandettini, P. A., Martin, A., Cox, R. W., & Saad, Z. S. (2013). Effective preprocessing procedures virtually eliminate distance-dependent motion artifacts in resting state FMRI. Journal of Applied Mathematics, 2013. http://doi.org/10.1155/2013/935154
#Kalcher, K., Boubela, R. N., Huf, W., Biswal, B. B., Baldinger, P., Sailer, U., … Windischberger, C. (2013). RESCALE: Voxel-specific task-fMRI scaling using resting state fluctuation amplitude. NeuroImage, 70, 80–88. http://doi.org/10.1016/j.neuroimage.2012.12.019
#Marcus, D. S., Harms, M. P., Snyder, A. Z., Jenkinson, M., Wilson, J. A., Glasser, M. F., … Van Essen, D. C. (2013). Human Connectome Project informatics: Quality control, database services, and data visualization. NeuroImage, 80, 202–219. http://doi.org/10.1016/j.neuroimage.2013.05.077
#Nichols, 2013 http://www2.warwick.ac.uk/fac/sci/statistics/staff/academic-research/nichols/scripts/fsl/StandardizedDVARS.pdf and https://www2.warwick.ac.uk/fac/sci/statistics/staff/academic-research/nichols/scripts/fsl/DVARS.sh
#Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. NeuroImage, 59(3), 2142–2154. http://doi.org/10.1016/j.neuroimage.2011.10.018

