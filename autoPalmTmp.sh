#!/bin/bash

##Auto construction of Palm commands for DBIS or DNS investigation
##Will need to be adapted for more covariates and 2 group analyses


###TO DO


subList=$1 #set to "all" if you want to include all subs with imaging data and a non-NA pheno value
#assume data is in /mnt/BIAC/munin2.dhe.duke.edu/Hariri/DBIS.01/Analysis/All_Imaging/FreeSurfer_AllSubs
imageType=$2 #For now this will always either be "thickness" or "area"
phenoFile=$3 #location of phenotype file, has to have full subject name DMHDSXXXX
phenoC=$4 #column in pheno file that contains main variable of interest, will always covary sex and assume first col is ID and second is sex. Can give two columns ex: 7,8 for group contrasts like comparing EMH to everyone else
smooth=$5 #$FWHM to smooth Data at
threads=$6 #number of parallel threads, speeds up preprocessing significantly
prefix=$7 #include img type, pheno contrast, and date. Ex: /mnt/BIAC/munin2.dhe.duke.edu/Hariri/DBIS.01/Analysis/Max/EMH/results/cortThick_arteriolCaliber_071717
export PATH=$PATH:/munin/DNS.01/Analysis/Max/scripts/huginBin/palm-alpha106/:
studyPre=$(tail -n1 $phenoFile | cut -d "," -f1 | cut -c1-3)
if [[ $studyPre == "DNS" ]];then
	study="DNS.01"
elif [[ $studyPre == "DMH" ]];then
	study="DBIS.01"
else
	echo "Subject names are flawed and do not start with DMHDS or DNS"
	exit
fi

if [[ $imageType == "area" ]];then
	echo "area has to be run locally, not on BIAC"
	export SUBJECTS_DIR=/munin/${study}/Analysis/All_Imaging/FreeSurfer_AllSubs/
else
	export SUBJECTS_DIR=/mnt/BIAC/munin2.dhe.duke.edu/Hariri/${study}/Analysis/All_Imaging/FreeSurfer_AllSubs/
fi
##only include sunjects without NA in pheno and with processed imaging
##eventually add a check for "good" imaging data from spenser QCs

echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "Assuming first column of pheno file is ID and second is Sex"
echo "Change something if that is not true"
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

if [[ $subList = "all" ]];then
	cut -d "," -f1 $phenoFile | tail -n +2 > ${prefix}.sublist.tmp
	subList=${prefix}.sublist.tmp
fi
sed ':a;N;$!ba;s/\n/ /g' $subList > ${prefix}.fsSubList.tmp
List=""
phenoCol=$(echo $phenoC | cut -d "," -f1)
phenoCol2=$(echo $phenoC | cut -d "," -f2)
numC=$(echo $phenoC | tr -cd , | wc -c)
if [[ $numC == 0 ]];then
	echo "1 covariate of interest.....Runing a brainWide ANCOVA"
	printf "0 1\n0 -1\n" > ${prefix}_contrasts.txt
elif [[ $numC == "1" ]];then
	echo "2 covariates of interest....Runing a two-sample t-test"
	printf "0 1 -1\n0 -1 1\n" > ${prefix}_contrasts.txt
else
	echo " more than 2 covariates of intersest....This script can't handle...exiting"
	exit
fi
echo "Subjects included in analysis:"
for i in $(cat $subList);do
	#echo $i
	pheno=$(grep $i $phenoFile | cut -d "," -f${phenoCol})
	pheno2tmp=$(grep $i $phenoFile | cut -d "," -f${phenoCol2})
	pheno2=$(echo " $pheno2tmp")
	img=${SUBJECTS_DIR}/${i}/surf/lh.$imageType
	sex=$(grep $i $phenoFile | cut -d "," -f2)
	if [[ $pheno != "NA" ]] && [[ -f $img ]];then
		##Make TCFE files		
		echo "$sex ${pheno}${pheno2}" >> ${prefix}_design.txt
		echo $i >> ${prefix}.goodSublist.tmp
		echo -n " $i"
	fi
done
echo ""
##Change sex back to 1 and 0
sed -i 's/female /0 /g' ${prefix}_design.txt
sed -i 's/male /1 /g' ${prefix}_design.txt
Text2Vest ${prefix}_design.txt ${prefix}_design.mat
Text2Vest ${prefix}_contrasts.txt ${prefix}_contrasts.con

###Prepare Suface data
outDir=$(echo $prefix | rev | cut -d "/" -f2- | rev)
cd $outDir

if [[ $imageType == "area" ]];then
	#register all subjects data to fsaverage and concatenate
	#citation: https://github.com/KirstieJane/NSPN_CODE/wiki/Group-Analysis-with-Freesurfer
	tfce_mediation step0-vertex -i ${prefix}.goodSublist.tmp $imageType -p $threads #https://github.com/trislett/TFCE_mediation
	tm_tools vertex-box-cox-transform -i lh.all.area.00.mgh $threads
	tm_tools vertex-box-cox-transform -i rh.all.area.00.mgh $threads
	mri_surf2surf --hemi lh --sval lh.all.area.00.boxcox.mgh --s fsaverage --fwhm $smooth --cortex --tval lh.all.area.sm10.boxcox.mgh #convert to gii for compatability with SUMA
	mri_surf2surf --hemi rh --sval rh.all.area.00.boxcox.mgh --s fsaverage --fwhm $smooth --cortex --tval rh.all.area.sm10.boxcox.mgh
	palm -i lh.all.area.sm10.boxcox.mgh -i rh.all.area.sm10.boxcox.mgh -d ${prefix}_design.mat -t ${prefix}_contrasts.con -o ${prefix} -n 1000 -T -tfce2d -s ${SUBJECTS_DIR}/fsaverage/surf/lh.white ${SUBJECTS_DIR}/fsaverage/surf/lh.white.avg.area.mgh -s ${SUBJECTS_DIR}/fsaverage/surf/rh.white ${SUBJECTS_DIR}/fsaverage/surf/rh.white.avg.area.mgh -nouncorrected -demean -corrcon -corrmod
	####Use Surf2surf to convert all output files into gii for reading with SUMA
	rm lh.all.area.00.mgh rh.all.area.00.mgh lh.all.area.03B.mgh rh.all.area.03B.mgh lh.all.area.00.boxcox.mgh lh.all.area.boxcox.03B.mgh rh.all.area.00.boxcox.mgh rh.all.area.boxcox.03B.mgh lh.all.area.sm$smooth.boxcox.gii lh.all.area.sm$smooth.boxcox.gii
elif [[ $imageType == "thickness" ]];then 
	mris_preproc --f ${prefix}.goodSublist.tmp --target fsaverage --hemi lh --meas $imageType --fwhm $smooth --out ${prefix}_tmp.lh.allData${imageType}.sm$smooth.mgh
	mris_preproc --f ${prefix}.goodSublist.tmp --target fsaverage --hemi rh --meas $imageType --fwhm $smooth --out ${prefix}_tmp.rh.allData${imageType}.sm$smooth.mgh
	palm -i ${prefix}_tmp.lh.allData${imageType}.sm$smooth.mgh -i ${prefix}_tmp.rh.allData${imageType}.sm$smooth.mgh -d ${prefix}_design.mat -t ${prefix}_contrasts.con -o ${prefix} -n 1000 -T -tfce2d -s ${SUBJECTS_DIR}/fsaverage/surf/lh.white ${SUBJECTS_DIR}/fsaverage/surf/lh.white.avg.area.mgh -s ${SUBJECTS_DIR}/fsaverage/surf/rh.white ${SUBJECTS_DIR}/fsaverage/surf/rh.white.avg.area.mgh -nouncorrected -demean -corrcon -corrmod
	#tfce_mediation step0-vertex -i $goodSubList $imageType -p $threads -f $smooth
	#mri_surf2surf --hemi lh --sval pFactor_ThicknessPalm6_tfce_tstat_m1_c1.mgz --s fsaverage 
else
	echo "bad image type, must be area or thickness"

fi


##Convert PALM output to gii for SUMA compatability
for i in $(ls ${prefix}_tfce*.mgz);do
	out=$(echo $i | sed 's/.mgh/.gii/g' | sed 's/.mgz/.gii/g' | sed 's/m1/lh/g' | sed 's/m2/rh/g')
	mri_surf2surf --hemi lh --sval $i --s fsaverage --fwhm 0 --cortex --tval $out
done

#remove DPV files, don't think I'll use them
rm ${prefix}_*.mgz
#rm ${prefix}_allImgs.nii.gz
#rm ${prefix}*.tmp
#rm ${prefix}_tmp*
