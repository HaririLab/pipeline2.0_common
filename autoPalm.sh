#!/bin/bash

##Auto construction of Palm commands for DBIS or DNS investigation
##Will need to be adapted for more covariates and 2 group analyses


###TO DO

#Save files that are registered to fsAverage, then just register files that haven't been previousl registered
#then just concatenate the relevant subjects and delete concatanated file after 

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
if [[ ! -f ${prefix}_tfce_tstat_mcfwep_lh_c1.gii ]];then
	rm -rf ${prefix}.goodSublist.tmp ${prefix}_design.txt ${prefix}_design.mat ${prefix}_contrasts.txt ${prefix}_contrasts.con
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
	numC=$(echo $phenoC | tr -cd , | wc -c)
	if [[ $numC == 0 ]];then
		echo "1 covariate of interest.....Runing a brainWide ANCOVA"
		printf "0 1\n0 -1\n" > ${prefix}_contrasts.txt
	elif [[ $numC == "1" ]];then
		echo "2 covariates of interest....Runing a two-sample t-test"
		printf "0 1 -1\n0 -1 1\n" > ${prefix}_contrasts.txt
		printf "0 1 -1\n0 -1 1\n"
		phenoCol2=$(echo $phenoC | cut -d "," -f2)
	else
		echo " more than 2 covariates of intersest....This script can't handle...exiting"
		exit
	fi
	echo "Subjects included in analysis:"
	lhList=""
	rhList=""
	for i in $(cat $subList);do
		#echo $i
		if [[ $numC == 0 ]];then
			pheno=$(grep $i $phenoFile | cut -d "," -f${phenoCol})
			img=${SUBJECTS_DIR}/${i}/surf/lh.$imageType
			sex=$(grep $i $phenoFile | cut -d "," -f2)
			if [[ $pheno != "NA" ]] && [[ -f $img ]];then
				##Make TCFE files		
				echo "$sex ${pheno}" >> ${prefix}_design.txt
				echo $i >> ${prefix}.goodSublist.tmp
				lhList=$(echo "$lhList ${SUBJECTS_DIR}/${i}/surf/lh.${imageType}.fsAverageReg.mgh")
				rhList=$(echo "$rhList ${SUBJECTS_DIR}/${i}/surf/rh.${imageType}.fsAverageReg.mgh")
				echo -n " $i"
			fi
		else
			pheno=$(grep $i $phenoFile | cut -d "," -f${phenoCol})
			pheno2tmp=$(grep $i $phenoFile | cut -d "," -f${phenoCol2})
			pheno2=$(echo " $pheno2tmp")
			img=${SUBJECTS_DIR}/${i}/surf/lh.$imageType
			sex=$(grep $i $phenoFile | cut -d "," -f2)
			if [[ $pheno != "NA" ]] && [[ -f $img ]];then
				##Make TCFE files		
				echo "$sex ${pheno}${pheno2}" >> ${prefix}_design.txt
				echo $i >> ${prefix}.goodSublist.tmp
				lhList=$(echo "$lhList ${SUBJECTS_DIR}/${i}/surf/lh.${imageType}.fsAverageReg.mgh")
				rhList=$(echo "$rhList ${SUBJECTS_DIR}/${i}/surf/rh.${imageType}.fsAverageReg.mgh")
				echo -n " $i"
			fi
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
		#register all images to fsAverage, save in allImaging dir, speeds up process so you don't reregister everytime
		for i in $(cat ${prefix}.goodSublist.tmp);do
			if [[ ! -f ${SUBJECTS_DIR}/${i}/surf/rh.$imageType.fsAverageReg.mgh ]] && [[ ! -f ${SUBJECTS_DIR}/${i}/surf/lh.$imageType.fsAverageReg.mgh ]];then
				mri_surf2surf --srcsubject $i --srchemi lh --srcsurfreg sphere.reg --trgsubject fsaverage --trghemi lh --trgsurfreg sphere.reg --tval ${SUBJECTS_DIR}/${i}/surf/lh.area.fsAverageReg.mgh --sval ${SUBJECTS_DIR}/${i}/surf/lh.area --jac --sfmt curv --noreshape --cortex
				mri_surf2surf --srcsubject $i --srchemi rh --srcsurfreg sphere.reg --trgsubject fsaverage --trghemi rh --trgsurfreg sphere.reg --tval ${SUBJECTS_DIR}/${i}/surf/rh.area.fsAverageReg.mgh --sval ${SUBJECTS_DIR}/${i}/surf/rh.area --jac --sfmt curv --noreshape --cortex
			else
				echo "$i already registered"
			fi
		done
		##Concatanate all the relevant files 
		mri_concat --o lh.all.area.00.mgh --i $lhList
		mri_concat --o rh.all.area.00.mgh --i $rhList
		tm_tools vertex-box-cox-transform -i lh.all.area.00.mgh $threads
		tm_tools vertex-box-cox-transform -i rh.all.area.00.mgh $threads
		mri_surf2surf --hemi lh --sval lh.all.area.00.boxcox.mgh --s fsaverage --fwhm $smooth --cortex --tval ${prefix}_tmp.lh.all.area.sm${smooth}.boxcox.mgh #convert to gii for compatability with SUMA
		mri_surf2surf --hemi rh --sval rh.all.area.00.boxcox.mgh --s fsaverage --fwhm $smooth --cortex --tval ${prefix}_tmp.rh.all.area.sm${smooth}.boxcox.mgh
		rm lh.all.${imageType}.00.mgh rh.all.${imageType}.00.mgh lh.all.${imageType}.00.boxcox.mgh rh.all.${imageType}.00.boxcox.mgh
		# Only correct over modalities becuause contrasts are orthogonal, ie both testing the same thing just looking for opposite sign.
		# Only detrend in ancova, messes up 2 group analyses gives an error about diminished rank
		if [[ $numC == 0 ]];then
			palm -i ${prefix}_tmp.lh.all.area.sm${smooth}.boxcox.mgh -i ${prefix}_tmp.rh.all.area.sm${smooth}.boxcox.mgh -d ${prefix}_design.mat -t ${prefix}_contrasts.con -o ${prefix} -n 1000 -T -tfce2d -s ${SUBJECTS_DIR}/fsaverage/surf/lh.white ${SUBJECTS_DIR}/fsaverage/surf/lh.white.avg.area.mgh -s ${SUBJECTS_DIR}/fsaverage/surf/rh.white ${SUBJECTS_DIR}/fsaverage/surf/rh.white.avg.area.mgh -nouncorrected -demean -corrmod -corrcon -demean -logp
		else
			palm -i ${prefix}_tmp.lh.all.area.sm${smooth}.boxcox.mgh -i ${prefix}_tmp.rh.all.area.sm${smooth}.boxcox.mgh -d ${prefix}_design.mat -t ${prefix}_contrasts.con -o ${prefix} -n 1000 -T -tfce2d -s ${SUBJECTS_DIR}/fsaverage/surf/lh.white ${SUBJECTS_DIR}/fsaverage/surf/lh.white.avg.area.mgh -s ${SUBJECTS_DIR}/fsaverage/surf/rh.white ${SUBJECTS_DIR}/fsaverage/surf/rh.white.avg.area.mgh -nouncorrected -corrmod -corrcon -logp
		fi
		####Use Surf2surf to convert all output files into gii for reading with SUMA
	elif [[ $imageType == "thickness" ]];then
		for i in $(cat ${prefix}.goodSublist.tmp);do
			if [[ ! -f ${SUBJECTS_DIR}/${i}/surf/rh.$imageType.fsAverageReg.mgh ]] && [[ ! -f ${SUBJECTS_DIR}/${i}/surf/lh.$imageType.fsAverageReg.mgh ]];then
				mri_surf2surf --srcsubject $i --srchemi lh --srcsurfreg sphere.reg --trgsubject fsaverage --trghemi lh --trgsurfreg sphere.reg --tval ${SUBJECTS_DIR}/${i}/surf/lh.thickness.fsAverageReg.mgh --sval ${SUBJECTS_DIR}/${i}/surf/lh.thickness --sfmt curv --noreshape --cortex
				mri_surf2surf --srcsubject $i --srchemi rh --srcsurfreg sphere.reg --trgsubject fsaverage --trghemi rh --trgsurfreg sphere.reg --tval ${SUBJECTS_DIR}/${i}/surf/rh.thickness.fsAverageReg.mgh --sval ${SUBJECTS_DIR}/${i}/surf/rh.thickness --sfmt curv --noreshape --cortex
			else
				echo "$i already registered"
			fi
		done
		##Concatanate all the relevant files 
		mri_concat --o ${prefix}_tmp.lh.all.thickness.00.mgh --i $lhList
		mri_concat --o ${prefix}_tmp.rh.all.thickness.00.mgh --i $rhList
		mri_surf2surf --hemi lh --sval ${prefix}_tmp.lh.all.thickness.00.mgh --s fsaverage --fwhm $smooth --cortex --tval ${prefix}_tmp.lh.all.thickness.sm${smooth}.mgh 
		mri_surf2surf --hemi rh --sval ${prefix}_tmp.rh.all.thickness.00.mgh --s fsaverage --fwhm $smooth --cortex --tval ${prefix}_tmp.rh.all.thickness.sm${smooth}.mgh
		if [[ $numC == 0 ]];then
			palm -i ${prefix}_tmp.lh.all.${imageType}.sm$smooth.mgh -i ${prefix}_tmp.rh.all.${imageType}.sm$smooth.mgh -d ${prefix}_design.mat -t ${prefix}_contrasts.con -o ${prefix} -n 1000 -T -tfce2d -s ${SUBJECTS_DIR}/fsaverage/surf/lh.white ${SUBJECTS_DIR}/fsaverage/surf/lh.white.avg.area.mgh -s ${SUBJECTS_DIR}/fsaverage/surf/rh.white ${SUBJECTS_DIR}/fsaverage/surf/rh.white.avg.area.mgh -nouncorrected -corrcon -corrmod -demean -logp
		else
			palm -i ${prefix}_tmp.lh.all.${imageType}.sm$smooth.mgh -i ${prefix}_tmp.rh.all.${imageType}.sm$smooth.mgh -d ${prefix}_design.mat -t ${prefix}_contrasts.con -o ${prefix} -n 1000 -T -tfce2d -s ${SUBJECTS_DIR}/fsaverage/surf/lh.white ${SUBJECTS_DIR}/fsaverage/surf/lh.white.avg.area.mgh -s ${SUBJECTS_DIR}/fsaverage/surf/rh.white ${SUBJECTS_DIR}/fsaverage/surf/rh.white.avg.area.mgh -nouncorrected -corrcon -corrmod -logp
		fi
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
	rm ${prefix}_allImgs.nii.gz
	rm ${prefix}*.tmp
	rm ${prefix}_tmp*
else
	echo " $prefix already completed"
fi
