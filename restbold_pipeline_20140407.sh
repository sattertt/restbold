#!/bin/bash
# ---------------------------------------------------------------
# restbold_single_subject_20140301.sh
#
# Run basic resting bold subject-level processing.
# Features:
#	-Converts dicoms niftis if desired. 
#		Dependency: sequence2nifti.sh
# 	-Performs distortion correction on timeseries if desired
#		Dependency: dico_correct_v2.sh
#	-Runs fsl-style prestats
#		Dependency: fsl5
#	-Coregister with T1-- FLRIT w/ or w/o BBR
#		Dependency: fsl5
#	-Creates design matricies consistent with a variety of confound regression techniques
#		Dependencies: R, make_conf_design.R, ANTS
#	-Runs simultanous confound regression and bandpass filtering
#		Dependency: AFNI
#	-Registers images to standard space-- using FLIRT, FNIRT, ANTS, or DRAMMS
#		Depeneencies: fsl5, ANTS, DRAMMS
#	-Extracts networks timeseries & calculates network measures
#		Dependencies: fsl, R, iGraph
#	-Makes seed correlation maps
#		Dependencies: AFNI
#
# Created: Ted Satterthwaite 
# 	    3/2014
# Contact: sattertt@upennn.edu 
#
# Please cite: Satterthwaite, T.D., Elliott, M.A., Gerraty, R.T., Ruparel, K., Loughead, J., Calkins, M.E., Eickhoff, S.B., Hakonarson, H., Gur, R.C., Gur, R.E., and Wolf, D.H. (2013). An improved framework for confound regression and filtering for control of motion artifact in the preprocessing of resting-state functional connectivity data. Neuroimage 64, 240-256.
#
# ---------------------------------------------------------------

#USAGE ---------------------------------------------------------------
Usage() {
    echo ""
    echo "Usage: `basename $0` [options] --subj=<subjid>  --dir=<output directory> --prestats_design=<prestats design.fsf file> --t1brain=<brain extracted t1> --t1wm=<T1 WM segment> --t1csf=<t1 CSF segment>"
    echo ""
    echo "Input arguments"
    echo "  --subj=<subjid>			: subjid (e.g. BBLID_SCNAID, ID_DATE) that gets appended to output files"
    echo "  --dir=<outdir>			: output directory.  Within this there are three subdirectories-- prestats, confound_regress, and coregistration"
    echo "  --prestats_design=<deisgn.fsf>      : design file for FSL prestats. This is to maintain backwards-compatibility; may be removed in future."
    echo "  --t1brain=<image>			: brain extracted T1 image"
    echo "  --t1seg=<image>			: Hard segmentation of T1 image"
    echo "  --t1seg_vals=<CSF, GM, WM>		: Intensity value in T1 image of CSF, GM, and WM segments"
    echo "  --coreg_method=<method>		: coregistration method; valid options at present are cost functions supported by flirt; default is BBR"
    echo "  --smooth=<sigma>			: smoothing in sigma mm (i.e. 6mm FWHM=2.54); applied in standard space"
    echo "  --config=<config>			: config file the specifies dramms, ants, R scripts directories"
    echo ""
    echo "Required input options"
    echo "  --input_nifti=<timeseries>		: input nifti timeseries"
    echo "**OR**"
    echo "  --input_dicoms=<dicom_dir>		: directory to find dicoms-- **BBL ONLY AT PRESENT**"
    echo ""
    echo "Normalization options-- at least one warp and template image is required"
    echo "   --ants_warp=<warp>			: use ANTs deformation"
    echo "   --ants_affine=<affine>             : path to ANTs affine"
    echo "   --ants_rigid=<rigid>               : path to ANTs rigid"
    echo "   --dramms=<warp>		 	: use DRAMMs deformation"
    echo "   --template=<image>                 : standard template that images will be registered to"
    echo "   --downsample=<mat>			: rigid body matrix for downsampling-- used only for example if normalization warp was to 1mm and desired output is 2mm"	
	
    echo ""
    echo "Optional distortion correction arguments--- **BBL ONLY AT PRESENT**, note unwrap sign is set to negative!"
    echo "   -d					: perform distortion correction"
    echo "   --example_dicom=<dcm>		: example dicom of restbold timeseries"
    echo "   --mag_map=<image>			: magnitude map for coregistraiton, from dico_b0_calc"
    echo "   --rps_map=<image>			: rps b0 map from dico_b0_calc"
    echo "   --b0_mask=<image>			: b0 mask file from dico_b0_calc"
    echo ""
    echo "Optional output"
    echo "   --seed=<image>			: run seed correlation on all unique labeled ROIs in specified image"
    echo "   --network=<image>			: extract timeseries from unique labled ROIs in specified image"
    echo "   -G					: extract graphical network measurs from network using iGraph-- **PENDING**"
    echo "   -vmhc=<warp>			: register images from template to specified asym temmplate using warp provided here-- **PENDING**"
    echo "   -vmhc=<sym template>		: symmetric template for VMHC-- **PENDING**"
    echo ""
    echo "" 
    echo "   -h		                        : display this help message"
    echo ""
    exit 1
}
# ---------------------------------------------------------------


# ---------------------------------------------------------------
# Functions for argument parsing
# ---------------------------------------------------------------
get_opt1() {
    arg=`echo $1 | sed 's/=.*//'`
    echo $arg
}

get_imarg1() {
    arg=`get_arg1 $1`;
    arg=`$FSLDIR/bin/remove_ext $arg`;
    echo $arg
}

get_arg1() {
    if [ X`echo $1 | grep '='` = X ] ; then
	echo "Option $1 requires an argument" 1>&2
	exit 1
    else
	arg=`echo $1 | sed 's/.*=//'`
	if [ X$arg = X ] ; then
	    echo "Option $1 requires an argument" 1>&2
	    exit 1
	fi
	echo $arg
    fi
}


# ---------------------------------------------------------------
# Make sure have minimal arguments before continuing (6 total)
# ---------------------------------------------------------------
if [ $# -lt 6 ] ; then Usage; exit 0; fi


# ---------------------------------------------------------------
# Set defaults
# ---------------------------------------------------------------

subj=""
dir=""
prestats_deisgn=""
t1brain=""
t1seg=""
t1seg_vals=""
coreg_method=bbr
use_nifti=no
use_dicoms=no
ants_warp=""
ants_rigid=""
ants_affine=""
dramms=""
template=""
downsample=""
roi=""
smooth=""
network=""
unwarp_sign=""
use_dico_nifti=no
example_dicom=""
rps_map=""
mag_map=""
b0_mask=""

# ---------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------

while [ $# -ge 1 ] ; do
    iarg=`get_opt1 $1`;
    case "$iarg"
	in
	--subj)
		subj=$(get_arg1 $1)
		shift;;
	--dir)
		dir=$(get_arg1 $1)
		shift;;
	--prestats_design)
		prestats_design=$(get_arg1 $1)        
		shift;;
	--t1brain)
		t1brain=$(get_arg1 $1)
		shift;;
	--t1seg)
		t1seg=$(get_imarg1 $1)
		shift;;
        --t1seg_vals)
                t1seg_vals=$(get_imarg1 $1)
                shift;;
	--coreg_method)
                coreg_method=$(get_arg1 $1)
                shift;;
	--smooth)
                smooth=$(get_arg1 $1)
                shift;;
	--config)
		config=$(get_arg1 $1)
		shift;;
	--input_nifti)
		input_nifti=$(get_imarg1 $1)
		use_nifti=yes
		shift;;
	--input_dicoms)
                dicom_dir=$(get_arg1 $1)
                use_dicoms=yes
		shift;;
	--ants_warp)
		ants_warp=$(get_arg1 $1)
		use_ants=yes
		shift;;
        --ants_affine)
                ants_affine=$(get_arg1 $1)
                shift;;
	--ants_rigid)
                ants_rigid=$(get_arg1 $1)
                shift;;
	--dramms)
		dramms=$(get_arg1 $1)
                use_dramms=yes
                shift;;
	--template)
		template=$(get_arg1 $1)
		shift;;
	--downsample)
		downsample=$(get_arg1 $1)
		use_downsample=yes
                shift;;
	--seed)
		roi=$(get_arg1 $1)
		use_seed=yes
		shift;;
	--network)
                network=$(get_arg1 $1)
		use_network=yes
		shift;;
	-d)	
		use_dico=yes
		shift;;
	--example_dicom)
		use_dico_nifti=yes
		example_dicom=$(get_arg1 $1)
		shift;;
	--mag_map)
		mag_map=$(get_arg1 $1)
		shift;;
	--rps_map)
		rps_map=$(get_arg1 $1)
		shift;;
	--b0_mask)
		b0_mask=$(get_arg1 $1)
		shift;;
	-h)
	  	Usage;
		exit 0;;
	*)
		echo "Unrecognised option $1" 1>&2
		exit 1;;
	esac
done

#----------------------------------------------------------------
# Check required input arguments
# ---------------------------------------------------------------

if [ X$subj = X ] ; then
  echo "The compulsory argument --subj to specify subject ID MUST be used"
  exit 1;
fi

if [ "$use_nifti" = no ] && [ "$use_dicoms" =  no ]; then 
	echo "nifti = $use_nifti"
	echo "dicom = $use_dicoms"
	echo "Must specify EITHER dicoms or nifti input timeseries"
	exit 1
fi

if [ X$dir = X ] ; then
  echo "The compulsory argument --dir to specify output directory MUST be used"
  exit 1;
fi

if [ X$prestats_design = X ] ; then
  echo "The compulsory argument --prestats_design to specify the prestats design .fsf template MUST be used"
  exit 1;
else
  prestats_design=$(ls $prestats_design 2> /dev/null)
  if [ ! -e "$prestats_design" ]; then
    echo "prestats design not found!"
    exit 1
  fi
fi


if [ X$config = X ] ; then
  echo "The compulsory argument --config to specify a configuration file MUST be specified"
  exit 1;
fi

if [ X$t1brain = X ] ; then
  echo "The compulsory argument --t1brain to specify a brain extracted t1 image MUST be specfied"
  exit 1;
else
  t1brain_test=$(imtest $t1brain)
  if [ "$t1brain_test" -eq 0 ]; then
    echo "t1brain image not found!"
    exit 1
  fi
fi


if [ X$t1seg = X ] ; then
  echo "The compulsory argument --t1seg to specify hard segmentation of t1 image MUST be specfied"
  exit 1;
else
  t1seg_test=$(imtest $t1seg)
  if [ "$t1seg_test" -eq 0 ]; then
    echo "t1 seg image not found!"
    exit 1
  fi
fi

if [ X$t1seg_vals = X ] ; then
  echo "The compulsory argument --t1seg_vals to specify intensity vals of t1 segments (GM, WM, CSF) MUST be used"
  exit 1;
fi

if [ "$use_ants" = no ] && [ "$use_dramms" =  no ]; then
        echo "Must specify EITHER a dramms or ants warp"
        exit 1
fi


if [ "$use_dramms" = yes ]; then
        dramms_warp=$(ls $dramms)
        if [ ! -e "$dramms_warp" ]; then
                echo "dramms warp not found!"
                exit 1
        fi
fi

if [ "$use_ants" = yes ]; then
	if [ ! -e "$ants_warp" ] || [ ! -e "$ants_affine" ] || [ ! -e "$ants_rigid" ]; then
		"at least one ants file is missing!"
		exit 1
	fi
fi


if [ X$template = X ] ; then
  echo "The compulsory argument --template to specify template image MUST be specfied"
  exit 1;
else
  template_test=$(imtest $template)
  if [ "$template_test" -eq 0 ]; then
    echo "template image not found!"
    exit 1
  fi
fi


if [ "$use_seed" = yes ]; then
        if [ ! -e "$roi" ]; then
                "roi file is missing"
                exit 1
        fi
	if [ X$smooth = X ] ; then
		  echo "The compulsory argument --smooth to specify standard-space smoothing MUST be used if a seed analysis is run"
	  exit 1;
	fi
fi

if [ "$use_network" = yes ]; then
        if [ ! -e "$network" ]; then
                echo "network file is missing!"
                exit 1
	fi
fi

if [ "$use_dico" = yes ]; then
	if [ ! -e "$mag_map" ] || [ ! -e "$rps_map" ] || [ ! -e $b0_mask ]; then
		echo "magnitude map, rpsmap or b0mask missing- need to specify these to run distortion correction"
		exit 1
	fi
	if [ "$use_dicoms" = no ] && [ "$use_dico_nifti" = no ]; then
		echo "need to either specify input dicoms or an example dicom to perform distortion correction"
		exit 1
	fi
	if [ "$use_nifti" = yes ] && [ "$use_dico_nifti" = no ]; then
                echo "need to provide an example dicom if input data is nifti timeseries"
                exit 1
        fi
        if [ "$use_nifti" = yes ] && [ "$use_dico_nifti" = no ]; then
                echo "need to provide an example dicom if input data is nifti timeseries"
                exit 1
        fi

fi
	
    echo ""

#----------------------------------------------------------------
# Display input arguments
# ---------------------------------------------------------------
#specify scriptdir
echo ""
echo "__________________________________________"
echo "Input arguments are:"
echo "subject id is $subj"
if [ "$use_nifti" == yes ]; then
	echo "input timeseries is $input_nifti"
fi
if [ "$use_dicoms" == yes ]; then
        echo "input dicom directory is $dicom_dir"
fi
echo "output directory is $dir"
echo "config file is $config"
echo "prestats design file is $prestats_design"
echo "t1 brain is $t1brain"
echo "t1seg is $t1seg"
echo "t1 segmentation values are $t1seg_vals"
echo "coreg method is $coreg_method"
echo "script directory is $scriptdir"
if [ "$use_dramms" = yes ]; then
        echo "dramms warp is $dramms_warp"
fi
if [ "$use_ants" = yes ]; then
	echo "ants warp is $ants_warp"
	echo "ants_affine is $ants_affine"
	echo "ants_rigid is $ants_rigid"
fi
if [ "$use_downsample" = yes ]; then
	echo "downsample is $downsample"
fi
if [ "$use_seed" = yes ]; then
	echo "seed roi is $roi"
fi
if [ "$use_network" = yes ]; then
        echo "network image is $network"
fi
if [ "$use_dico" = yes ]; then
        echo "will perform distortion correction"
	echo "magnitude image for coregistration is $mag_map"
	echo "field map is $rps_map"
	echo "mask is $b0_mask"
fi
 
echo "__________________________________________"

#----------------------------------------------------------------
# Setup paths-- from config file
# ---------------------------------------------------------------
source $config


if [ "$use_dramms" = yes ]; then
        echo "dramms directory is $DRAMMSDIR"
fi

if [ "$use_ants" = yes ]; then
        echo "ants directory is $ANTSDIR"
fi

echo "R dir is $RDIR"
echo "R scripts directory is $RSCRIPTDIR"
echo "FSL directory is $FSLDIR"


#prevent prestats from auto-submitting to grid
export SGE_ROOT=""


#Make output directory
if [  ! -d "$dir" ]; then
	echo "making output directory $dir"
	mkdir $dir
	dir=$(ls -d $dir)
fi

#----------------------------------------------------------------
# Run distortion correction
# ---------------------------------------------------------------

if [ "$use_dico" = yes ]; then
	echo ""
	echo "__________________________________________"

	#check if output is present
        if [ -e "$dir/dico/${subj}_restbold_dico.nii.gz" ]; then
                echo "dico already done!"
        else
		echo "set to run distoriton correction"
	        echo "Note that will assume negative sign!!"
		#make dico directory
		if [ ! -d "$dir/dico" ]; then
        		echo "making dico directory"
                	mkdir $dir/dico
	        fi
		dicodir=$dir/dico
		echo "output director is  $dicodir"
		cd $dicodir #may be necessary for the output directory to work?

        	#first check if input is nifti or example dicom
		if [ "$use_dico_nifti" = yes ]; then
			echo "will distortion correct data from specified nifti, using the exmple dicom"
			$DICOSCRIPTSDIR/dico_correct_v2.sh -n -FS -e $example_dicom -f $mag_map $dicodir/${subj}_restbold $rps_map $b0_mask $input_nifti
		fi

		if [ "$use_dico_nifti" = yes ]; then
                	echo "will distortion correct data from specified nifti, using the example dicom"
	                $DICOSCRIPTSDIR/dico_correct_v2.sh -n -FS -e $example_dicom -f $mag_map $dicodir/${subj}_restbold $rps_map $b0_mask $input_nifti
        	fi

	        if [ "$use_dico_nifti" = no ]; then
        	        echo "will distortion correct data from dicoms within directory $dicom_dir"
			dicoms=$(ls $dicom_dir/*.dcm)
	                $DICOSCRIPTSDIR/dico_correct_v2.sh -n -FS -f $mag_map ${subj}_restbold $rps_map $b0_mask $dicoms
        	fi

		#convert back to .nii.gz format
		echo "converting timeseries back to zipped nifti"
		fslchfiletype NIFTI_GZ $dicodir/${subj}_restbold
        	fslchfiletype NIFTI_GZ $dicodir/${subj}_restbold_dico

		#check that output was made correctly
		if [ ! -e "$dicodir/${subj}_restbold_dico.nii.gz" ]; then
			echo "expected distortion-corrected timeseries not present!"
			exit 1
		fi

	fi
        
	#set variables for pipeline steps	
        input_nifti=$dir/dico/${subj}_restbold_dico.nii.gz
        use_dicoms=no  #set this so do not try to convet dicoms to nifti w/o dico in next step
fi

#----------------------------------------------------------------
# Convert dicoms to nifti as needed -- BBL option only at present
# ---------------------------------------------------------------

if [ "$use_dicoms" == yes ]; then
	input_nifti=$(ls $dir/nifti/${subj}_restbold.nii.gz)
	if  [ ! -e "$input_nifti" ]; then
		echo "converting dicoms to nifti"
		if [ ! -d "$dir/nifti" ]; then
			echo "making nifti directory"
			mkdir $dir/nifti	
		fi
		$DICOSCRIPTSDIR/sequence2nifti.sh BOLD $dir/nifti/${subj}_restbold.nii $dicom_dir/*.dcm   #this is hard coded
		echo "changing file type"
		fslchfiletype NIFTI_GZ $dir/nifti/${subj}_restbold.nii
		input_nifti=$(ls $dir/nifti/${subj}_restbold.nii.gz)
		if [ ! -e "$input_nifti" ]; then
			echo "expected nifti not present! something went wrong!"
			exit 1
		fi
	else
		echo "output nifti in fact already present! will not convert dicoms to nifti!"
	fi
fi


#----------------------------------------------------------------
# Run prestats
# ---------------------------------------------------------------

#run prestats-- this can probably be improved by taking it off self-submisison to qsub
echo ""
echo "__________________________________________"
prestats_name=$(basename $prestats_design | cut -d. -f1)
prestats_out=$(ls $dir/${prestats_name}.feat/${subj}_filtered_func_data.nii.gz 2> /dev/null)
if [ ! -e "$prestats_out" ]; then
	echo "running prestats"
	prestats_name=$(basename $prestats_design | cut -d. -f1)
	###note that this convention requires output name design file name to be congrunet
	nvol=$(fslinfo $input_nifti |grep dim4 |head -n1 | awk '{print $2}') #calculate number of volx
	echo "input nifti image is $input_nifti"
	echo "output directory is $dir/${prestats_name}.feat"
        echo "nvol is $nvol"
	cp $prestats_design $dir/"$subj"_prestats.fsf
        sed -i "s%XINPUT_NIFTIX%$input_nifti%g" $dir/"$subj"_prestats.fsf
	sed -i "s%XDIRX%$dir%g" $dir/"$subj"_prestats.fsf
        sed -i s/XNVOLX/$nvol/g $dir/"$subj"_prestats.fsf
	feat $dir/"$subj"_prestats.fsf

	k=0
	while [ ! -e $dir/${prestats_name}.feat/absbrainthresh.txt -a $k -lt 31 ]; do
		echo "***** prestats still running"
		k=$((k+1))
		sleep 120
	done

	echo "***** prestats has finished running"

	if [ $k == 31 ];then
		echo "feat prestats timed out!!"
		exit 1
	fi
	
	#rename output images to include subject id
	if [ ! -e "$dir/${prestats_name}.feat/${subj}_mean_func.nii.gz" ]; then
		echo "renaming images"
		prestats_imgs=$(ls ${dir}/${prestats_name}.feat/*.nii.gz)
		for i in $prestats_imgs; do 
			imgname=$(basename $i)
			mv $i $dir/${prestats_name}.feat/${subj}_${imgname}
		done
	fi

else 
	echo "prestats already run"
fi

#define prestats outputs that will be needed later
prestatsdir=$(ls -d ${dir}/${prestats_name}.feat)
example_func=$(ls $prestatsdir/${subj}_example_func.nii.gz)
filtered_func=$(ls $prestatsdir/${subj}_filtered_func_data.nii.gz)
mask=$(ls $prestatsdir/${subj}_mask.nii.gz)
if [ ! -e "$filtered_func" ]; then
	echo "expected prestats output not present-- something went wrong"
	exit 1
fi



#----------------------------------------------------------------
# Coregister functional and structural data
# ---------------------------------------------------------------

echo ""
echo "__________________________________________"


#check if output is present
output=$(ls $dir/coregistration/${subj}_ep2struct.mat 2> /dev/null)
if [ -e "$output" ]; then
	echo "coregistration already run and complete"
else
        #make output directory if not already present
        if [ ! -d "$dir/coregistration" ]; then
                echo "making coregistration directory"
                mkdir ${dir}/coregistration
        fi
        coregdir=$(ls -d ${dir}/coregistration)

	#bet the example_func image if not done
	if [ ! -e "$prestatsdir/${subj}_example_func_brain.nii.gz" ]; then
		echo "running bet on example_func"
		bet $example_func $prestatsdir/${subj}_example_func_brain -f 0.3
	fi
         example_func_brain=$(ls $prestatsdir/${subj}_example_func_brain.nii.gz)

        #make structural segments if not done yet
	if [ ! -e "$coregdir/${subj}_t1wm.nii.gz" ]; then
		echo "making structural WM and CSF segments"
	        csfval=$(echo $t1seg_vals | cut -d, -f1)
	        wmval=$(echo $t1seg_vals | cut -d, -f3)
	        echo "csf value is $csfval, wm val is $wmval"
		fslmaths $t1seg -thr $csfval -uthr $csfval -bin $coregdir/${subj}_t1csf
                fslmaths $t1seg -thr $wmval -uthr $wmval -bin $coregdir/${subj}_t1wm
	fi
	t1wm=$(ls $coregdir/${subj}_t1wm.nii.gz)
	t1csf=$(ls $coregdir/${subj}_t1csf.nii.gz)


	#run flirt
	echo "running coregistration using method $coreg_method"
	echo "in $example_func_brain"
	echo "ref $t1brain"
	echo "out $coregdir/${subj}_ep2struct"
	echo "omat $coregdir/${subj}_ep2struct.mat"
	flirt -in $example_func_brain -ref $t1brain -dof 6 -out $coregdir/${subj}_ep2struct -omat $coregdir/${subj}_ep2struct.mat -cost $coreg_method -wmseg $t1wm 	 
        convert_xfm -omat $coregdir/${subj}_struct2ep.mat -inverse $coregdir/${subj}_ep2struct.mat

	cp $example_func_brain $coregdir  #copy to coregdir for help viewing subject-space rois
	cp $mask $coregdir
fi

#general coregistration-resultant dependencies for later use
coregdir=$(ls -d $dir/coregistration)
coregmat=$(ls $coregdir/${subj}_ep2struct.mat)
coregmat_inv=$(ls $coregdir/${subj}_struct2ep.mat)
t1wm=$(ls $coregdir/${subj}_t1wm.nii.gz)
t1csf=$(ls $coregdir/${subj}_t1csf.nii.gz)
example_func_brain=$(ls $prestatsdir/${subj}_example_func_brain.nii.gz)
if [ ! -e "$coregmat" ]; then
	echo "coregistration not present as expected-- something went wrong!!"
	exit 1
fi


#----------------------------------------------------------------
# Get confound timecourses
# ---------------------------------------------------------------

echo ""
echo "__________________________________________"

#check if output is present
confound_output=$(ls $prestatsdir/confound_timecourses/${subj}_mask_timecourse.txt 2> /dev/null)
if [ -e "$confound_output" ]; then
        echo "confound timecourses already extracted"
else
        echo "extracting confound timecourses"

	#make output directory if not present
	if [ ! -d "$prestatsdir/confound_timecourses" ]; then
		mkdir $prestatsdir/confound_timecourses
	        tcdir=$(ls -d $prestatsdir/confound_timecourses)
	fi
        tcdir=$(ls -d $prestatsdir/confound_timecourses)

	#move t1wm and t1csf to epi image
        flirt -in $t1wm -out $coregdir/${subj}_t1wm_2ep -ref $example_func_brain -applyxfm -init $coregmat_inv
        flirt -in $t1csf -out $coregdir/${subj}_t1csf_2ep -ref $example_func_brain -applyxfm -init $coregmat_inv

	#threshold/binarize
	echo "thresholding/binarizing wm/csf segments"
	fslmaths $coregdir/${subj}_t1wm_2ep -thr 0.5 -bin $coregdir/${subj}_t1wm_2ep
        fslmaths $coregdir/${subj}_t1csf_2ep -thr 0.5 -bin $coregdir/${subj}_t1csf_2ep

	#extract timecourses
	echo "now extracting timecourses"
	fslmeants -i $filtered_func -m $coregdir/${subj}_t1wm_2ep -o $tcdir/${subj}_wm_timecourse.txt
        fslmeants -i $filtered_func -m $coregdir/${subj}_t1csf_2ep -o $tcdir/${subj}_csf_timecourse.txt
        fslmeants -i $filtered_func -m $coregdir/${subj}_mask -o $tcdir/${subj}_mask_timecourse.txt
fi

#----------------------------------------------------------------
# Make confound matrix
# ---------------------------------------------------------------

echo ""
echo "__________________________________________"

confmat_output=$(ls ${prestatsdir}/confound_regress_36EV/confmat36EV.txt 2> /dev/null)
if [ -e "$confmat_output" ]; then
        echo "confound matrix already presnt"
else
        echo "making 24 and 36 EV confound matricies"
	$RDIR/R  --slave --file=${RSCRIPTDIR}/design_matrix_36EV_20140122.R --args "$prestatsdir"  #note hard coded
fi

if [ ! -e "${prestatsdir}/confound_regress_36EV/confmat36EV.txt" ]; then
	echo "confound matrix not present as expected-- exiting"
	exit 1
fi

#----------------------------------------------------------------
# Run confound regression
# ---------------------------------------------------------------

echo ""
echo "__________________________________________"

#check if output is present
confreg_output=$(ls ${prestatsdir}/confound_regress_36EV/${subj}_filtered_func_data_36EV.nii.gz 2> /dev/null)
if [ -e "$confreg_output" ]; then
        echo "regressed timeseries already present"
else
        echo "running 24 and 36 EV confound regressions"
	echo ""
	echo "36EV"
	3dBandpass -prefix ${prestatsdir}/confound_regress_36EV/${subj}_filtered_func_data_36EV.nii.gz -ort ${prestatsdir}/confound_regress_36EV/confmat36EV.txt 0.01 0.08 $filtered_func 
	echo ""
	echo "24EV"
	3dBandpass -prefix ${prestatsdir}/confound_regress_24EV/${subj}_filtered_func_data_24EV.nii.gz -ort ${prestatsdir}/confound_regress_24EV/confmat24EV.txt 0.01 0.08 $filtered_func
fi

filtered_func_36ev=$(ls ${prestatsdir}/confound_regress_36EV/${subj}_filtered_func_data_36EV.nii.gz)
filtered_func_24ev=$(ls ${prestatsdir}/confound_regress_24EV/${subj}_filtered_func_data_24EV.nii.gz)
if [ ! -e "$filtered_func_36ev" ] || [ ! -e "$filtered_func_24ev" ]; then
	echo "confound regression did not complete as expected"
	exit 1
fi


#----------------------------------------------------------------
# Move regressed timecourses to standard space
# ---------------------------------------------------------------

echo ""
echo "__________________________________________"

#check if output present
if [ ! -e "${prestatsdir}/confound_regress_36EV/${subj}_filtered_func_data_36EV_std.nii.gz" ]; then
	echo "moving timecourses to standard space"
	
	#Check which deformation-- ANTS or DRAMMS-- is specified
	if [ "$use_dramms" = yes ]; then 
		echo "using dramms"

		#combine warps
		echo "combining warps"
		$DRAMMSDIR/dramms-combine -c -f $example_func_brain -t $t1brain $coregmat $dramms_warp $coregdir/ep2mni_warp.nii.gz

		echo "36ev"
		$DRAMMSDIR/dramms-warp $filtered_func_36ev $coregdir/ep2mni_warp.nii.gz ${prestatsdir}/confound_regress_36EV/${subj}_filtered_func_data_36EV_std.nii.gz

		echo "24ev"
		$DRAMMSDIR/dramms-warp $filtered_func_36ev $coregdir/ep2mni_warp.nii.gz ${prestatsdir}/confound_regress_24EV/${subj}_filtered_func_data_24EV_std.nii.gz

		echo "moving example func and mask to standard space"
                $DRAMMSDIR/dramms-warp $mask $coregdir/ep2mni_warp.nii.gz ${prestatsdir}/${subj}_mask_std.nii.gz
                $DRAMMSDIR/dramms-warp $example_func_brain $coregdir/ep2mni_warp.nii.gz ${prestatsdir}/${subj}_example_func_brain_std.nii.gz


	fi

        if [ "$use_ants" = yes ]; then
		echo "using ants"
	
		#convert the coregistration to ANTS format 		
                coregtxt=$(ls -d $coregdir/${subj}_ep2struct.txt 2> /dev/null)
		if [ ! -e "$coregtxt" ]; then
			echo "converting mat to ants format"
			c3d_affine_tool -src $example_func_brain -ref $t1brain $coregmat -fsl2ras -oitk $coregdir/${subj}_ep2struct.txt
	                coregtxt=$(ls -d $coregdir/${subj}_ep2struct.txt)
		fi
		
		if [ "$use_downsample" = no ]; then
			${ANTSDIR}/antsApplyTransforms -e 3 -d 3 -i $filtered_func_24ev -o ${prestatsdir}/confound_regress_24EV/${subj}_filtered_func_data_24EV_std.nii.gz  -r $template  -t $ants_warp -t $ants_affine -t $ants_rigid -t $coregtxt 
	
			${ANTSDIR}/antsApplyTransforms -e 3 -d 3 -i $filtered_func_36ev -o ${prestatsdir}/confound_regress_36EV/${subj}_filtered_func_data_36EV_std.nii.gz  -r $template  -t $ants_warp -t $ants_affine -t $ants_rigid -t $coregtxt 

			${ANTSDIR}/antsApplyTransforms -e 3 -d 3 -i $mask -o ${prestatsdir}/${subj}_mask_std.nii.gz  -r $template  -t $ants_warp -t $ants_affine -t $ants_rigid -t $coregtxt
			${ANTSDIR}/antsApplyTransforms -e 3 -d 3 -i $example_func_brain -o ${prestatsdir}/${subj}_example_func_brain_std.nii.gz  -r $template  -t $ants_warp -t $ants_affine -t $ants_rigid -t $coregtxt

		fi

                if [ "$use_downsample" = yes ]; then
                        ${ANTSDIR}/antsApplyTransforms -e 3 -d 3 -i $filtered_func_24ev -o ${prestatsdir}/confound_regress_24EV/${subj}_filtered_func_data_24EV_std.nii.gz  -r $template  -t $downsample -t $ants_warp -t $ants_affine -t $ants_rigid -t $coregtxt 

                        ${ANTSDIR}/antsApplyTransforms -e 3 -d 3 -i $filtered_func_36ev -o ${prestatsdir}/confound_regress_36EV/${subj}_filtered_func_data_36EV_std.nii.gz  -r $template  -t $downsample -t $ants_warp -t $ants_affine -t $ants_rigid -t $coregtxt 

			${ANTSDIR}/antsApplyTransforms -e 3 -d 3 -i $mask -o ${prestatsdir}/${subj}_mask_std.nii.gz  -r $template -t $downsample -t $ants_warp -t $ants_affine -t $ants_rigid -t $coregtxt

                        ${ANTSDIR}/antsApplyTransforms -e 3 -d 3 -i $example_func_brain -o ${prestatsdir}/${subj}_example_func_brain_std.nii.gz  -r $template -t $downsample -t $ants_warp -t $ants_affine -t $ants_rigid -t $coregtxt

                fi

	fi
else
	echo "timecourses already in standard space"
fi

data24=${prestatsdir}/confound_regress_24EV/${subj}_filtered_func_data_24EV_std.nii.gz
data36=${prestatsdir}/confound_regress_36EV/${subj}_filtered_func_data_36EV_std.nii.gz

if [ ! -e "$data24" ] || [ ! -e "$data36" ]; then
	echo "standard space timecourses not present as expected"
	exit 1
fi

#----------------------------------------------------------------
# Seed correlation map
# ---------------------------------------------------------------

echo ""
echo "__________________________________________"

if [ "$use_seed" == yes ]; then
	echo "set to run seed correlation analysis"

	#smooth first
	if [ ! -e "${prestatsdir}/confound_regress_36EV/${subj}_filtered_func_data_36EV_std_sigma${smooth}.nii.gz" ]; then
		echo "smoothing data at sigma $smooth"
		echo "24ev"
		fslmaths $data24 -s $smooth ${prestatsdir}/confound_regress_24EV/${subj}_filtered_func_data_24EV_std_sigma${smooth}
		echo "36ev"
	        fslmaths $data36 -s $smooth ${prestatsdir}/confound_regress_36EV/${subj}_filtered_func_data_36EV_std_sigma${smooth}
	else
		echo "data already smoothed at sigma $smooth"
	fi

	data24sm=${prestatsdir}/confound_regress_24EV/${subj}_filtered_func_data_24EV_std_sigma${smooth}.nii.gz
	data36sm=${prestatsdir}/confound_regress_36EV/${subj}_filtered_func_data_36EV_std_sigma${smooth}.nii.gz


	#now get seed TS
	roiname=$(basename $roi | cut -d. -f1)
	if [ ! -e "${prestatsdir}/confound_regress_36EV/roi_timecourses/${subj}_${roiname}_tc_sigma${smooth}.txt" ]; then
		echo "extracting TS"
		if [ ! -d "$prestatsdir/confound_regress_36EV/roi_timecourses" ]; then
			echo "making tc output directories"
			mkdir $prestatsdir/confound_regress_24EV/roi_timecourses
			mkdir $prestatsdir/confound_regress_36EV/roi_timecourses
		fi
		fslmeants -i $data24sm -m $roi -o ${prestatsdir}/confound_regress_24EV/roi_timecourses/${subj}_${roiname}_tc_sigma${smooth}.txt
        	fslmeants -i $data36sm -m $roi -o ${prestatsdir}/confound_regress_36EV/roi_timecourses/${subj}_${roiname}_tc_sigma${smooth}.txt
	else
		echo "TS present"
	fi

	#run and transform correlations

	if [ ! -e "${prestatsdir}/confound_regress_36EV/seed_maps/" ]; then
		echo "making seed map directories"
		mkdir ${prestatsdir}/confound_regress_36EV/seed_maps/
		mkdir ${prestatsdir}/confound_regress_24EV/seed_maps/
	fi

	zrout=$(ls -d ${prestatsdir}/confound_regress_36EV/seed_maps/${subj}_${roiname}_sigma${smooth}_zr.nii.gz 2> /dev/null)
	if [ ! -e "$zrout" ]; then
		echo ""
		echo "Computing Correlations"
		3dfim+ -input $data24sm  -ideal_file ${prestatsdir}/confound_regress_24EV/roi_timecourses/${subj}_${roiname}_tc_sigma${smooth}.txt -out Correlation -bucket  ${prestatsdir}/confound_regress_24EV/seed_maps/${subj}_${roiname}_sigma${smooth}_r.nii.gz
        	3dfim+ -input $data36sm  -ideal_file ${prestatsdir}/confound_regress_36EV/roi_timecourses/${subj}_${roiname}_tc_sigma${smooth}.txt -out Correlation -bucket  ${prestatsdir}/confound_regress_36EV/seed_maps/${subj}_${roiname}_sigma${smooth}_r.nii.gz

		echo "Z-transforming correlations"
		3dcalc -a  ${prestatsdir}/confound_regress_24EV/seed_maps/${subj}_${roiname}_sigma${smooth}_r.nii.gz -expr 'log((a+1)/(a-1))/2' -prefix  ${prestatsdir}/confound_regress_24EV/seed_maps/${subj}_${roiname}_sigma${smooth}_zr.nii.gz
	        3dcalc -a  ${prestatsdir}/confound_regress_36EV/seed_maps/${subj}_${roiname}_sigma${smooth}_r.nii.gz -expr 'log((a+1)/(a-1))/2' -prefix  ${prestatsdir}/confound_regress_36EV/seed_maps/${subj}_${roiname}_sigma${smooth}_zr.nii.gz
	else
		echo "maps already present!"
	fi
fi

#----------------------------------------------------------------
# Network timecourses
# ---------------------------------------------------------------

echo ""
echo "__________________________________________"

if [ "$use_network" == yes ]; then
        echo "set to extract network timecourses"

	outname=$(remove_ext $network)
	outname=$(basename $outname)

	echo "network name is $outname"	

	if [ ! -e "${prestatsdir}/confound_regress_36EV/networks/" ]; then
                echo "making network directories"
                mkdir ${prestatsdir}/confound_regress_36EV/networks/
                mkdir ${prestatsdir}/confound_regress_24EV/networks/
        fi

	#find out how many rois/clusters in image
	netmax=$(fslstats $network -R  | awk '{print $2}' | cut -d. -f1)
	netmin=$(fslstats $network -l 0 -R  | awk '{print $1}' | cut -d. -f1)  #this prevents it from looking at voxels w value zero
	echo "extracting timecourses for ROIs $netmin to $netmax"

	echo $prestatsdir
	#remove old files in case incomplete
	rm -f ${prestatsdir}/confound_regress_36EV/networks//${outname}.txt
        rm -f ${prestatsdir}/confound_regress_24EV/networks//${outname}.txt

	for (( r=$netmin; r<=$netmax; r++)); do
		echo "roi $r"
		#make temprorary mask
		fslmaths $network -thr $r -uthr $r -bin $prestatsdir/masktmp
		fslmeants -i $data24  -m $prestatsdir/masktmp -o $prestatsdir/tctmpA.txt #note that are not using smoothed image
                fslmeants -i $data36  -m $prestatsdir/masktmp -o $prestatsdir/tctmpB.txt 
		if [ $r == $netmin ];then
			mv $prestatsdir/tctmpA.txt ${prestatsdir}/confound_regress_24EV/networks/${outname}_timecourse.txt
                        mv $prestatsdir/tctmpB.txt ${prestatsdir}/confound_regress_36EV/networks/${outname}_timecourse.txt

		else
			paste ${prestatsdir}/confound_regress_24EV/networks/${outname}_timecourse.txt  $prestatsdir/tctmpA.txt >>  $prestatsdir/tctmp2A.txt
                        paste ${prestatsdir}/confound_regress_36EV/networks/${outname}_timecourse.txt  $prestatsdir/tctmpB.txt >>  $prestatsdir/tctmp2B.txt

			mv $prestatsdir/tctmp2A.txt ${prestatsdir}/confound_regress_24EV/networks/${outname}_timecourse.txt
			mv $prestatsdir/tctmp2B.txt ${prestatsdir}/confound_regress_36EV/networks/${outname}_timecourse.txt
		fi
		rm -f $prestatsdir/masktmp.nii.gz
		rm -f $prestatdir/tctmp.txt
	done
fi

