#!/bin/env bash

#SBATCH --job-name=fmriprep
#SBATCH --account=kg98
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=linden.parkes@monash.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=8000
#SBATCH --qos=normal
#SBATCH -A kg98
#SBATCH --array=1-171%30
# IMPORTANT! set the array range above to exactly the number of people in your SubjectIDs.txt file. e.g., if you have 90 subjects then array should be: --array=1-90

# Assign input args
WhichProject=CNP

case $WhichProject in
	# input a text file here that lists the subjects in your bids directory to be processed.
	# Note, the 'WhichProject' stuff can all be deleted if your only processing a single dataset.
	CNP) SUBJECT_LIST="/home/lindenmp/kg98_scratch/Linden/ResProjects/NormModelling/data/CNP/CNP_SubjectIDs.txt" ;;
esac

subject=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SUBJECT_LIST})
echo -e "\t\t\t --------------------------- "
echo -e "\t\t\t ----- ${SLURM_ARRAY_TASK_ID} ${subject} ----- "
echo -e "\t\t\t --------------------------- \n"

# paths
workdir=/tmp/linden2/work/${subject} # this directory is just a temp directory where intermediate outputs/files from fmriprep are dumped. It must be subject specific though!
bidsdir=/home/lindenmp/kg98_scratch/Linden/ResProjects/NormModelling/data/CNP # path to a valid BIDS dataset
derivsdir=/home/lindenmp/kg98_scratch/Linden/ResProjects/GSR/data/CNP/derivatives # path to whether derivatives will go. This can be anywhere
fslicense=/home/lindenmp/license.txt # path and freesurfer licence .txt file. Download your own from freesurfer website and store. Or just leave it and use mine!

# other things to consider below:
# 1) -t flag. I have it set to rest but you might want to change this depending on your needs.
# 2) --use-syn-dyc is on. This is a map-free distortion correction method that is currently listed as experimental by fmriprep developers.
# 3) this variant of the code runs freesurfer.

# --------------------------------------------------------------------------------------------------
# MASSIVE modules
module load fmriprep/1.1.1
unset PYTHONPATH

# Run fmriprep
fmriprep \
--participant-label ${subject} \
--output-space {T1w,template,fsnative,fsaverage5} \
--use-aroma \
--mem_mb 80000 \
-t rest \
-n-cpus 8 \
--fs-license-file ${fslicense} \
--use-syn-sdc \
-w ${workdir} \
${bidsdir} \
${derivsdir} \
participant
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# perform non-aggressive ICA-AROMA on unsmoothed MNI/T1w fmriprep outputs
module purge
module load fsl/5.0.9

derivsdir_subj=${derivsdir}/fmriprep/${subject}/func/

if [ -f "${derivsdir_subj}${subject}_task-rest_bold_space-MNI152NLin2009cAsym_preproc.nii.gz" ]; then
    echo "AROMA denoising: space-MNI152NLin2009cAsym_preproc"
    cd ${derivsdir_subj}
    # clean MNI
    fsl_regfilt -i ${subject}_task-rest_bold_space-MNI152NLin2009cAsym_preproc.nii.gz \
        -f $(cat ${subject}_task-rest_bold_AROMAnoiseICs.csv) \
        -d ${subject}_task-rest_bold_MELODICmix.tsv \
        -o ${subject}_task-rest_bold_space-MNI152NLin2009cAsym_variant-AROMAnonaggr_preproc.nii.gz
fi

if [ -f "${derivsdir_subj}${subject}_task-rest_bold_space-T1w_preproc.nii.gz" ]; then
    echo "AROMA denoising: space-T1w_preproc"
    cd ${derivsdir_subj}
    # clean T1w
    fsl_regfilt -i ${subject}_task-rest_bold_space-T1w_preproc.nii.gz \
        -f $(cat ${subject}_task-rest_bold_AROMAnoiseICs.csv) \
        -d ${subject}_task-rest_bold_MELODICmix.tsv \
        -o ${subject}_task-rest_bold_space-T1w_variant-AROMAnonaggr_preproc.nii.gz
fi
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
if [ -f "${derivsdir_subj}${subject}_task-rest_bold_space-MNI152NLin2009cAsym_preproc.nii.gz" ] || [ -f "${derivsdir_subj}${subject}_task-rest_bold_space-T1w_preproc.nii.gz" ]; then
    cd ${derivsdir_subj}
    # perform non-aggressive ICA-AROMA+2P on unsmoothed MNI/T1w fmriprep outputs 
    # first, get wm/csf columns out of confounds file
    awk -v OFS='\t' '{if(NR>1)print $1,$2}' ${subject}_task-rest_bold_confounds.tsv > wmcsf.tsv
    # append to melodic mix
    paste ${subject}_task-rest_bold_MELODICmix.tsv wmcsf.tsv > conf.tsv
    mv conf.tsv conf.tsv.in
    sed -e $'s/\t\t/\t/g' conf.tsv.in > conf.tsv # fix double tab issue following paste
    rm conf.tsv.in
    # build custom confound files and indices
    numCol=$(awk '{print NF}' ${subject}_task-rest_bold_MELODICmix.tsv | sort -nu | tail -n 1) # number of MELODIC columns
    expr $numCol + 1 > csf_idx.txt # idx of appended csf column
    expr $numCol + 2 > wm_idx.txt # idx of appended wm column
    paste -d , csf_idx.txt wm_idx.txt > idx.csv # idx of appended csf/wm columns
    rm csf_idx.txt wm_idx.txt
    paste -d , ${subject}_task-rest_bold_AROMAnoiseICs.csv idx.csv > conf_idx.csv # append csf/wm indices to AROMA noise indices
fi

if [ -f "${derivsdir_subj}${subject}_task-rest_bold_space-MNI152NLin2009cAsym_preproc.nii.gz" ]; then
    echo "AROMA+2P denoising: space-MNI152NLin2009cAsym_preproc"
    cd ${derivsdir_subj}
    # clean MNI
    fsl_regfilt -i ${subject}_task-rest_bold_space-MNI152NLin2009cAsym_preproc.nii.gz \
        -f $(cat conf_idx.csv) \
        -d conf.tsv \
        -o ${subject}_task-rest_bold_space-MNI152NLin2009cAsym_variant-AROMAnonaggr+2P_preproc.nii.gz
fi

if [ -f "${derivsdir_subj}${subject}_task-rest_bold_space-T1w_preproc.nii.gz" ]; then
    echo "AROMA+2P denoising: space-T1w_preproc"
    cd ${derivsdir_subj}
    # clean T1w
    fsl_regfilt -i ${subject}_task-rest_bold_space-T1w_preproc.nii.gz \
        -f $(cat conf_idx.csv) \
        -d conf.tsv \
        -o ${subject}_task-rest_bold_space-T1w_variant-AROMAnonaggr+2P_preproc.nii.gz
fi
# --------------------------------------------------------------------------------------------------

echo -e "\t\t\t ----- DONE ----- "
