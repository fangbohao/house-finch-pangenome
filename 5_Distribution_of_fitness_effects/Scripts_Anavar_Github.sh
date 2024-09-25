#!/bin/bash
#
#SBATCH -p shared # partition (queue)
#SBATCH -c 1 # number of cores
#SBATCH --cpus-per-task=20
#SBATCH --mem 90000 # memory pool for all cores
#SBATCH -t 3-00:00 # time (D-HH:MM)

#SBATCH -o Anavar_Jan18_2024.out # STDOUT
#SBATCH -e Anavar_Jan18_2024.err # STDERR
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bfang@fas.harvard.edu

cd Pangenome/Anavar

$anavar All_selINDEL_CDS_neuINDEL.control All_selINDEL_CDS_neuINDEL.out All_selINDEL_CDS_neuINDEL.log & \
$anavar All_selINDEL_intron_neuINDEL.control All_selINDEL_intron_neuINDEL.out All_selINDEL_intron_neuINDEL.log & \
$anavar All_selSNP_CDS_neuSNP.control All_selSNP_CDS_neuSNP.out All_selSNP_CDS_neuSNP.log & \
$anavar All_selSNP_intron_neuSNP.control All_selSNP_intron_neuSNP.out All_selSNP_intron_neuSNP.log & \
$anavar All_selSV_CDS_neuINDEL.control All_selSV_CDS_neuINDEL.out All_selSV_CDS_neuINDEL.log & \
$anavar All_selSV_intron_neuINDEL.control All_selSV_intron_neuINDEL.out All_selSV_intron_neuINDEL.log & \
$anavar East_selINDEL_CDS_neuINDEL.control East_selINDEL_CDS_neuINDEL.out East_selINDEL_CDS_neuINDEL.log & \
$anavar East_selINDEL_intron_neuINDEL.control East_selINDEL_intron_neuINDEL.out East_selINDEL_intron_neuINDEL.log & \
$anavar East_selSNP_CDS_neuSNP.control East_selSNP_CDS_neuSNP.out East_selSNP_CDS_neuSNP.log & \
$anavar East_selSNP_intron_neuSNP.control East_selSNP_intron_neuSNP.out East_selSNP_intron_neuSNP.log & \
$anavar East_selSV_CDS_neuINDEL.control East_selSV_CDS_neuINDEL.out East_selSV_CDS_neuINDEL.log & \
$anavar East_selSV_intron_neuINDEL.control East_selSV_intron_neuINDEL.out East_selSV_intron_neuINDEL.log & \
$anavar West_selINDEL_CDS_neuINDEL.control West_selINDEL_CDS_neuINDEL.out West_selINDEL_CDS_neuINDEL.log & \
$anavar West_selINDEL_intron_neuINDEL.control West_selINDEL_intron_neuINDEL.out West_selINDEL_intron_neuINDEL.log & \
$anavar West_selSNP_CDS_neuSNP.control West_selSNP_CDS_neuSNP.out West_selSNP_CDS_neuSNP.log & \
$anavar West_selSNP_intron_neuSNP.control West_selSNP_intron_neuSNP.out West_selSNP_intron_neuSNP.log & \
$anavar West_selSV_CDS_neuINDEL.control West_selSV_CDS_neuINDEL.out West_selSV_CDS_neuINDEL.log & \
$anavar West_selSV_intron_neuINDEL.control West_selSV_intron_neuINDEL.out West_selSV_intron_neuINDEL.log

wait