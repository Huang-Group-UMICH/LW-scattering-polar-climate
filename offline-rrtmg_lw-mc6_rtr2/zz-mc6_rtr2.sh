#!/bin/bash 
#===================================
#  Bourne-Again shell script at flux @ U of Michigan
#
#  Description:
#
#
#  Usage:
#
#
#  History:
#
#  Author:
#    Yi-Hsuan Chen, CLaSP, U of Michigan
#    yihsuan@umich.edu
#===================================

###################
# user setting
###################

# homepath
homepath="/home/yihsuan/"

# wroking directory
wrkdir="./"

# temp variable
temp=`date +%Y%m%d%H%M%S`

#choice_file=("sas" "saw")  # sas: sub-Arctic summer, saw: sub-Arctic winter
choice_file=("sas")
#choice_file=("saw")
k_cloud=(42 46 48)  # 42=360hPa=8km, 46=637hPa=4km, 48=832hPa=1.km
#k_cloud=(46)  # 42=360hPa=8km, 46=637hPa=4km, 48=832hPa=1.km
iciwp=("0.02" "0.05" "0.1") # ice water path (kg/m2). 100 g/m2 is used in ICRCCM and Kuo et al. (2017).
#iciwp=("0.02")
#flag_scat=("Scat" "ScatEmis_ice")
flag_scat=("Scat" "noScat") # Scat and noScat cases
#flag_scat=("ScatEmis_ice")
#q_factor=("0.")
rei=("5." "10." "20." "50." "100.")  # ice effective radius (microns)
#rei=("5." "100.")  # microns
q_factor=("0.25" "0.5" "1.0" "1.25" "1.5") # scale water vapor amount. q_factor * Q
#q_factor=("1.0")
#filehead="xx3_all_cld_atm_surf"
filehead="data_plot01-xx3_all_cld_atm_surf"  # output files will start by the $filehead

fileout_head="data_plot01-bbxx3_all_cld_atm_surf"    # file head of merged files 
out_suffix="out"

filename_param="./zz-input-parameter.txt"

exe="./test111"  # name of executable

############################
# output files description
############################

#-------------------
#  Direct outputs
#-------------------
#  $filename starts by "$filehead-atm_{choice_file}-kcld_${k_cloud}-rei_{rei}-iciwp_{iciwp}-qfac_{q_factor}-{flag_scat}-"
#    $filename-flux.txt		: T & Q profiles, and all-sky & clear-sky LW fluxes and heating rates at each level will be in this file.
#    $filename-spectral.txt 	: Same as $filename-flux.txt, with additional all-sky band-by-band fluxes.
#    $filename-vars1.txt 	: Band-by-band, surface upward LW flux (FLUS), surface downward LW flux (FLDS), and TOA upward LW flux (FLUT)
#    $filename-out.txt 		: Selected output: location of the cloud layer, rei, iciwp, total column water vapor, band-avg cloud optical depth, surface downward flux, TOA upward flux
#					play(1,k_cloud),rei(1,k_cloud-1),iciwp(1,k_cloud),tmq,taucavg,dflx(1,1),uflx(1,nlay+1)

#-------------------
#  Merged outputs
#-------------------
#  merge $filename-out.txt that has the same rei & q_factor into a new file. 
#  The new file name is "$fileout_head-atm_{choice_file}-kcld_${k_cloud}-iciwp_{iciwp}-rei_q-{flag_scat}-.txt"
#  These files are used to plot Figure 1.


##################
# program start
##################

#n1=${#choice_file[@]}
#n2=${#k_cloud[@]}
#n3=${#rei[@]}
#n4=${#iciwp[@]}
#n5=${#flag_scat[@]}
#n6=${#q_factor[@]}

make

for n1 in ${choice_file[@]} ; do
for n2 in ${k_cloud[@]} ; do
for n3 in ${rei[@]} ; do
for n4 in ${iciwp[@]} ; do
for n5 in ${flag_scat[@]} ; do
for n6 in ${q_factor[@]} ; do

  #echo $n1,$n2,$n3,$n4,$n5,$n6
  n33=`printf "%3.3i" $n3`
  n66=`printf "%0.2f" $n6`

  filename="$filehead-atm_$n1-kcld_$n2-rei_$n33-iciwp_$n4-qfac_$n66-$n5"
  #echo $filename

  filename_flux=$filename-flux.txt
  filename_spectral=$filename-spectral.txt
  filename_out=$filename-out.txt
  filename_vars1=$filename-vars1.txt

  #echo $n1 $n2 $n3 $n4 $n5 $n6 $filename_flux $filename_spectral $filename_out > $filename_param \
  #  && echo "Done. create [$filename]"  || exit 1

  cat > $filename_param << EOF1
$n1
$n2
$n3
$n4
$n5
$n6
$filename_flux
$filename_spectral
$filename_out
$filename_vars1
EOF1
  echo "Done. create [$filename]"  || exit 1

  $exe 

done; done ; done ; done ; done ; done

# merge files
for n1 in ${choice_file[@]} ; do
for n2 in ${k_cloud[@]} ; do
for n4 in ${iciwp[@]} ; do
for n50 in ${flag_scat[@]} ; do
   n5="-"$n50"-"
   #ls $filehead*$n1*$n2*$n5*${out_suffix}.txt
 
   filename1="${filehead}*atm_${n1}*kcld_${n2}*iciwp_${n4}*$n5*${out_suffix}.txt"
   ls $filename1

   file1="$fileout_head-atm_$n1-kcld_${n2}-iciwp_${n4}-rei_q$n5.txt"
   cat $filename1 >> $file1 || exit 1
   #cat $filehead*$n1*$n2*$n5*${out_suffix}.txt >> $file1 || exit 1
   echo $file1
   echo "--------------"
done
done
done
done


exit 0

