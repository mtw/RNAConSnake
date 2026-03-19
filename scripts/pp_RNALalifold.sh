#!/usr/bin/env bash
# Last changed Time-stamp: <2024-11-28 10:45:48 mtw>

# This can be used as a RNALalifold post-processor:
# this script post-processes a set of Stockholm files by doing
#   - determine maximum covariation level
#   - run RNAz on the alignment
#   - create covariance model from mlocarna stk (optional)


source $RNACONSERVATIONDIR/bin/RNAconservation_resources.sh

function check_infernal {
    RC_cmbuild=$(which cmbuild)
    RC_cmcalibrate=$(which cmcalibrate)

    if [ ! -f "${RC_cmbuild}" ]
    then
	echo "cmbuild (Infernal) not not found, exiting ..."
	exit 3
     fi

    if [ ! -f "${RC_cmcalibrate}" ]
    then
	echo "cmcalibrate (Infernal) not not found, exiting ..."
	exit 3
     fi
}

wd=$(pwd)
CM=OFF         # flag for building a CM
STRIP=OFF      # flag for stripping redundant lines from MSA
REMGAP=OFF     # flag for removing gap-only/high-gap fraction lines from MSA
LOCARNATE=OFF  # flag for making RNAz LocARNA-aware
DUMPFIRST=OFF   # flag for dumping the structure of the the first sequence in the input alignment
log="${wd}/pp_RNALalifold.log"
logcsv="${wd}/pp_RNALalifold.log.csv"
touch $log
touch $logcsv

for i in "$@"
do
    case $i in
      -f)             # dump first structure of alignment to output
    DUMPFIRST=ON
    shift
    ;;

    -g)             # remove gap-only/high-gap lines from MSA
	REMGAP=ON
	shift
	;;

    -l)             # make RNAz LocARNA-aware
	LOCARNATE=ON
	shift
	;;
    -s)
	STRIP=ON    # strip MSA, i.e. make non-redundant
	shift
	;;
    -cm)
	CM=ON     # special aptamer treatment if OFF by default
	check_infernal
        shift
        ;;
    *)
        # unknown option
        ;;
  esac
done

for w in $(ls *stk)
do
  wbn=$(basename $w .stk)
  echo "[I pp_RNALalifold] processing $wbn"

  cp ${wbn}.stk ${wbn}.orig.stk

  # optional: remove gaps
  if [[ "$REMGAP" == "ON" ]]; then
      ${RC_removegaponly} -a ${wbn}.stk -i stockholm -r 0.5 > ${wbn}.remgap.stk 2> /dev/null
      cp -f ${wbn}.remgap.stk ${wbn}.stk
  fi
  #  $reformataln -a $w.clean.aln -i clustal -o clustal > $w.reformat.aln

  # optional: remove duplicate sequence from alignment
  if [[ "$STRIP" == "ON" ]]; then
      ${RC_stripaln} -a ${wbn}.stk -f S --nosingle > ${wbn}.strip.stk
      cp -f ${wbn}.strip.stk ${wbn}.stk
  fi

  ${RC_eslreformat} clustal ${wbn}.stk > ${wbn}.aln

  # optional: realign stripped MSA with mlocarna here

  # run RNAz
  RNAZ_OPTIONS="-d"
  if [[ "$LOCARNATE" == "ON" ]]; then
      RNAZ_OPTIONS="${RNAZ_OPTIONS} -l"
  fi
  ${RC_rnaz} ${RNAZ_OPTIONS} ${wbn}.aln > ${wbn}.rnaz.txt
  if [[ -s "${wbn}.rnaz.txt" ]]; then # grep SVM class probability
      rnazprob=$(grep probability ${wbn}.rnaz.txt | perl -ne 's/[a-zA-Z\-\:\s]+//g;print')
  else
      rnazprob="-1"
  fi

  # run alifoldz.pl
  ${RC_alifoldz} -f -t 0.0 < ${wbn}.aln > ${wbn}.alifoldz.txt
  alifoldzscore=$(tail -n 1 ${wbn}.alifoldz.txt | grep -m 1 '.')

  # run RNAalifold
  RNAALIFOLD_OPTIONS="-t4 --aln --color -r --cfactor 0.6 --nfactor 0.5 -p --aln-EPS-cols=200 --aln-stk -f S"
  ${RC_rnaalifold} ${RNAALIFOLD_OPTIONS} < ${wbn}.stk > ${wbn}.alifold.out
  mv RNAalifold_results.stk ${wbn}.RNAalifold_results.stk
  #mv alidot.ps ${wbn}.alidot.ps
  #mv alifold.out ${wbn}.alifold.out
  #mv alirna.ps ${wbn}.alirna.ps
  #mv aln.ps ${wbn}.aln.ps
  magick ${wbn}_aln.ps ${wbn}_aln.png
  magick ${wbn}_ss.ps ${wbn}_ss.png


  # convert RNAalifold output to clustal format
  ${RC_eslreformat} clustal ${wbn}.RNAalifold_results.stk > ${wbn}.RNAalifold_results.aln
  cp ${wbn}.RNAalifold_results.aln ${wbn}.RNAalifold_results.aln~
  # remove conservation line from clustal files
  perl -i -ne 'print unless /\*+/' ${wbn}.RNAalifold_results.aln
  # remove white-space line from the end of the  clustal file
  perl -0777 -i -pe 's/\s+\z/\n/' ${wbn}.RNAalifold_results.aln
  # refold RNAalifold output
  ${RC_refold} ${wbn}.RNAalifold_results.aln ${wbn}_dp.ps | ${RC_rnafold} --noPS -C > ${wbn}_refold.out

  # extract first sequence/structure from refold output file
  refold_firstseq=$(cat ${wbn}_refold.out | head -2 | tail -1)
  refold_firststruc=$(cat ${wbn}_refold.out | head -3 | tail -1)
  alifold_consstruc=$(fgrep SS_cons ${wbn}.RNAalifold_results.stk | perl -ane 'print $F[2]')

  # optional: create CMs
  if [[ "$CM" == ON ]] && [[ $( echo "$rnazprob >= 0.9" | bc -l) -eq 1  ]]
  then
    echo "Building CM for ${wbn} (rnazprob=$rnazprob)"
    ${RC_cmbuild} ${wbn}.cm ${wbn}.RNAalifold_results.stk > ${wbn}.cmbuild.out 2> ${wbn}.cmbuild.err
    ${RC_cmcalibrate} ${wbn}.cm > ${wbn}.cmcalibrate.out 2> ${wbn}.cmcalibrate.err
    touch HAVE_CM
  fi

  if [[ "$CM" == ON ]] && [[ $( echo "$alifoldzscore <= -2.0" | bc -l) -eq 1  ]]
  then
    if [[ ! -e HAVE_CM ]]
    then
      echo "Building CM for ${wbn} (alifoldzscore=$alifoldzscore)"
      ${RC_cmbuild} ${wbn}.cm ${wbn}.RNAalifold_results.stk > ${wbn}.cmbuild.out 2> ${wbn}.cmbuild.err
      ${RC_cmcalibrate} ${wbn}.cm > ${wbn}.cmcalibrate.out 2> ${wbn}.cmcalibrate.err
    fi
  fi
  if [[ -e HAVE_CM ]]
  then
    rm HAVE_CM
  fi
  #RNAplot -a ${wbn}.mlocarna.stk --aln --covar --aln-EPS-cols=300 -t 4 --auto-id --id-prefix ${wbn}.mlocarna

  # extract info from RNAalifold output
  ${RC_alifoldmaxcovar} < ${wbn}_ali.out >> ${wbn}_alifoldmaxcovar.log
  maxcovarval=$?
  nrseq=$(head -1 ${wbn}_ali.out | perl -ane 'print "$F[0]"')
  alilen=$(head -1 ${wbn}_ali.out | perl -ane 'print "$F[5]"')
  maxcovarcount=$(head -1 ${wbn}_alifoldmaxcovar.log |  perl -ane 'print "$F[1]"')

  # dump results
  echo "X $wbn $nrseq $alilen maxcovar $maxcovarval $maxcovarcount $rnazprob $alifoldzscore ${refold_firstseq} ${refold_firststruc} ${alifold_consstruc}" >> $log
  echo "$wbn;$nrseq;$alilen;$maxcovarval;$maxcovarcount;$rnazprob;$alifoldzscore;${refold_firstseq};${refold_firststruc};${alifold_consstruc}" >> $logcsv
done
cd ..
