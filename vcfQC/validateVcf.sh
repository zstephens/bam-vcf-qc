#########################################################################
###
###		validateVcf.sh
###
###		Written by: Zach Stephens
###		Date:		09/04/2014
###		For:		Secondary Analysis Quality Control Tools
###					Mayo Clinic, Summer 2014
###
#########################################################################
###
###		./validateVcf.sh configFile.txt input.vcf outReportDir/
###
###		This script exits 0 if VCF is ok, 1 if failed.
###
###		If no third argument is provided, no report will be created.
###

CONFIG=$1
IN_VCF=$2

if [ $# -eq 3 ]
then
	REPORT_OUT=$3
else
	REPORT_OUT="None"
fi

VCF_PERL=$( cat $CONFIG | grep -w '^VCF_PERL' | cut -d '=' -f2 )
VCF_VALIDATE=$( cat $CONFIG | grep -w '^VCF_VALIDATE' | cut -d '=' -f2 )
VCF_STATS=$( cat $CONFIG | grep -w '^VCF_STATS' | cut -d '=' -f2 )

VCF_REPORT=$( cat $CONFIG | grep -w '^VCF_REPORT' | cut -d '=' -f2 )

PYTHON_EXEC=$( cat $CONFIG | grep -w '^PYTHON_EXEC' | cut -d '=' -f2 )

TEMP_DIR=$( cat $CONFIG | grep -w '^TEMP_DIR' | cut -d '=' -f2 )

TEMP_V="${TEMP_DIR}vcfTools_validate.log"
TEMP_S="${TEMP_DIR}vcfTools_stats.dat"
TEMP_E="${TEMP_DIR}vcfTools_stats.errors"

## file exists, non-zero size
if [ ! -s $IN_VCF ]
then
	echo -e "$IN_VCF : doesn't exist"
	exit 1;
fi


## vcftools vcf-validate
export PERL5LIB=$VCF_PERL
$VCF_VALIDATE $IN_VCF 2> $TEMP_V

if [ $? -eq 0 ]
then
	## vcftools vcf-stats
	$VCF_STATS $IN_VCF > $TEMP_S 2> $TEMP_E
else
	echo -n "" > $TEMP_S
fi


## generate report
$PYTHON_EXEC $VCF_REPORT $IN_VCF $TEMP_V $TEMP_S $REPORT_OUT

## if everything went as planned, copy useful logs into report directory and delete the rest
if [ $? -eq 0 ] && [ $REPORT_OUT != "None" ]
then
	mv -f $TEMP_V $REPORT_OUT/logs/vcfTools_validate.log
	mv -f $TEMP_E $REPORT_OUT/logs/vcftools_stats.errors
	rm $TEMP_S
fi
## if we had a fatal error, exit 1
if [ $? -eq 1 ]
then
	exit 1;
fi

