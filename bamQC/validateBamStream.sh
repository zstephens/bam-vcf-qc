#########################################################################
###
###		validateBamStream.sh
###
###		Written by: Zach Stephens
###		Date:		09/04/2014
###		For:		Secondary Analysis Quality Control Tools
###					Mayo Clinic, Summer 2014
###
#########################################################################
###
###		cat input.bam | ./validateBamStream.sh opMode configFile.txt outReportDir/
###
###		This script exits 0 if BAM is ok, 1 if failed.
###
###		Op-modes:
###			"lenient"	-	only check EOF, header and first record
###			"moderate"	-	check EOF + every record for basic errors
###			"strict"	-	moderate + Picard ValidateSam
###
###		If no third argument is provided, no report will be created.
###
###		example command:
###
###		cat myBam.bam | ./validateBamStream.sh moderate validateBam_config.txt outReport/	
###
###		input bam can presumably be streamed from anywhere (e.g. aligner/samtools output)
###
###

OP_MODE=$1
CONFIG=$2

if [ $# -eq 3 ]
then
	REPORT_OUT=$3
else
	REPORT_OUT="None"
fi

SAMTOOLS_EXEC=$( cat $CONFIG | grep -w '^SAMTOOLS_EXEC' | cut -d '=' -f2 )
PYTHON_EXEC=$( cat $CONFIG | grep -w '^PYTHON_EXEC' | cut -d '=' -f2 )
JAVA_EXEC=$( cat $CONFIG | grep -w '^JAVA_EXEC' | cut -d '=' -f2 )
PICARD_VALIDATESAM=$( cat $CONFIG | grep -w '^PICARD_VALIDATESAM' | cut -d '=' -f2 )

EOF_CHECK=$( cat $CONFIG | grep -w '^EOF_CHECK' | cut -d '=' -f2 )
EOF_STREAM=$( cat $CONFIG | grep -w '^EOF_STREAM' | cut -d '=' -f2 )
BAM_QC=$( cat $CONFIG | grep -w '^BAM_QC' | cut -d '=' -f2 )
PIC_FAIL=$( cat $CONFIG | grep -w '^PIC_FAIL' | cut -d '=' -f2 )
AGG_ERRORS=$( cat $CONFIG | grep -w '^AGG_ERRORS' | cut -d '=' -f2 )
GEN_REPORT=$( cat $CONFIG | grep -w '^GEN_REPORT' | cut -d '=' -f2 )

TEMP_DIR=$( cat $CONFIG | grep -w '^TEMP_DIR' | cut -d '=' -f2 )

TEMP_H="${TEMP_DIR}temp.header"
TEMP_L1="${TEMP_DIR}temp.log1"
TEMP_L2="${TEMP_DIR}temp.log2"
TEMP_L3="${TEMP_DIR}temp.log3"
TEMP_BE="${TEMP_DIR}temp.bamQCErrors"
TEMP_PE="${TEMP_DIR}temp.picardErrors"
TEMP_PL="${TEMP_DIR}temp.picardLog"
TEMP_F="${TEMP_DIR}temp.data"
TEMP_A="${TEMP_DIR}temp.allErrors"

function basicChecks {
	## check for EOF marker
	$PYTHON_EXEC $EOF_STREAM 1> $TEMP_L1
	if [ ! $? -eq 0 ]
	then
		echo -e "\nBasic checks failed: EOF marker not found\n" >> $TEMP_L1
	fi
}

function basicChecks2 {
	## check validity of header and first record
	zcat | $BAM_QC -H 1> $TEMP_L2 2> $TEMP_BE
	if [ ! $? -eq 0 ]
	then
		echo -e "\nBasic checks failed: header and/or first record corrupted\n" >> $TEMP_L2
	fi
}

function bamQCStream {
	zcat | $BAM_QC -o $TEMP_F 1> $TEMP_L2 2> $TEMP_BE
	if [ ! $? -eq 0 ]
	then
		echo -e "\nBamQC failed\n" >> $TEMP_L2
	fi
}

function picardValidateStream {
	$JAVA_EXEC -jar $PICARD_VALIDATESAM I=/dev/stdin MAX_OUTPUT=null O=$TEMP_PE MODE=SUMMARY 1> $TEMP_L3 2> $TEMP_PL
	if [ ! $? -eq 0 ]
	then
		$PYTHON_EXEC $PIC_FAIL $TEMP_PL $TEMP_PE
	fi
}

if [ $OP_MODE == "lenient" ]
then
	cat | tee >(basicChecks2) | basicChecks
	echo -n "" > $TEMP_PE
	cat $TEMP_L1
	cat $TEMP_L2

elif [ $OP_MODE == "moderate" ]
then
	cat | tee >(bamQCStream) | basicChecks
	echo -n "" > $TEMP_PE
	cat $TEMP_L1
	cat $TEMP_L2

elif [ $OP_MODE == "strict" ]
then
	cat | tee >(picardValidateStream) >(bamQCStream) | basicChecks
	cat $TEMP_L1
	cat $TEMP_L2
	cat $TEMP_L3
else
	echo "Unknown Op-Mode."
	exit 1;
fi

## fail out immediately on really bad errors
grep -q "Basic checks failed" $TEMP_L2
if [ $? -eq 0 ]
then
	echo "Basic checks failed."
	exit 1;
fi
grep -q "BamQC failed" $TEMP_L2
if [ $? -eq 0 ]
then
	echo "BamQC failed."
	exit 1;
fi

if [ $REPORT_OUT != "None" ] && [ $OP_MODE != "lenient" ]
then
	## parse error logs from BamQC and PicardValidateBam, and perform some additional checks
	$PYTHON_EXEC $AGG_ERRORS $TEMP_F $TEMP_BE $TEMP_PE $TEMP_A

	DID_WE_FAIL=$?

	## generate output report
	$PYTHON_EXEC $GEN_REPORT $TEMP_F $TEMP_A $REPORT_OUT

	## if everything went as planned, copy useful logs into report directory and delete the rest
	if [ $? -eq 0 ]
	then
		mv -f $TEMP_BE $REPORT_OUT/logs/bamQC_output.log
		mv -f $TEMP_F $REPORT_OUT/logs/raw_bamQC_metrics.dat
		rm $TEMP_A
		rm $TEMP_PE
		if [ $OP_MODE == "strict" ]
		then
			mv -f $TEMP_PL $REPORT_OUT/logs/picardValidateSam_output.log
		fi
	fi
	## if we had a fatal error, exit 1
	if [ $? -eq 1 ]
	then
		exit 1;
	fi

	## if we encountered any errors worthy of failure, exit 1
	if [ $DID_WE_FAIL -eq 1 ]
	then
		exit 1;
	else
		exit 0;
	fi
fi

