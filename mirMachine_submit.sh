#!/bin/sh


usage()
{
cat << EOF
usage: bash mirMachine_submit.sh -f fasta_file -i input_file [OPTIONS] 
Basic parameters
-f    | --fasta_file         (FASTA file.                    The reference FASTA file to search for miRNAs
                             required unless -db provided)
-db   | --blast_database     (Blast database title.          Blast database title
                             ignored if -g provided)
-i    | --input_file         (mature miRNA FASTA file        Fasta file for known mature miRNA sequences
                             Required)
-m    | --mismatches         (default to 1)                  Number of mismatches for mir_find hits
-n    | --number_of_hits     (default to 20. Set 0 for 
                             infinite number of hits)        Number of hits to eliminate as siRNA in mir_fold
-long                        (Optional)                      To assess secondary structure of the suspect list

For sRNA-seq based predictions:
-s    | --sRNAseq            (Optional)                      Indicate
-lmax                        (default to 24)
-lmin                        (default to 20)
-rpm                         (default to 10)


-h    | --help                                               Brings up this menu
EOF
}


mismatches=1
number_of_hits=20
run="short"
sRNAseq="false"
lmax=24
lmin=20
rpm=10

BASEDIR="${BASH_SOURCE[0]%/mirMachine_submit.sh}"

if [ $# -eq 0 ]; then
    usage
    exit
fi

while [ "$1" != "" ]; do
    case $1 in
        -f | --fasta_file )
            shift        
            fasta_file=$1
            blast_database=$1"_db" 
        ;;
        -db | --blast_database )
            shift        
            blast_database=$1
        ;;
        -i | --input_file )
            shift        
            input_file=$1
        ;;
        -n | --number_of_hits )
            shift        
            if [ $1 -eq 0]; then 
            	$1=""
            fi
            number_of_hits=$1
        ;;
        -m | --mismatches )
            shift
            mismatches=$1
		;;
        -long )
            run="long"
		;;
        -s | --sRNAseq )       
            sRNAseq="true"
		;;
        -lmax )
            shift        
            lmax=$1
        ;;
        -lmin )
            shift        
            lmin=$1
        ;;
        -rpm )
            shift        
            rpm=$1
        ;;        
		-h | --help )    usage
            exit
        ;;
        * )              usage
            exit 1
    esac
    shift
done


if [ -z $blast_database ]; then
    echo "Provide either a reference FASTA file (-f, --fasta_file) or a blast database (-db, --blast_database)"
    usage
    exit
fi

if [ -z $input_file ]; then
    echo "Mature miRNA sequences is required, provide it with the flag: -i input_file"
    usage
    exit
fi



echo "Mature miRNA sequences = ${input_file}"
echo "Number of mismatches   = ${mismatches}"
echo "Number of hits allowed = ${number_of_hits}"


if [ -z $fasta_file ]; then
    echo "Blast DB               = ${blast_database}"
    echo "==============================================================================="
    echo
else
    echo "Reference FASTA file   = ${fasta_file}"
    echo "Blast DB               = ${blast_database}"
    echo "==============================================================================="
    echo
    echo "Running makeblastdb"
    echo "==============================================================================="
    makeblastdb -dbtype nucl -parse_seqids -in $fasta_file -out $blast_database -title $blast_database
fi

if [ $sRNAseq == "true" ]; then
    mismatches=0;
    echo "==============================================================================="
    echo "sRNAseq preprocessing..."
    echo "==============================================================================="
    $BASEDIR/mir_filter_reads.pl $input_file $lmin $lmax $rpm

	echo
	echo "Running mir_find..."
	echo "==============================================================================="
	echo $mismatches | $BASEDIR/mir_find.pl $input_file.filtered.fasta $blast_database

	echo
	echo "Running mir_fold..."
	echo "==============================================================================="
	echo $number_of_hits | $BASEDIR/mir_fold.pl $input_file.nr.fasta $input_file.filtered.fasta.results.tbl $blast_database $run $sRNAseq
else
	echo
	echo "Running mir_find..."
	echo "==============================================================================="
	echo $mismatches | $BASEDIR/mir_find.pl $input_file $blast_database


	echo
	echo "Running mir_fold..."
	echo "==============================================================================="
	echo $number_of_hits | $BASEDIR/mir_fold.pl $input_file $input_file.results.tbl $blast_database $run $sRNAseq
fi

echo
echo "All finished!"
