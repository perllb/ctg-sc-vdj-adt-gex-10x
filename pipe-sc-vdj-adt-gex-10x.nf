#!/usr/bin/env nextFlow

// Base params
runfolder = params.runfolder
basedir = params.basedir
metaid = params.metaid

// Output dirs
outdir = params.outdir
fqdir = params.fqdir
qcdir = params.qcdir
countdir = params.countdir
aggdir = params.aggdir
ctgqc = params.ctgqc
metadir = params.metadir

// Input feature reference
featureref = params.features

// Demux args
b2farg = params.bcl2fastqarg
demux = params.demux

// Delivery email (for ctg-delivery.info.csv)
email = params.email
// Read and process CTG samplesheet 
sheet = file(params.sheet)

println "============================="
println ">>> sc-cite-multi-10x pipeline >>>"
println ""
println "> INPUT: "
println ""
println "> run-meta-id		: $metaid "
println "> basedir		: $basedir "
println "> runfolder		: $runfolder "
println "> sample-sheet		: $sheet "
println "> feature-ref          : $featureref "
println ""
println " - demux settings " 
println "> bcl2fastq-arg	: '${b2farg}' "
println "> demux                : $demux " 
println "> metadir              : $metadir "
println ""
println " - output directory structure "
println "> outdir               : $outdir "
println "> fastq                : $fqdir "
println "> qc                   : $qcdir "
println "> count                : $countdir " 
println "> aggregated           : $aggdir "
println "> ctg-qc               : $ctgqc "
println "> metadata             : $metadir "
println ""
println "============================="

// set new samplesheets for channels and demux
chsheet = file("$basedir/sample-sheet.nf.channel.csv")
demuxsheet = file("$basedir/sample-sheet.nf.demux.csv")

// Read and process sample sheet                                                                           
all_lines = sheet.readLines()
write_b = false // if next lines has sample info
chsheet.text=""

for ( line in all_lines ) { 
    if ( write_b ) {
        chsheet.append(line + "\n")
	}
    if (line.contains("[Data]")){
        write_b = true
    }
}   

// extract samplesheet info
Channel
    .fromPath(chsheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.Sample_Project, row.Sample_Species, row.Sample_Lib, row.Sample_Pair ) }
    .tap{infoall}
    .into { crlib_ch; cragg_ch; fqc_ch }

// Projects
Channel
    .fromPath(chsheet)
    .splitCsv(header:true)
    .map { row -> row.Sample_Project }
    .unique()
    .tap{infoProject}
    .into { count_summarize ; deliveryInfo }


println " > Samples to process: "
println "[Sample_ID,Sample_Project,Sample_Species,Sample_Lib,pair]"
infoall.subscribe { println "Info: $it" }

println " > Projects to process : "
println "[Sample_Project]"
infoProject.subscribe { println "Info Projects: $it" }

// Create delivery info file
process delivery_info {

	tag "${metaid}-${projid}"

	input:
	val projid from deliveryInfo

	"""
	mkdir -p ${outdir}	
	deliveryinfo="${outdir}/ctg-delivery.info.csv"
	echo "projid,${projid}" > \$deliveryinfo
	echo "email,${email}" >> \$deliveryinfo
	echo "pipeline,sc-vdj-adt-gex-10x" >> \$deliveryinfo
	"""
}

// Parse samplesheet for demux
process parsesheet {

	tag "$metaid"

	input:
	val chsheet

	output:
	val demuxsheet into demux_sheet

	// Create metadir output 
	File directory = new File(metadir);
    	if (! directory.exists()){
           directory.mkdir();
	}

	"""
	python $basedir/bin/ctg-parse-samplesheet.10x.py -s $chsheet -o $demuxsheet -i 'dual'
	"""
}

// generate count libraries csv
process gen_libraries_csv {

    	tag "${sid}_${projid}"

	input:
	val sheet
	set sid, projid, ref, lib, pair from crlib_ch

 	output:
	set sid, projid, ref, lib, pair into count_lib_csv

	when:
	lib == 'gex'

	script:
	// Get species to set references
	if ( ref == "Human" || ref == "human" ){
	   genome=params.genome_h
	   vdjref=params.vdj_ref_h
}
	else if ( ref == "Mouse" || ref == "mouse"){
	   genome=params.genome_m
	   vdjref=params.vdj_ref_m
}
	else if ( ref == "custom" || ref == "Custom") {
	   genome=params.genome_custom
	   vdjref=params.vdj_ref_custom
}
	else {
	     print "> Species not recognized - check samplesheet" 
	     genome="ERR"
}

	// Create metadir output 	 
	File directory = new File(metadir);
    	if (! directory.exists()){
           directory.mkdir();
	}

	libcsv=metadir + "/" + projid + "_" + sid + "_libraries.csv"

	"""

# Gene expression
echo '[gene-expression]' > $libcsv
echo "reference,$genome" >> $libcsv

if grep "${sid}_FB" $chsheet; then 
   # Feature 
   echo '[feature]' >> $libcsv
   echo "reference,$featureref" >> $libcsv
fi

# VDJ
echo '[vdj]' >> $libcsv
echo "reference,$vdjref" >> $libcsv


# libraries
echo '[libraries]' >> $libcsv
echo 'fastq_id,fastqs,lanes,feature_types' >> $libcsv
# Add gex
echo '$sid,${fqdir}/${projid},any,Gene Expression' >> $libcsv

if grep "${sid}_TCR" $chsheet; then 
   # Add VDJ-T
   echo '${sid}_TCR,${fqdir}/${projid},any,VDJ-T' >> $libcsv
fi

if grep "${sid}_BCR" $chsheet; then 
   # Add VDJ-B
   echo '${sid}_BCR,${fqdir}/${projid},any,VDJ-B' >> $libcsv
fi

if grep "${sid}_FB" $chsheet; then 
   # Add Feature Barcode
   echo '${sid}_FB,${fqdir}/${projid},any,Antibody Capture' >> $libcsv
fi
        """
}

// Run mkFastq
process mkfastq {

	tag "${metaid}"

	input:
	val sheet from demux_sheet
	
	output:
	val "count" into count
	val "go" into fastqc_go

	when:
	demux == 'y'

	"""

cellranger mkfastq \\
	   --id=${metaid} \\
	   --run=$runfolder \\
	   --samplesheet=$sheet \\
	   --jobmode=local \\
	   --localmem=115 \\
	   --output-dir ${fqdir} \\
	   $b2farg

#Remove Undetermined fastq
rm -r -f ${fqdir}/Undetermined*
"""
}

// count multi
process count {

	tag "${sid}-${projid}"
	publishDir "${countdir}/", mode: "copy", overwrite: true

	input: 
	val ready from count
	set sid, projid, ref, lib, pair from count_lib_csv

	output:
	file "${sid}/outs/" into samplename
	val "${qcdir}/cellranger/${sid}.metrics_summary.csv" into count_metrics

	when:
	lib == 'gex'

	// Create metadir output 	 
	File directory = new File(countdir);
    	if (! directory.exists()){
           directory.mkdir();
	}

	libcsv=metadir + "/" + projid + "_" + sid + "_libraries.csv"
	"""

cellranger multi \\
	--id=$sid \\
	--csv=$libcsv \\
	--localmem=115 \\
	--jobmode=local \\
	--localcores=${task.cpus} 

## Copy h5 file for aggregation
# NA 

## Copy metrics file for qc
# Remove if it exists
if [ -f ${qcdir}/cellranger/${sid}.metrics_summary.csv ]; then
	rm -r ${qcdir}/cellranger/${sid}.metrics_summary.csv
fi
mkdir -p ${qcdir}
mkdir -p ${qcdir}/cellranger/
cp ${sid}/outs/per_sample_outs/${sid}/metrics_summary.csv ${qcdir}/cellranger/${sid}.metrics_summary.csv

## Copy to delivery folder 
mkdir -p ${outdir}/summaries
mkdir -p ${outdir}/summaries/cloupe
mkdir -p ${outdir}/summaries/web-summaries
cp ${sid}/outs/per_sample_outs/${sid}/web_summary.html ${outdir}/summaries/web-summaries/${sid}.web_summary.html

# GEX cloupe
cp ${sid}/outs/per_sample_outs/${sid}/count/cloupe.cloupe ${outdir}/summaries/cloupe/${sid}_cloupe.cloupe

# VDJ-T vloupe
if [ -e ${sid}/outs/per_sample_outs/${sid}/vdj_t/vloupe.vloupe ]; then
   cp ${sid}/outs/per_sample_outs/${sid}/vdj_t/vloupe.vloupe ${outdir}/summaries/cloupe/${sid}_vdj_t_vloupe.vloupe
fi

# VDJ-B vloupe
if [ -e ${sid}/outs/per_sample_outs/${sid}/vdj_b/vloupe.vloupe ]; then
   cp ${sid}/outs/per_sample_outs/${sid}/vdj_b/vloupe.vloupe ${outdir}/summaries/cloupe/${sid}_vdj_b_vloupe.vloupe
fi

## Copy to CTG-QC dir 
mkdir -p ${ctgqc}
mkdir -p ${ctgqc}/web-summaries
cp ${sid}/outs/per_sample_outs/${sid}/web_summary.html ${ctgqc}/web-summaries/${sid}.web_summary.html

	"""

}

process summarize_count {

	tag "${metaid}"

	input:
	val metrics from count_metrics.collect()
	val projid from count_summarize

	output:
	val "y" into mqc_count 	
	val "x" into run_summarize

	"""
cd $outdir

mkdir -p ${qcdir}
mkdir -p ${qcdir}/cellranger

# summaries
python $basedir/bin/ctg-sc-vdj-adt-gex-count-metrics-concat.py -i ${outdir} -o ${qcdir}/cellranger -p ${projid}

	"""
}

// aggregation
process gen_aggCSV {

	tag "${sid}_${projid}"

	input:
	set sid, projid, ref, lib, pair from cragg_ch

	output:
	val projid into craggregate

	when:
	lib == 'aa'

	"""
mkdir -p ${aggdir}

aggcsv=${aggdir}/${projid}_libraries.csv

if [ -f \${aggcsv} ]
then
	if grep -q $sid \$aggcsv
	then
		echo ""
	else
		echo "${sid},${aggdir}/${sid}.molecule_info.h5" >> \$aggcsv
	fi
else
	echo "sample_id,molecule_h5" > \$aggcsv
	echo "${sid},${aggdir}/${sid}.molecule_info.h5" >> \$aggcsv
fi


	"""
}

process aggregate {

	publishDir "${outdir}/aggregate/", mode: 'move', overwrite: true
	tag "$projid"
  
	input:
	val projid from craggregate.unique()
	//val moleculeinfo from count_agg.collect()

	output:
	file "${projid}_agg/outs" into doneagg
	val projid into md5_agg_go

	when:
	agg == 'aa' 

	"""
cellranger aggr \
   --id=${projid}_agg \
   --csv=${aggdir}/${projid}_libraries.csv \
   --normalize=mapped

## Copy to delivery folder 
cp ${projid}_agg/outs/web_summary.html ${outdir}/summaries/web-summaries/${projid}_agg.web_summary.html
cp ${projid}_agg/outs/count/cloupe.cloupe ${outdir}/summaries/cloupe/${projid}_agg_cloupe.cloupe

## Copy to CTG QC dir 
cp ${outdir}/summaries/web-summaries/${projid}_agg.web_summary.html ${ctgqc}/web-summaries/

## Remove the molecule_info.h5 files that are stored in the aggregate folder (the original files are still in count-cr/../outs 
rm ${aggdir}/*h5

	"""
}

process fastqc {

	tag "${sid}-${projid}"

	input:
	val go from fastqc_go
	set sid, projid, ref, lib, pair from fqc_ch	
		
	output:
	val projid into mqc_cha
	val projid into md5_fastqc_go

	"""
mkdir -p ${qcdir}
mkdir -p ${qcdir}/fastqc

for file in ${fqdir}/${projid}/${sid}*fastq.gz
	do fastqc -t ${task.cpus} \$file --outdir=${qcdir}/fastqc
done
	"""
}

process multiqc_count_run {

	tag "${metaid}"

	input:
	val x from run_summarize.collect()
	val projid from mqc_cha.collect()
		
	output:
	val "x" into summarized

	"""
cd ${outdir}
mkdir -p ${qcdir}/multiqc
multiqc -f ${fqdir} ${qcdir}/fastqc/ ${qcdir}/cellranger/ --outdir ${qcdir}/multiqc -n ${metaid}_sc-vdj-adt-gex-10x_summary_multiqc_report.html

cp -r ${qcdir} ${ctgqc}

	"""

}

// Final process, when all is done: md5 recursively from output root folder
process md5sum {

	input:
	val v from summarized.collect()
	val projid from md5_fastqc_go.unique()
	val x from md5_agg_go.collect()
	
	"""
# Remove Undetermined files!
rm ${fqdir}/adt/Undetermined*
rm ${fqdir}/rna/Undetermined*

cd ${outdir}
find -type f -exec md5sum '{}' \\; > ctg-md5.${projid}.txt
        """ 

}
