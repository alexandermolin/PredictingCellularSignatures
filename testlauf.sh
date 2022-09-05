#!bin/bash






#working direktory und dinge zum anpassen:
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
OUTPUT_DIR=${SCRIPT_DIR}


DOCKER_CS_USERNAME=alexander.molin@stud.fh-campuswien.ac.at
DOCKER_CS_TOKEN=2f1140dbb675b9d8bb2679aaeca7e9e1



#download raw data sc:
echo "sc daten muessen immer noch manuell runtergeladen werden. und jeweils ein Tabulatorzeichen (\t) an den anfang der .txt setzen"


#download Homo_sapiens.GRCh38.cdna.all.fa.gz for Kallisto Index:
if [ ! -f $SCRIPT_DIR/Homo_sapiens.GRCh38.cdna.all.fa.gz ]
then
	wget ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
else
    echo "Homo_sapiens.GRCh38.cdna.all.fa.gz already found, if you want to download the data again, delete it."
fi


#download raw data bulk:
echo "Starte Daten-download"

for i in \
SRR12363244 SRR12363245 SRR12363246 SRR12363247 SRR12363248 \
SRR13632931 SRR13632932 SRR13632933 SRR13632934 SRR13632935 \
SRR16760002 SRR16760003 SRR16760004 SRR16760005 SRR16760006 SRR16760007 SRR16760008 SRR16760009 SRR16760010 SRR16760012
do

if [ ! -f $SCRIPT_DIR/output_kallisto_$i/abundance.tsv ]
then
	fastq-dump --gzip --split-3 $i
else
    echo "output_kallisto_$i already found, if you want to download the data again, delete the output_kallisto-folder."
fi

done

echo "Bulk-Daten healthy fertig"

for i in \
SRR13632931 SRR13632932 SRR13632933 SRR13632934 SRR13632935 \
SRR14788899 SRR14788897 SRR14788895 SRR14788893 SRR14788891 
do

if [ ! -f $SCRIPT_DIR/output_kallisto_$i/abundance.tsv ]
then
	fastq-dump --gzip --split-3 $i
else
    echo "output_kallisto_$i already found, if you want to download the data again, delete the output_kallisto-folder."
fi

done

echo "Bulk-Daten unhealthy fertig"

echo "Starte mapping"

#Indexfile for Kallisto mapping
if [ ! -f ${SCRIPT_DIR}/transcripts.idx ]
then
	kallisto index -i transcripts.idx Homo_sapiens.GRCh38.cdna.all.fa.gz
else
    echo "Kallisto Indexfile already found, if you want a new one, delete the old."
fi


#fastqc all .fastq in wd:
find ${SCRIPT_DIR} -name "*.fastq.gz" | while read READ
do
    SAMPLE_NAME_FASTQ_GZ=$(echo $(basename $READ))
    SRR_NAME=$(echo ${SAMPLE_NAME_FASTQ_GZ}| cut "-c1-11")

    if [ ! -f ${SRR_NAME}_1.fastq.gz ]
    then 
        if [ ! -f ${OUTPUT_DIR}/${SRR_NAME}_fastqc.html ]
        then
            fastqc ${SAMPLE_NAME_FASTQ_GZ} -o ./
        else
            echo "$SRR_NAME fastqc-file already found, del file for new output"
        fi
    else
        if [ ! -f ${OUTPUT_DIR}/${SRR_NAME}_2_fastqc.html ]
        then
            fastqc ${SRR_NAME}_1.fastq.gz ${SRR_NAME}_2.fastq.gz -o ./
        else
            echo "$SRR_NAME fastqc-file already found, del file for new output"
        fi
    fi


    #mapping all .fastq in wd:
    if [ ! -f output_kallisto_${SRR_NAME}/abundance.tsv ]
    then
        if [ ! -f ${SRR_NAME}_1.fastq.gz ]
        then
            kallisto quant -i transcripts.idx -o output_kallisto_${SRR_NAME} -b 100 --single -l 100 -s 5 ${SAMPLE_NAME_FASTQ_GZ}
        else
            kallisto quant -i transcripts.idx -o output_kallisto_${SRR_NAME} -b 100 ${SRR_NAME}_1.fastq.gz ${SRR_NAME}_2.fastq.gz
        fi
    else 
        echo "$SRR_NAME kallisto-file already found, del folder for new output"
    fi

done

#MultiQC quality control:
if [ ! -f $SCRIPT_DIR/multiqc_report.html ]
then
	multiqc . -f ${OUTPUT_DIR}
else
    echo "MultiQC already found, if you want a new one, delete the old."
fi


#preparing files for cibersort:
if [ ! -f $SCRIPT_DIR/cs_in/sc_for_signature_matrix.txt ]
then
	Rscript code_mapping_bulk_data_and_prep_sc_for_sig_matrix.R
    #restore order in the filesystem:
    mkdir ./cs_in
    mkdir ./cs_out
    mv gene_subset_file.txt ./cs_in/
    mv mixture* ./cs_in/
    mv sc_for_signature_matrix.txt ./cs_in
    perl -p -i -e 's/\n/\t/g;' merged_classes.txt
    sed -i -e '$a\' merged_classes.txt
    sed -i '$ s/.$//' merged_classes.txt
    mv merged_classes.txt ./cs_in
    mv sc_for_signature_matrix_merged.txt ./cs_in
else
    echo "sc_for_signature_matrix.txt already found, if you want a new one, delete the old."
fi


#cs via docker
#first run for signature matrix and fractions:

DOCKER_INPUT_DIR=${SCRIPT_DIR}/cs_in
DOCKER_OUTPUT_DIR=${SCRIPT_DIR}/cs_out
DOCKER_REFERENCE_FILE=sc_for_signature_matrix.txt
DOCKER_MIXTURE_FILE_NAME=mixture_all.txt

if [ ! -f $SCRIPT_DIR/cs_in/fraction_plots.pdf ]
then
    sudo docker run \
        -v ${DOCKER_INPUT_DIR}:/src/data \
        -v ${DOCKER_OUTPUT_DIR}:/src/outdir \
        cibersortx/fractions \
        --username ${DOCKER_CS_USERNAME} \
        --token ${DOCKER_CS_TOKEN} \
        --single_cell TRUE \
        --refsample ${DOCKER_REFERENCE_FILE} \
        --mixture ${DOCKER_MIXTURE_FILE_NAME} \
        --perm 500 \
        --rmbatchBmode TRUE

    #mv cs_out/CIBERSORTx_Results.txt ./cs_in/cs_fractions_mixture_all_table_first.txt
    mv cs_out/CIBERSORTx_Adjusted.txt ./cs_in/cs_fractions_mixture_all_table_first.txt
    Rscript code_create_fraction_plot.R
    mv cs_out/CIBERSORTx_sc_for_signature_matrix_inferred_phenoclasses.CIBERSORTx_sc_for_signature_matrix_inferred_refsample.bm.K999.txt ./cs_in/cs_signatur_matrix.txt
    mv cs_out/CIBERSORTx_sc_for_signature_matrix_inferred_phenoclasses.CIBERSORTx_sc_for_signature_matrix_inferred_refsample.bm.K999.pdf ./cs_in/cs_signature_matrix_heatmap.pdf
    mv cs_out/CIBERSORTx_sc_for_signature_matrix_inferred_phenoclasses.txt ./cs_in/cs_class_labels_file.txt
    mv cs_out/CIBERSORTx_sc_for_signature_matrix_inferred_refsample.txt ./cs_in/cs_reference_sample_file
    mv fraction_plots.pdf ./cs_in/
    mv cs_out/CIBERSORTx_cell_type_sourceGEP.txt ./cs_in/
    mv cs_out/CIBERSORTx_Mixtures-Adjusted.txt ./cs_in/cs_mixture_all_adjusted_first.txt
    mv cs_out/CIBERSORTx_Fractions-Adjusted.txt ./cs_in/cs_fractions_all_adjusted_merged_first.txt
else
    echo "fraction_plots.pdf file already found, if you want a new one, delete the old."
fi

#second run for signature matrix_merged and fractions:

DOCKER_INPUT_DIR=${SCRIPT_DIR}/cs_in
DOCKER_OUTPUT_DIR=${SCRIPT_DIR}/cs_out
DOCKER_REFERENCE_FILE=sc_for_signature_matrix_merged.txt
DOCKER_MIXTURE_FILE_NAME=mixture_all.txt

if [ ! -f $SCRIPT_DIR/cs_in/fraction_plots_merged.pdf ]
then
    sudo docker run \
        -v ${DOCKER_INPUT_DIR}:/src/data \
        -v ${DOCKER_OUTPUT_DIR}:/src/outdir \
        cibersortx/fractions \
        --username ${DOCKER_CS_USERNAME} \
        --token ${DOCKER_CS_TOKEN} \
        --single_cell TRUE \
        --refsample ${DOCKER_REFERENCE_FILE} \
        --mixture ${DOCKER_MIXTURE_FILE_NAME} \
        --perm 500 \
        --rmbatchBmode TRUE

    #mv cs_out/CIBERSORTx_Results.txt ./cs_in/cs_fractions_mixture_all_table_first.txt
    mv cs_out/CIBERSORTx_Adjusted.txt ./cs_in/cs_fractions_mixture_all_table_first_merged.txt
    Rscript code_create_fraction_plot.R
    mv cs_out/CIBERSORTx_sc_for_signature_matrix_merged_inferred_phenoclasses.CIBERSORTx_sc_for_signature_matrix_merged_inferred_refsample.bm.K999.txt ./cs_in/cs_signatur_matrix_merged.txt
    mv cs_out/CIBERSORTx_sc_for_signature_matrix_merged_inferred_phenoclasses.CIBERSORTx_sc_for_signature_matrix_merged_inferred_refsample.bm.K999.pdf ./cs_in/cs_signature_matrix_heatmap_merged.pdf
    mv cs_out/CIBERSORTx_sc_for_signature_matrix_merged_inferred_phenoclasses.txt ./cs_in/cs_class_labels_file_merged.txt
    mv cs_out/CIBERSORTx_sc_for_signature_matrix_merged_inferred_refsample.txt ./cs_in/cs_reference_sample_file_merged.txt
    mv fraction_plots.pdf ./cs_in/fraction_plots_merged.pdf
    mv cs_out/CIBERSORTx_cell_type_sourceGEP.txt ./cs_in/CIBERSORTx_cell_type_sourceGEP_merged.txt
    mv cs_out/CIBERSORTx_Mixtures-Adjusted.txt ./cs_in/cs_mixture_all_adjusted_merged.txt
    
else
    echo "fraction_plots_merged.pdf file already found, if you want a new one, delete the old."
fi

#third run for Gene express profiles:

DOCKER_SIGNATURE_MATRIX=cs_signatur_matrix_merged.txt
#DOCKER_CLASSES_FILE=merged_classes.txt

if [ ! -f $SCRIPT_DIR/something.txt ]
then
    find ${SCRIPT_DIR}/cs_in/ -name "mixture*" | while read READ
    do
    DOCKER_MIXTURE_FILE_NAME=$(echo $(basename $READ))
   
        sudo docker run \
            -v ${DOCKER_INPUT_DIR}:/src/data \
            -v ${DOCKER_OUTPUT_DIR}:/src/outdir \
            cibersortx/gep \
            --username ${DOCKER_CS_USERNAME} \
            --token ${DOCKER_CS_TOKEN} \
            --mixture ${DOCKER_MIXTURE_FILE_NAME} \
            --sigmatrix ${DOCKER_SIGNATURE_MATRIX} \
            --rmbatchBmode TRUE 

            #--classes ${DOCKER_CLASSES_FILE} \
            #--threads 8
            
        mv cs_out/*_SM_GEPs_Filtered.txt ./cs_in/cs_sm_geps_filtered_${DOCKER_MIXTURE_FILE_NAME}
        mv cs_out/*_Fractions-Adjusted.txt ./cs_in/cs_fractions_table_${DOCKER_MIXTURE_FILE_NAME}
        mv cs_out/*_GEPs.txt ./cs_in/cs_geps_${DOCKER_MIXTURE_FILE_NAME}
        mv cs_out/*_GEPs_StdErrs.txt ./cs_in/cs_geps_stderrs_${DOCKER_MIXTURE_FILE_NAME}
        mv cs_out/*_GEPs_Pvals.txt ./cs_in/cs_geps_pvalues_${DOCKER_MIXTURE_FILE_NAME}
        mv cs_out/*_GEPs_Qvals.txt ./cs_in/cs_geps_qvalues_${DOCKER_MIXTURE_FILE_NAME}
        mv cs_out/*_GEPs_CVs.txt ./cs_in/cs_geps_coefvariation_${DOCKER_MIXTURE_FILE_NAME}
        mv cs_out/*_GEPs_ThresholdPlots.pdf ./cs_in/cs_geps_threshold_${DOCKER_MIXTURE_FILE_NAME}_plot.pdf
        mv cs_out/CIBERSORTxGEP_Weights.txt ./cs_in/cs_geps_weights_${DOCKER_MIXTURE_FILE_NAME}
        mv cs_out/*_GEPs_Filtered.txt ./cs_in/cs_geps_filtered_${DOCKER_MIXTURE_FILE_NAME}
        mv cs_out/*GEP_Mixture-Adjusted.txt ./cs_in/cs_adjusted_merged_${DOCKER_MIXTURE_FILE_NAME}
    done
    Rscript code_compare_geps.R
    mv GEP_heatmap_from_r.pdf ./cs_in/
else
    echo "GEP_heatmap_from_r.pdf already found, if you want a new one, delete the old."
fi

#fourth run for hi-res:


DOCKER_INPUT_DIR=${SCRIPT_DIR}/cs_in
DOCKER_OUTPUT_DIR=${SCRIPT_DIR}/cs_out
DOCKER_REFERENCE_FILE=sc_for_signature_matrix.txt
DOCKER_MIXTURE_FILE_NAME=mixture_all_sort.txt
DOCKER_SIGNATURE_MATRIX=cs_signature_matrix.txt
DOCKER_CLASSES_FILE=merged_file.txt


if [ ! -f $SCRIPT_DIR/something.txt ]
then
    sudo docker run \
        -v ${DOCKER_INPUT_DIR}:/src/data \
        -v ${DOCKER_OUTPUT_DIR}:/src/outdir \
        cibersortx/hires \
        --username ${DOCKER_CS_USERNAME} \
        --token ${DOCKER_CS_TOKEN} \
        --mixture ${DOCKER_MIXTURE_FILE_NAME} \
        --sigmatrix ${DOCKER_SIGNATURE_MATRIX} \
        --rmbatchBmode TRUE \
        --classes ${DOCKER_CLASSES_FILE} 

         #--subsetgenes ${DOCKER_GENES_SUBSET}

else
    echo "something.txt file already found, if you want a new one, delete the old."
fi

