process analysis {
    label 'r' 
    storeDir "${params.tmp_storage}/${id}/${fasta_type}/R_analysis/"
    maxForks 100
    input:
        //tuple val(id), path(prokka_faa), path(prokka_gff), path(bakta_faa), path(bakta_gff), path(eggnog_faa), path(eggnog_gff), path(pgap_faa), path(pgap_gff)
        tuple val(id), val(type1), path(prokka_1), path(prokka_2), val(type2), path(bakta_1), path(bakta_2), val(type3), path(eggnog_1), path(eggnog_2), val(type4), path(pgap_1), path(pgap_2)
        val(fasta_type)
    output:
        //Rscript R_comparison.R --args ${id} ${prokka_faa} ${prokka_gff} ${bakta_faa} ${bakta_gff} ${eggnog_faa} ${eggnog_gff} ${pgap_faa} ${pgap_gff}
        path("comparison.csv")
        publishDir "${params.output}/${id}/analysis_report_${fasta_type}/", mode: 'copy'
    script:
        """
        R_comparison.R ${id} ${prokka_1} ${prokka_2} ${bakta_1} ${bakta_2} ${eggnog_1} ${eggnog_2} ${pgap_1} ${pgap_2} ${fasta_type}
        """
    }
