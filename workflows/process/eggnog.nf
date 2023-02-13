process eggnog_database {
    label 'eggnog'
    storeDir "${params.databases}/eggnog"
   
    output: 
        path("eggnog_db")
    script:
        """
        mkdir eggnog_db
        download_eggnog_data.py -y --data_dir eggnog_db/
        """  
    stub:
        """
        touch eggnog_db
        """

}

process eggnog {
    label 'eggnog'
    input:
        tuple val(name), val(species), path(fasta), val(type)
        path(eggnog_db_dir)
    output: 
    	tuple val(name), val(type),val("eggnog"),path("${type}_${name}_eggnog.emapper.genepred.fasta"), path("${type}_${name}_eggnog.emapper.decorated.gff"), emit: annotation_eggnog
//        tuple val(name), val(type), path("${type}_${name}_eggnog.emapper.decorated.gff"), emit: gff
        publishDir "${params.output}/${name}/eggnog", mode: 'copy' 
    script:
        """
        emapper.py -i ${fasta} \
            --cpu ${task.cpus} \
            -m diamond \
            --data_dir ${eggnog_db_dir} \
            --itype genome \
            --genepred prodigal \
            --dmnd_ignore_warnings \
            --translate \
            --go_evidence non-electronic \
            --pfam_realign none \
            --report_orthologs \
            --decorate_gff yes \
            --evalue 0.001 \
            --score 60 \
            --pident 40 \
            --query_cover 20 \
            --subject_cover 20 \
            --tax_scope auto \
            --target_orthologs all \
            -o ${type}_${name}_eggnog

        """ 
    stub: 
        """
        touch ${type}_${name}_eggnog.emapper.genepred.fasta \
            ${type}_${name}_eggnog.emapper.decorated.gff
        """
}