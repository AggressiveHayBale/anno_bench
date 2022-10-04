process eggnog_db_download {
    label 'eggnog'
    storeDir "${params.databases}/eggnog"
   
    output: 
        path("eggnog_db")
    script:
        """
        mkdir eggnog_db
        download_eggnog_data.py -y --data_dir eggnog_db/
        """  

}

process eggnog_emapper {
    label 'eggnog'
    publishDir "${params.output}/${name}/${params.eggnogdir}", mode: 'copy', pattern: "*"
    
    input:
        tuple val(name), path(fasta)
        path(eggnog_db_dir)
    output: 
        tuple val(name), path("${name}_eggnog*"), path("eggnog_version.txt")
    script:
        """
        emapper.py -i ${fasta} \
            --cpu 14 \
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
            -o ${name}_eggnog
        
        """  
}