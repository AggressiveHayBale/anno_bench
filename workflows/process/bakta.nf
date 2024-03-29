process bakta_database {
    label 'ubuntu'
    storeDir "${params.databases}/bakta"   

    output:
        path("db.tar.gz")
    script:
        """
        wget --no-check-certificate https://zenodo.org/record/7669534/files/db.tar.gz
        """
    stub:
        """
        touch db.tar.gz
        """
    }

process bakta {
    label 'bakta'
    maxForks 200
    publishDir "${params.output}/${name}/bakta", mode: 'copy' 
    storeDir "${params.tmp_storage}/${type}/${name}/bakta/"
    errorStrategy 'retry'
        maxRetries 2
    input: 
        tuple val(name), val(species), path(fasta), val(type)
        path(bakta_db_dir)
    output: 
    //    	tuple val(name), val(type), val("bakta"), path("${type}_${name}_bakta.faa"),path("${type}_${name}_bakta.gff3"), emit: annotation_bakta

    	tuple val(name), val(type), path("${type}_${name}_bakta.faa"),path("${type}_${name}_bakta.gff3"), emit: annotation_bakta
    script: 
        """
        tar xzf ${bakta_db_dir}
        rm ${bakta_db_dir}
        amrfinder_update --force_update --database db/amrfinderplus-db/
        bakta --prefix "${type}_${name}_bakta" --db \$PWD/db --keep-contig-headers --threads ${task.cpus} ${fasta}
        rm -rf db
        
        """
    stub: 
        """
        touch ${type}_${name}_bakta.faa \
            ${type}_${name}_bakta.gff3
        """
}
