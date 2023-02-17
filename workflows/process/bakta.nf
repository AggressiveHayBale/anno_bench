process bakta_database {
    label 'ubuntu'
    storeDir "${params.databases}/bakta"   

    output:
        path("db.tar.gz")
    script:
        """
        wget --no-check-certificate https://zenodo.org/record/7025248/files/db.tar.gz
        """
    stub:
        """
        touch db.tar.gz
        """
    }

process bakta {
    label 'bakta'
    storeDir "${params.tmp_storage}/bakta"
    maxForks 200
    input: 
        tuple val(name), val(species), path(fasta), val(type)
        path(bakta_db_dir)
    output: 
    	tuple val(name), val(type), val("bakta"), path("${type}_${name}_bakta.faa"),path("${type}_${name}_bakta.gff3"), emit: annotation_bakta
        publishDir "${params.output}/${name}/bakta", mode: 'copy' 
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
