process bakta_database {
    label 'ubuntu'
    storeDir "${params.databases}/bakta"   

    output:
        path("db.tar.gz")
    script:
        """
        wget --no-check-certificate https://zenodo.org/record/5215743/files/db.tar.gz
        """
    stub:
        """
        touch db.tar.gz
        """
    }

process bakta {
    label 'bakta'
    errorStrategy 'ignore'
    input: 
        tuple val(name), val(species), path(fasta), val(type)
        path(bakta_db_dir)
    output: 
    	tuple val(name), path("${type}_${name}_bakta.faa"), emit: faa
        tuple val(name), path("${type}_${name}_bakta.gbff"), emit: gbff
        publishDir "${params.output}/${name}/bakta", mode: 'copy' 
    script: 
        """
        tar xzf ${bakta_db_dir}
        rm ${bakta_db_dir}
        amrfinder_update --force_update --database db/amrfinderplus-db/
        bakta --output output --prefix "${type}_${name}_bakta" --db \$PWD/db --keep-contig-headers --threads ${task.cpus} ${fasta}
        rm -rf db

        mv output/${type}_${name}_bakta.faa ${type}_${name}_bakta.faa
        mv output/${type}_${name}_bakta.gbff ${type}_${name}_bakta.gbff

        """
    stub: 
        """
        touch ${type}_${name}_bakta.faa \
            ${type}_${name}_bakta.gbff
        """
}
