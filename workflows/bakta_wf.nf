process bakta {
    label 'bakta'

    take:
        file(fasta_input)
        file(database)
    script: 
        """
        tar xzf ${database}
        rm ${database}
        amrfinder_update --force_update --database db/amrfinderplus-db/
        bakta --output \$PWD --prefix ${name}_bakta --db \$PWD/db --keep-contig-headers --threads ${task.cpus} ${fasta}
        # reduce fingerprint on local systems
        rm -rf db
        """
    output:
        path("${name}_bakta*") emit: bakta_ch
}

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