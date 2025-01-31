process bakta_database {
    label 'bakta'
    storeDir "${params.databases}/bakta"   

    output:
        path("bakta_db")
    script:
        """
        mkdir bakta_db
        wget --no-check-certificate https://zenodo.org/record/7669534/files/db.tar.gz
        tar xzf db.tar.gz -C bakta_db

        amrfinder_update --force_update --database /db/amrfinderplus-db/
	rm /db/amrfinderplus/latest
        cp -r /db/amrfinderplus/20* /db/amrfinderplus/latest
        chown -R mambauser:users /db/amrfinderplus 
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
   // errorStrategy 'retry'
   //     maxRetries 2
    input: 
        tuple val(name), val(species), path(fasta), val(type)
        path(bakta_db_dir)
    output: 
    //    	tuple val(name), val(type), val("bakta"), path("${type}_${name}_bakta.faa"),path("${type}_${name}_bakta.gff3"), emit: annotation_bakta
    	tuple val(name), val(type), path("${type}_${name}_bakta.faa"),path("${type}_${name}_bakta.gff3"), emit: annotation_bakta
    //amrfinder_update --force_update --database ${bakta_db_dir}/db/amrfinderplus-db/
    script: 
       """
       bakta --prefix "${type}_${name}_bakta" --db ${bakta_db_dir}/db --keep-contig-headers ${fasta}
       """
       stub: 
        """
        touch ${type}_${name}_bakta.faa \
       ${type}_${name}_bakta.gff3
        """
}
