process pgap_database {
    label 'ubuntu'
    storeDir "${params.databases}/pgap"   
    input: 
        val pgap_v
    output:
        path("pgap_db")

    script:
        """
        wget --no-check-certificate https://s3.amazonaws.com/pgap/input-${params.pgap_v}.tgz
        tar xzf input-${params.pgap_v}.tgz -C pgap_db
        mv pgap_db/input*/* pgap_db/.
        """
    stub:
        
        touch input-${params.pgap_v}
        
    }

process pgap {
    label 'pgap'
    storeDir "${params.tmp_storage}/${type}/${name}/pgap/"
    maxForks 300
    input: 
        tuple val(name), val(species), path(fasta), val(type)
    path(pgap_db)

    output: 
        tuple val(name), val(type), path("${type}_${name}_pgap.faa"), path("${type}_${name}_pgap.gff"), emit: annotation_pgap
        publishDir "${params.output}/${name}/pgap", mode: 'copy' 
    script:
        """
        bash yaml_creator.sh ${fasta} "${species}"
        echo '${params.pgap_v}' > VERSION
        
        mv ${pgap_db} /pgap/pgap/input
        cwltool /pgap/pgap/pgap.cwl --fasta ${fasta} --ignore_all_errors --report_usage --submol meta.yaml
        mv annot.faa ${type}_${name}_pgap.faa
        mv annot.gff ${type}_${name}_pgap.gff

        """
    stub:
        """ 
        touch ${type}_${name}_pgap.faa \
              ${type}_${name}_pgap.gff
        """
}
