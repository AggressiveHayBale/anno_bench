process pgap_database {
    label 'ubuntu'
    storeDir "${params.databases}/pgap"   
    input: 
        val pgap_v
    output:
        path("input-${params.pgap_v}.tgz")

    script:
        """
        wget --no-check-certificate https://s3.amazonaws.com/pgap/input-${params.pgap_v}.tgz

        """
    stub:
        
        touch input-${params.pgap_v}.tgz
        
    }

process pgap {
    label 'pgap'
    storeDir "${params.tmp_storage}/pgap"
    maxForks 300
    input: 
        tuple val(name), val(species), path(fasta), val(type)
        path(pgap_db)

    output: 
        tuple val(name), val(type),val("pgap"), path("${type}_${name}_pgap.faa"), path("${type}_${name}_pgap.gff"), emit: annotation_pgap
        publishDir "${params.output}/${name}/pgap", mode: 'copy' 
    script:
        """
        bash yaml_creator.sh ${fasta} "${species}"
        echo '${params.pgap_v}' > VERSION
        
        mkdir /pgap/pgap/input
        tar xzf ${pgap_db} --directory "/pgap/pgap/input" --strip-components=1
        rm ${pgap_db}
        
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
