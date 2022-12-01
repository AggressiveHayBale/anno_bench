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
    errorStrategy 'ignore'
    input: 
        tuple val(name), val(species), path(fasta), val(type)
        path(pgap_db)

    output: 
        tuple val(name), path("${type}_${name}_pgap.faa"), emit: gff
        tuple val(name), path("${type}_${name}_pgap.gbk"), emit: gbk
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
        mv annot.gbk ${type}_${name}_pgap.gbk

        """
    stub:
        """ 
        touch ${type}_${name}_pgap.faa \
              ${type}_${name}_pgap.gbk
        """
}
