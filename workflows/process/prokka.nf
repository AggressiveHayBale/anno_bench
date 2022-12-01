process prokka {
    label 'prokka'
    errorStrategy 'ignore'
    input: 
        tuple val(name), val(species), path(fasta), val(type)
    output: 
        tuple val(name), path("${type}_${name}_prokka.faa"), emit: faa
        tuple val(name), path("${type}_${name}_prokka.gbk"), emit: gbk
        publishDir "${params.output}/${name}/prokka", mode: 'copy' 
    script:
        """
        prokka --compliant --fast \
            --outdir output \
            --force  \
            --prefix ${type}_${name}_prokka  \
            --cpus ${task.cpus}  \
            --quiet ${fasta}

        mv output/${type}_${name}_prokka.faa ${type}_${name}_prokka.faa
        mv output/${type}_${name}_prokka.gbk ${type}_${name}_prokka.gbk  
        """
    stub: 
        """
        touch ${type}_${name}_prokka.faa \
            ${type}_${name}_prokka.gbk
        """
}