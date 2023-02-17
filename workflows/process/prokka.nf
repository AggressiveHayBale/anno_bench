process prokka {
    label 'prokka'
    storeDir "${params.tmp_storage}/prokka"
    maxForks 100
    input: 
        tuple val(name), val(species), path(fasta), val(type)
    output: 
        tuple val(name), val(type), val("prokka"), path("${type}_${name}_prokka.faa"),path("${type}_${name}_prokka.gff"), emit: annotation_prokka
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
        mv output/${type}_${name}_prokka.gff ${type}_${name}_prokka.gff  

        """
    stub: 
        """
        touch ${type}_${name}_prokka.faa \
            ${type}_${name}_prokka.gff
        """
}