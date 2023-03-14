process prokka {
    label 'prokka'
    storeDir "${params.tmp_storage}/${type}/${name}/prokka/"
    maxForks 100
    publishDir "${params.output}/${name}/prokka", mode: 'copy'
    input: 
        tuple val(name), val(species), path(fasta), val(type)
    output: 
        tuple val(name), val(type), path("${type}_${name}_prokka.faa"),path("${type}_${name}_prokka.gff"), emit: annotation_prokka
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