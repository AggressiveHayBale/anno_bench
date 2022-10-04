process prokka {
    label 'prokka'
    input: 
        file(fasta_input)
    output: 
    	tuple val(name), path("${name}_prokka.*"), emit: prokka_ch
    script:
        """
        prokka --compliant --fast\
            --outdir \$PWD \
            --force  \
            --prefix "${name_prokka}"  \
            --quiet  ${dir} \
            --cpus ${task.cpus}
        """
}