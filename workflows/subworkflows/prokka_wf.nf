process prokka {
    label 'prokka'
    input: 
        tuple val(name), path(dir),  
    output: 
    	tuple val(name), path("${name}_prokka.*"), emit: prokka_ch
    script:
        """
        echo ${name}
        echo ${dir}
        prokka --compliant --fast\
            --outdir \$PWD \
            --force  \
            --prefix "${name}_prokka"  \
            --cpus ${task.cpus}  \
            --quiet ${dir}
        """
}