process fasta_mod_wf {
    label 'ubuntu'
    input: 
        tuple val(name), path(dir), val(contig), val(noise)

    output: 
        tuple val(name), path(dir), path("${name}_contigsplit.fasta"), path("${name}_noise.fasta"), emit: fasta_mod_ch
    script: 
    """
        splitter.sh ${name} ${dir} ${contig}
        noise.sh ${name} ${dir} ${noise}
    """
    stub:
    """
        touch "${name}"_contigsplit.fasta
        touch "${name}"_noise.fasta
    """

