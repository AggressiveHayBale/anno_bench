process fasta_mod {
    label 'ubuntu'
    input: 
        tuple val(name), val(species), path(dir), val(contig), val(noise)

    output: 
        tuple val(name), val(species), path("original_${name}.fasta"), val("original"), emit: original
        tuple val(name), val(species), path("split_${name}.fasta"), val("split"), emit: split
        tuple val(name), val(species), path("noise_${name}.fasta"), val("noise"), emit: noise
    publishDir "${params.output}/${name}/fasta", mode: 'copy'
    script: 
    """ 
        cat ${dir} > original_${name}.fasta
        bash splitter.sh ${name} ${dir} ${contig}
        bash noise.sh ${name} ${dir} ${noise}
    """
    stub:
    """
        touch "${name}"_contigsplit.fasta
        touch "${name}"_noise.fasta
    """
}
