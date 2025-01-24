process fasta_mod {
    label 'ubuntu'
    storeDir "${params.tmp_storage}/${name}/fasta_mod/"
    maxForks 100
    publishDir "${params.output}/${name}/fasta", mode: 'copy'
    storeDir "${params.tmp_storage}${name}/fasta/"
    input: 
    // [GCA000012145, Rickettsia felis, cont_noise/anno_bench/bin/GCA_000012145.1.fna, 0.005, 0.01, 0.02]
        tuple val(name), val(species), path(dir)//, val(noise), val(noise2), val(noise3)

    output: 
        // tuple val(name), val(species), path("original_${name}.fasta"), val("original"), emit: original
        // tuple val(name), val(species), path("noise_${noise}_${name}.fasta"), val("noise"), emit: noise
        // tuple val(name), val(species), path("noise_${noise2}_${name}.fasta"), val("noise2"), emit: noise2
        // tuple val(name), val(species), path("noise_${noise3}_${name}.fasta"), val("noise3"), emit: noise3
        tuple val(name), val(species), path("original_${name}.fasta"), emit: original
      //  tuple val(name), val(species), path("noise_${noise}_${name}.fasta"), emit: noise
      //  tuple val(name), val(species), path("noise_${noise2}_${name}.fasta"), emit: noise2
      //  tuple val(name), val(species), path("noise_${noise3}_${name}.fasta"), emit: noise3
    //Needed if one line fastas are included 
    //bash fasta_normalisation.sh ${dir} > original_${name}.fasta
    //bash splitter.sh ${name} ${dir} ${contig}
    script: 
    """ 
        cat ${dir} > original_${name}.fasta
    """
    stub:
    """
        touch "original_${name}.fasta"
    """
}
