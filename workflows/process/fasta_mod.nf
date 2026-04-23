process fasta_mod {
    label 'ubuntu'
    //for anno_bench_tmp 
    //For the new ones
 //   storeDir "${params.tmp_storage}${name}/fasta/"
    //Dunno what this one is
   // storeDir "${params.tmp_storage}/${name}/fasta_mod/"
    // Old one for anno_bench_tmp
    storeDir "${params.tmp_storage}/fasta/${name}"
    maxForks 100
    publishDir "${params.output}/${name}/fasta", mode: 'copy'
    input: 
    // [GCA000012145, Rickettsia felis, cont_noise/anno_bench/bin/GCA_000012145.1.fna, 0.005, 0.01, 0.02]
    // path(dir) to val(dir) to ignore the need for local storing the fastas and use GCS instead
        tuple val(name), val(species), val(dir), val(noise), val(noise2), val(noise3)

    output: 
        // tuple val(name), val(species), path("original_${name}.fasta"), val("original"), emit: original
        // tuple val(name), val(species), path("noise_${noise}_${name}.fasta"), val("noise"), emit: noise
        // tuple val(name), val(species), path("noise_${noise2}_${name}.fasta"), val("noise2"), emit: noise2
        // tuple val(name), val(species), path("noise_${noise3}_${name}.fasta"), val("noise3"), emit: noise3
        tuple val(name), val(species), path("original_${name}.fasta"), emit: original
        tuple val(name), val(species), path("noise_${noise}_${name}.fasta"), emit: noise
        tuple val(name), val(species), path("noise_${noise2}_${name}.fasta"), emit: noise2
        tuple val(name), val(species), path("noise_${noise3}_${name}.fasta"), emit: noise3
    //Needed if one line fastas are included 
    //bash fasta_normalisation.sh ${dir} > original_${name}.fasta
    //bash splitter.sh ${name} ${dir} ${contig}
    script: 
    """ 
        cat ${dir} > original_${name}.fasta
        bash noise.sh "${name}" "${dir}" "${noise}"
        bash noise.sh "${name}" "${dir}" "${noise2}"
        bash noise.sh "${name}" "${dir}" "${noise3}"
    """
    stub:
    """
        touch "original_${name}.fasta"
        touch "noise_${noise}_${name}.fasta"
        touch "noise_${noise2}_${name}.fasta"
        touch "noise_${noise3}_${name}.fasta"
    """
}
