include { fasta_mod } from './process/fasta_mod.nf' 

workflow fasta_mod_wf{
    take: 
        fasta 
    main: 
    
    fasta_mod(fasta)

    original = fasta_mod.out.original.map {name, species, path -> [name, species, path,"original"] }
    noise = fasta_mod.out.noise.map {name, species, path -> [name, species, path,"noise"] }
    noise2 = fasta_mod.out.noise2.map {name, species, path -> [name, species, path,"noise2"] }
    noise3 = fasta_mod.out.noise3.map {name, species, path -> [name, species, path,"noise3"] }

    //combined_fasta_ch = fasta_mod.out.original.mix(fasta_mod.out.noise).mix(fasta_mod.out.noise2).mix(fasta_mod.out.noise3)
    combined_fasta_ch = original.mix(noise).mix(noise2).mix(noise3)

    emit: 
    combined_fasta_ch

}

