include { fasta_mod } from './process/fasta_mod.nf' 

workflow fasta_mod_wf{
    take: 
        fasta 
    main: 
    
    fasta_mod(fasta)

    combined_fasta_ch = fasta_mod.out.original.mix(fasta_mod.out.noise).mix(fasta_mod.out.noise2).mix(fasta_mod.out.noise3)

    emit: 
    combined_fasta_ch

}

