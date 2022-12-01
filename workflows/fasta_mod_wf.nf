include { fasta_mod } from './process/fasta_mod.nf' 

workflow fasta_mod_wf{
    take: 
        fasta 
    main: 
    
    fasta_mod(fasta)

    combined_fasta_ch = fasta_mod.out.original.mix(fasta_mod.out.split).mix(fasta_mod.out.noise)

    emit: 
    combined_fasta_ch

}

