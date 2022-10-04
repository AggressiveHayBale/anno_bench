include { prokka } from './workflows/subworkflows/prokka_wf.nf' 
//include{prokka_database} from './workflows/prokka_wf.nf'
//include{bakta} from './workflows/bakta_wf.nf'
//include{bakta_database} from './workflows/bakta_wf.nf'
//include{eggnog} from './workflows/eggnog_wf.nf'
//include{eggnog_database} from './workflows/eggnog_wf.nf'

workflow annotation_wf{
    take: 
        fasta_mod_ch // val(name), path(dir), path(contigsplit), path(noise)

    main: 
    //prokka

    //if (params.prokka_db) { prokka_db = file(params.prokka_db) }
    //    else { prokka_db = prokka_database() }
    prokka(name, dir)
    //bakta

    //if (params.bakta_db) { bakta_db = file(params.bakta_db) }
    //        else { bakta_db = bakta_database() }
    //bakta(fasta, bakta_db)

    //eggnog

    //if (params.eggnog_database) { eggnog_db = file(params.bakta_db) }
    //    else { eggnog_db = eggnog_database() } 
    //eggnog(fasta, eggnog_db)
    
    //pgap
             
    
    emit: 
    annotation = prokka_ch
}

