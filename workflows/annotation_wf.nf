include { prokka } from './process/prokka.nf' 
include{ bakta } from './process/bakta.nf'
include{ bakta_database } from './process/bakta.nf'
include{ eggnog } from './process/eggnog.nf'
include{ eggnog_database } from './process/eggnog.nf'
include{ pgap } from './process/pgap.nf'
include{ pgap_database } from './process/pgap.nf'
workflow annotation_wf{
    take: 
        combined_fasta_ch 
    main: 
    //prokka
    prokka(combined_fasta_ch)

    prokka_outs = prokka.out.annotation_prokka.map { it -> [it[0],[it[1], it[2], it[3], it[4]]]}.groupTuple(by:0)

    //bakta

    if (params.bakta_db) { bakta_db = file(params.bakta_db) }
        else { bakta_db = bakta_database() }
    bakta(combined_fasta_ch, bakta_db)


    bakta_outs = bakta.out.annotation_bakta.map { it -> [it[0],[it[1], it[2], it[3], it[4]]]}.groupTuple(by:0)
    //eggnog

    if (params.eggnog_db) { eggnog_db = file(params.egg_db) }
        else { eggnog_db = eggnog_database() } 
    eggnog(combined_fasta_ch, eggnog_db)

    eggnog_outs = eggnog.out.annotation_eggnog.map { it -> [it[0],[it[1], it[2], it[3], it[4]]]}.groupTuple(by:0)
   
    //pgap
    if (params.pgap_db) { pgap_db = file(params.pgap_db) }
        else { pgap_db = pgap_database(params.pgap_v) }

    pgap(combined_fasta_ch, pgap_db)

    pgap_outs = pgap.out.annotation_pgap.map { it -> [it[0],[it[1], it[2], it[3], it[4]]]}.groupTuple(by:0)
 
    combined=prokka_outs.join(bakta_outs).join(eggnog_outs).join(pgap_outs)
    emit: 
    combined
}
