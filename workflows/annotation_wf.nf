include { prokka } from './process/prokka.nf' 
include{ bakta } from './process/bakta.nf'
include{ bakta_database } from './process/bakta.nf'
include{ eggnog } from './process/eggnog.nf'
include{ eggnog_database } from './process/eggnog.nf'
include{ pgap } from './process/pgap.nf'
include{ pgap_database } from './process/pgap.nf'
workflow annotation_wf{
    take: 
        combined_fasta_ch // val(name), val(species) val(name), val(species),path(dir), path(contigsplit), path(noise)
    main: 
    //prokka

    prokka(combined_fasta_ch)

    //prokka_outs = prokka.out.gbk.groupTuple(by: 0).join(prokka.out.faa.groupTuple(by: 0))  // val(name), [path(fasta1),path(fasta2),path(fasta3)]
    prokka_outs = prokka.out.gbk.mix(prokka.out.faa).groupTuple(by: 0)
    //bakta

    if (params.bakta_db) { bakta_db = file(params.bakta_db) }
        else { bakta_db = bakta_database() }
    bakta(combined_fasta_ch, bakta_db)

    //bakta_outs = bakta.out.gbff.groupTuple(by: 0).join(bakta.out.faa.groupTuple(by: 0))  // val(name), [path(fasta1),path(fasta2),path(fasta3)]
    bakta_outs = bakta.out.gbff.mix(bakta.out.faa).groupTuple(by: 0)

    //eggnog

    if (params.eggnog_db) { eggnog_db = file(params.egg_db) }
        else { eggnog_db = eggnog_database() } 
    eggnog(combined_fasta_ch, eggnog_db)
    eggnog_outs = eggnog.out.gff.groupTuple(by: 0).mix(eggnog.out.genepred_fasta.groupTuple(by: 0))  // val(name), [path(fasta1),path(fasta2),path(fasta3)]

    //pgap
    if (params.pgap_db) { pgap_db = file(params.pgap_db) }
        else { pgap_db = pgap_database(params.pgap_v) }

    pgap(combined_fasta_ch, pgap_db)
    pgap_outs = pgap.out.gbk.groupTuple(by: 0).mix(pgap.out.gff.groupTuple(by: 0))  // val(name), [path(fasta1),path(fasta2),path(fasta3)]

    combined_annotation_ch = prokka_outs.join(bakta_outs, by:0).join(eggnog_outs, by:0).join(pgap_outs, by:0).view()// val(name), [prokka], [bakta]

    emit: 
    combined_annotation_ch
}

