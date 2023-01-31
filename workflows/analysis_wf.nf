include { analysis} from './process/analysis.nf' 
include { analysis as analysis_org} from './process/analysis.nf' 
include { analysis as analysis_noise} from './process/analysis.nf' 
include { analysis as analysis_split} from './process/analysis.nf' 
workflow analysis_wf{
    take: 
    combined
    main: 
    //OTHER APPROACH
    // branched_ch = combined.branch {
    //          original: it[1][0] == "original"
    //          noise: it[1][0] == "noise"
    //          split: it[1][0] == "split"
    //      }

    // original_grouped = branched_ch.original.groupTuple(by:0).view()
    // noise_grouped = branched_ch.noise.groupTuple(by:0)
    // split_grouped = branched_ch.split.groupTuple(by:0)
    // analysis(original_grouped)
    //END OF OTHER APPROACH
    //org_ch=combined.filter{ it -> it[1][0][0] == 'noise'}.view()
    //gives prokka
    //combined.map { it -> tuple (it[1]) }.view()
    //gives bakta
    //combined.map { it -> tuple (it[2]) }.view()
        
    original=combined.map { it ->
        tuple( it[0], it[1].findAll { it.first() == "original" },it[2].findAll { it.first() == "original" },
        it[3].findAll { it.first() == "original" },it[4].findAll { it.first() == "original" } )
    }

    noise=combined.map { it ->
        tuple( it[0], it[1].findAll { it.first() == "noise" },it[2].findAll { it.first() == "noise" },
        it[3].findAll { it.first() == "noise" },it[4].findAll { it.first() == "noise" } )
    }

    split=combined.map { it ->
        tuple( it[0], it[1].findAll { it.first() == "split" },it[2].findAll { it.first() == "split" },
        it[3].findAll { it.first() == "split" },it[4].findAll { it.first() == "split" } )
    }
    original_proc=original.map { it -> tuple(it[0],
    it[1].flatten()[1],it[1].flatten()[2],it[1].flatten()[3],
    it[2].flatten()[1],it[2].flatten()[2],it[2].flatten()[3],
    it[3].flatten()[1],it[3].flatten()[2],it[3].flatten()[3],
    it[4].flatten()[1],it[4].flatten()[2],it[4].flatten()[3])}
    
    noise=noise.map { it -> tuple(it[0],
    it[1].flatten()[1],it[1].flatten()[2],it[1].flatten()[3],
    it[2].flatten()[1],it[2].flatten()[2],it[2].flatten()[3],
    it[3].flatten()[1],it[3].flatten()[2],it[3].flatten()[3],
    it[4].flatten()[1],it[4].flatten()[2],it[4].flatten()[3])}
    
    split=split.map { it -> tuple(it[0],
    it[1].flatten()[1],it[1].flatten()[2],it[1].flatten()[3],
    it[2].flatten()[1],it[2].flatten()[2],it[2].flatten()[3],
    it[3].flatten()[1],it[3].flatten()[2],it[3].flatten()[3],
    it[4].flatten()[1],it[4].flatten()[2],it[4].flatten()[3])}
    
    report_original=analysis_org(original_proc,'original')
    report_noise=analysis_noise(noise,'noise')
    report_split=analysis_split(split,'split')

    report_original.mix(report_noise).mix(report_split).collectFile(name: 'collected_reports.csv', keepHeader: true ,storeDir: '${params.output}/final_report/')

    emit: 
    combined

}

