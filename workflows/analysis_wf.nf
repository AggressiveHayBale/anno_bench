include { analysis} from './process/analysis.nf' 
include { analysis as analysis_org} from './process/analysis.nf' 
include { analysis as analysis_noise} from './process/analysis.nf' 
include { analysis as analysis_noise2} from './process/analysis.nf' 
include { analysis as analysis_noise3} from './process/analysis.nf' 
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

    noise2=combined.map { it ->
        tuple( it[0], it[1].findAll { it.first() == "noise2" },it[2].findAll { it.first() == "noise2" },
        it[3].findAll { it.first() == "noise2" },it[4].findAll { it.first() == "noise2" } )
    }

    noise3=combined.map { it ->
        tuple( it[0], it[1].findAll { it.first() == "noise3" },it[2].findAll { it.first() == "noise3" },
        it[3].findAll { it.first() == "noise3" },it[4].findAll { it.first() == "noise3" } )
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
    
    noise2=noise2.map { it -> tuple(it[0],
    it[1].flatten()[1],it[1].flatten()[2],it[1].flatten()[3],
    it[2].flatten()[1],it[2].flatten()[2],it[2].flatten()[3],
    it[3].flatten()[1],it[3].flatten()[2],it[3].flatten()[3],
    it[4].flatten()[1],it[4].flatten()[2],it[4].flatten()[3])}

    noise3=noise3.map { it -> tuple(it[0],
    it[1].flatten()[1],it[1].flatten()[2],it[1].flatten()[3],
    it[2].flatten()[1],it[2].flatten()[2],it[2].flatten()[3],
    it[3].flatten()[1],it[3].flatten()[2],it[3].flatten()[3],
    it[4].flatten()[1],it[4].flatten()[2],it[4].flatten()[3])}
    
    report_original=analysis_org(original_proc,'original')
    report_noise=analysis_noise(noise,'noise')
    report_noise2=analysis_noise2(noise2,'noise2')
    report_noise3=analysis_noise3(noise3,'noise3')
    
    report_original.collectFile(name: 'collected_original_fasta_report.csv', keepHeader: true ,storeDir: "${params.output}/final_report/")
    report_noise.collectFile(name: 'collected_noise_fasta_report.csv', keepHeader: true ,storeDir: "${params.output}/final_report/")
    report_noise2.collectFile(name: 'collected_noise2_fasta_report.csv', keepHeader: true ,storeDir: "${params.output}/final_report/")
    report_noise3.collectFile(name: 'collected_noise3_fasta_report.csv', keepHeader: true ,storeDir: "${params.output}/final_report/")
    
    
    emit: 
    combined

}

