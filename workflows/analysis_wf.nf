include { analysis } from './process/analyis.nf' 

workflow analysis{
    take: 
        combined_annotation_ch 
    main: 
    
    analysis(combined_annotation_ch)


    //emit: 


}

