include { analysis } from './process/analyis.nf' 

workflow analysis{
    take: 
        combined_annotation_ch 
    main: 

    //WIP filter channel to only contain original, noise, split for analysis 
    combined_annotation_ch.filter(

    analysis(combined_annotation_ch)


    //emit: 


}

