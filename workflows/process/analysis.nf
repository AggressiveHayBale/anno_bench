process analysis {
    label 'ubuntu'
    storeDir "${params.databases}/analysis"   

    input:


    output:
        path("")
    script:
        """
        r_comparison.R --args 
        """
    stub:
        """

        """
    }