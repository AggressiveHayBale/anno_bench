process contig_splitter {
        
        input: 
        tuple val(name), path(dir)

        output:

        //ADD file not *
        tuple val(name), path("*"), emit: split_ch

        shell:
        '''

        cat !{dir} |  grep -v '^>' | sed -e '1i\\>contig_1' > !{name}_contigsplit.fasta

        line_count = wc -l < !{name}_contigsplit.fasta 

        ct=1
        contig_lines=()
        while [ !{ct} -le 5 ]
        do 
            contig_lines+=(${RANDOM:0:!{line_count}})
            ((ct++)) 
        done
        '''
}