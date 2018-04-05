#!/usr/bin/env nextflow

//////////////////////////////////////////////////////////////////////////////////
// parsing functions

String findValue(String inputStr) {
configFile = file('config_polyCRACKER.txt')
allLines = configFile.readLines()
for( line in allLines ) {
    if(line.startsWith(inputStr)){
        println (line - inputStr - '= ')
        return (line - inputStr - '= ');
    }
}

}


//////////////////////////////////////////////////////////////////////////////////
// parsing functions

vcf_file = params.vcf_file;
tab_file = params.tab_file;
snps_interval = params.snps_interval;
iqtree = params.iqtree.asType(Integer);
work_dir = params.work_dir;


//////////////////////////////////////////////////////////////////////////////////
// vcf to tab

vcf_channel = Channel.fromPath(vcf_file)

process vcf2tab {

    echo false

    input:
    file vcf_channel

    output:
    file tab_channel


    script:
    """python ../../../CNS_commandline.py vcf2tab $vcf_channel tab_channel"""

}



//////////////////////////////////////////////////////////////////////////////////
// partition tab file into chunks

process tab2chunks {

    echo false

    input:
    file tab_channel

    output:
    file fastas_channel

    script:
        """python ../../../CNS_commandline.py tab2chunks $tab_channel ${snps_interval} fastas_channel"""


}


// split fasta chunks into files for tree generation
fastas_channel.splitText()
                .filter {it.toString().size() > 1}
                .set {fastas_check}
fastas_check.map {it -> tuple(it.split())}
            .into {final_fastas}



//////////////////////////////////////////////////////////////////////////////////
// generate trees

process generate_tree {

    echo false
    cpus 3

    input:
        set interval, fasta from final_fastas

    output:
        set interval, file(trees) into final_trees

    script:
        if(iqtree == 1)
        """echo $interval && echo $fasta
           python ../../../CNS_commandline.py run_iqtree --n-threads 2 $fasta trees"""
        else
        """python ../../../CNS_commandline.py run_fasttree $fasta trees"""

}

//////////////////////////////////////////////////////////////////////////////////
// collect the tree writing output

process write_trees2files {

    echo false

    input:
        set interval, file(trees) from final_trees

    output:
        file write_file

    script:
    """python ../../../CNS_commandline.py write_trees_intervals $interval $trees write_file"""

}

write_file.collectFile(name: 'local_trees.nwk', newLine: true)
    .set {local_trees}


process final_output {

    echo false

    input:
    file local_trees

    script:
    """python ../../../CNS_commandline.py local_trees2final_output $local_trees ${work_dir}"""

    }