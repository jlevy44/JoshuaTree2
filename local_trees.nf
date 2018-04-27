#!/usr/bin/env nextflow


//////////////////////////////////////////////////////////////////////////////////
// parsing functions

vcf_file = params.vcf_file;
tab_file = params.tab_file;
snps_interval = params.snps_interval;
phylogeny = params.phylogeny;
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
    """python ../../../JoshuaTree2.py vcf2tab -vcf $vcf_channel -tab tab_channel"""

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
        """python ../../../JoshuaTree2.py tab2chunks -tab $tab_channel -i ${snps_interval} -o fastas_channel"""


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
        """python ../../../JoshuaTree2.py generate_phylogeny -f $fasta -p ${phylogeny} -t trees"""

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