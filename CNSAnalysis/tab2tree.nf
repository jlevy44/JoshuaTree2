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

tab_file = params.tab_file;
tree_out = params.tree_out;
n_samples = params.n_samples;
n_bootstrap = params.n_bootstrap.asType(Integer);
tree_final = params.tree_final;
phyml = params.phyml.asType(Integer);
fasta2phylip = params.fasta2phylip;
starting_tree = params.starting_tree;


//////////////////////////////////////////////////////////////////////////////////
// generate bootstrap fasta samples of SNP matrix

tab_bootstraps = Channel.fromPath(tab_file)

process generate_bootstrap_fastas {

    echo false

    input:
    each x from 1..n_bootstrap
    file tab_in from tab_bootstraps

    output:
    file fasta_bootstraps


    script:
    """python ../../../CNS_commandline.py tab2fasta $tab_in ${n_samples} fasta_bootstraps"""

}



//////////////////////////////////////////////////////////////////////////////////
// generate trees from SNP fasta samples

process generate_bootstrap_trees {

    echo false
    cpus 3

    input:
    file fasta_bootstraps

    output:
    file trees

    script:
    if(phyml == 1)
        """python ../../../CNS_commandline.py run_phyml --tree-out trees --fasta2phylip ${fasta2phylip} --starting-tree ${starting_tree} $fasta_bootstraps 1"""
    else
        """python ../../../CNS_commandline.py run_fasttree $fasta_bootstraps trees"""


}



//////////////////////////////////////////////////////////////////////////////////
// generate consensus tree


trees.collectFile(name: 'trees.nh')
    .into {final_tree; support_tree}

process generate_consensus_tree {

    input:
    file trees from final_tree

    script:
    """python ../../../CNS_commandline.py generate_consensus_tree $trees ${tree_out}"""

}

process get_support {

    input:
    file trees from support_tree

    script:
    """python ../../../CNS_commandline.py get_support $trees ${tree_final}"""

}