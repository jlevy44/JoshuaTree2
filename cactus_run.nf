#!/usr/bin/env nextflow


//////////////////////////////////////////////////////////////////////////////////
// parsing functions

cactus_run_files = Channel.from(params.cactus_run_files.split(','));
work_dir = params.work_dir;
environment = params.environment;


//////////////////////////////////////////////////////////////////////////////////
// run cactus
process run_cactus{

    input:
    val cactus_run_files

    script:
    """
    cd ${work_dir}
    set +u
    source ${environment}
    set -u
    sh $cactus_run_files
    deactivate
    """

}