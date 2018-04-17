#!/usr/bin/env nextflow


//////////////////////////////////////////////////////////////////////////////////
// parsing functions

cactus_run_files = Channel.create(params.cactus_run_files.split(","));
work_dir = params.work_dir;
environment = params.environment;


//////////////////////////////////////////////////////////////////////////////////
// run cactus
process run_cactus{

    input:
    file cactus_run_files

    script:
    """
    cd ${work_dir}
    source ${environment}
    sh $cactus_run_files
    """

}