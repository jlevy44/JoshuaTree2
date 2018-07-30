import click, os, sys, subprocess

#######################
#### RUN CLI GROUP ####

CONTEXT_SETTINGS = dict(help_option_names=['-h','--help'], max_content_width=90)

@click.group(context_settings= CONTEXT_SETTINGS)
@click.version_option(version='0.01')
def scaffolder():
    pass

##############################   # fixme just throw this into another python script dedicated to scaffolding pipeline, not worth having in main joshuaTree code. It looks really bad, play it off as old script
#### SCAFFOLD VIA SYNTENY ####

@scaffolder.command()
@click.option('-i', '--scaffolding_inputs_dir', default = './scaffolding_inputs', show_default=True, help='Path containing fasta file one.', type=click.Path(exists=False))
@click.option('-o', '--scaffolding_outputs_dir', default = './scaffolding_outputs', show_default=True, help='Path containing fasta file two.', type=click.Path(exists=False))
@click.option('-n', '--new_genome_name', default = 'Bdist', show_default=True, help='New genome name.', type=click.Path(exists=False))
@click.option('-w', '--weights_file', default = './weights.txt', show_default=True, help='Weights file.', type=click.Path(exists=False))
@click.option('-p', '--primary_proteome_id', default = 314, show_default=True, help='Primary proteome id.', type=click.Path(exists=False))
def scaffold_assemblies(scaffolding_inputs_dir,scaffolding_outputs_dir, new_genome_name, weights_file, primary_proteome_id):
    """Scaffold assemblies based on synteny to references."""
    cwd = os.getcwd()
    scaffolding_bin = os.path.abspath('scaffolding_tool_bin')+'/'
    scaffolding_inputs_dir = os.path.abspath(scaffolding_inputs_dir)
    scaffolding_outputs_dir = os.path.abspath(scaffolding_outputs_dir)
    query_dir = scaffolding_inputs_dir+'/query/'
    references_dir = scaffolding_inputs_dir+'/references/'
    subprocess.call('python %s/renameGenomes.py %s %s'%(scaffolding_bin,new_genome_name, query_dir),shell=True)
    subprocess.call('cp %s/pipelineFunctions.py %s'%(scaffolding_bin,query_dir))
    # add build references and weights file
    os.chdir(query_dir)
    subprocess.call('python %s/scaffoldPipeline.sh --write_sh 1 --version v0 --cds %s '
                    '--gene_name_old %s --build_sample 1 --bb 0 --nuc 0 --com1_2 1 --allmaps 1 --reference_genomes_dir %s'
                    '--output_dir %s --query_genome_dir %s --bin %s --weights_file %s'%(scaffolding_bin,primary_proteome_id, new_genome_name, references_dir,
                                        scaffolding_outputs_dir, query_dir, scaffolding_bin, os.path.abspath(weights_file)), shell=True) # fixme write nextflow config
    os.chdir(cwd)
    """writeSh = params.write_sh.asType(Integer);
buildRef = params.build_ref.asType(Integer);
version = params.version;
CDS = params.cds;
CDSFasta = params.cds_fasta;
geneNameOld = params.gene_name_old;
buildSamp = params.build_sample.asType(Integer);
BB = params.bb.asType(Integer);
nuc = params.nucmer.asType(Integer);
com1_2 = params.com1_2.asType(Integer);
allmaps = params.allmaps.asType(Integer);
nextVersion = 'v' + ((version - 'v').asType(Integer)+1).asType(String);
chanBuildSamples = Channel.fromPath(version + '/*'+version,type: 'dir', relative: true)
workingDir = new File('').getAbsolutePath()
reference_genome_dir = params.reference_genomes_dir
query_genome_dir = params.query_genome_dir
output_dir = params.output_dir
bin = params.bin
""" # fixme add reference and output directories to scaffold pipeline, dockerize scaffold pipeline, make docker images

#### RUN CLI ####

if __name__ == '__main__':
    scaffolder()