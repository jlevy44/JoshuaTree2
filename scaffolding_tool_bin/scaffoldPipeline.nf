#!/usr/bin/env nextflow

//////////////////////////////////////////////////////////////////////////////////
// parsing functions

params.folder = './'
writeSh = params.write_sh.asType(Integer);
buildRef = params.build_ref.asType(Integer);
version = params.version;
CDS = params.cds;
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
weights_file = params.weights_file
nextDir = file(output_dir+'/'+nextVersion)

result = nextDir.mkdir()
chanBuildRef = Channel.fromPath(reference_genome_dir + '/*',type: 'dir', relative: false)
chanNextSamples = Channel.fromPath(version + '/*'+version,type: 'dir', relative: true)
                        .map { nextVersion + '/' + it.name - version + nextVersion}
                        .map{ file(it) }
chanNextSamples2 = Channel.fromPath(version + '/*'+version,type: 'dir', relative: true)
                        .map { nextVersion + '/' + it.name - version + nextVersion}
                        .map{ file(it) }
                        .subscribe { println it }

linkageChannel = Channel.create()
process writeShFiles {

executor = 'local'

output:
    file 'done' into linkageChannel

script:
if(writeSh)
"""
#!/bin/bash
touch done
module load bedtools/2.25.0
cd ${query_genome_dir}
python ${bin}/writeShFiles.py ${version} ${buildSample} ${buildRef} ${constructSample} ${CDSspecies} ${CDSOld} ${CDSgeneNaming} ${BB} ${nuc} ${weights_file} ${reference_genome_dir} ${query_genome_dir} ${output_dir} ${bin}
"""
else {
    """
    #!/bin/bash
    echo Not writing analysis files...
    touch done
    """
    }
}


linkageChannel = linkageChannel.take(1)
linkageChannel15= Channel.create()
refStr = Channel.create()

process findRef {
executor = 'local'
input:
    file 'done' from linkageChannel

output:
    file 'done' into linkageChannel15
    stdout refStr

script:
"""
#!/usr/bin/python
open('touch','w').close()
with open('${query_genome_dir}/references.txt','r') as f:
    print f.read()
"""
}


linkageChannel15 = linkageChannel15.take(1)

linkageChannel1 = Channel.create()

process buildReference {

clusterOptions = { buildRef == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 6 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
file 'done' from linkageChannel15
each ref from chanBuildRef.toList()

output:
file 'done' into linkageChannel1


script:
if(buildRef)
"""
#!/bin/bash
touch done
cd ${reference_genome_dir}/${ref}
pwd
sh buildRef.sh
"""
else
"""touch done"""
}




linkageChannel1 = linkageChannel1.last()


linkageChannel2 = Channel.create()

process buildSample{
clusterOptions = {buildSamp == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 6 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt' }

input:
each sample from chanBuildSamples.toList()
file 'done' from linkageChannel1

output:
val sample into linkageChannel2

script:
if(buildSamp)
"""
#!/bin/bash
touch done
cd ${query_genome_dir}/${version}/${sample}
pwd
sh build.sh
"""
else
"""#!/bin/bash
touch done"""
}
//linkageChannel2 = linkageChannel2.take(1)

linkageChannel3 = Channel.create()
process nucmerfy {

clusterOptions = {nuc == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 6 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt'}
input:
val sample from linkageChannel2



output:
val sample into linkageChannel3

script:
if(nuc)
"""
#!/bin/bash
touch done
cd ${query_genome_dir}/${version}/${sample}
sh nucCommand.sh
"""
else
"""touch done"""
}


linkageChannel3_5 = Channel.create()
process BBfy {

clusterOptions = {BB == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 8 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt'}
input:
val sample from linkageChannel3



output:
val sample into linkageChannel3_5


script:
if(BB)
"""
#!/bin/bash
export _JAVA_OPTIONS='-Xmx15G'
touch done
module load bedtools/2.25.0
cd ${query_genome_dir}/${version}/${sample}
sh BB_build.sh > nohup.out
"""
else
"""touch done"""
}


linkageChannel4 = Channel.create()
linkageChannel5 = Channel.create()
refStrList = refStr.take(1)
    .map { it - '\n' }
    .toList()
process com_1_2 {
clusterOptions = { com1_2 == 1 ? '-P plant-analysis.p -cwd -q normal.q -pe pe_slots 3 -e OutputFile.txt' : '-P plant-analysis.p -cwd -l high.c -pe pe_slots 1 -e OutputFile.txt'}

input:
val sample from linkageChannel3_5
each refText from refStrList

output:
val sample into linkageChannel4


script:
if(com1_2)
"""
#!/bin/bash
touch done
module load bedtools/2.25.0
cd ${query_genome_dir}
python ${bin}/com1_2.py ${sample} ${refText} ${CDS} ${version} ${reference_genome_dir}
"""
else
"""touch done"""
}
process Allmaps {
clusterOptions = '-P plant-analysis.p -cwd -l h_rt=24:00:00 -pe pe_slots 16 -e OutputFile.txt'
queue = 'long'
errorStrategy = 'ignore'

input:
val sample from linkageChannel4


output:
val sample into linkageChannel5



script:
if(allmaps)
"""
#!/bin/bash
touch done
cd ${query_genome_dir}/${version}/${sample}
sh qsub_build.sh
"""
else
    error "Unable to build new genome: ${sample}"
}
