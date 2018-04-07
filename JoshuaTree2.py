#################
#### IMPORTS ####

import click, os, pandas as pd, numpy as np, subprocess
from pybedtools import BedTool
from pyfaidx import Fasta
from itertools import combinations
import scipy.sparse as sps
import glob, re
from random import randint
from Bio import SeqIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
import multiprocessing as mp
from time import sleep

#################
#### CLASSES ####


#########################
#### SYNTENY CLASSES ####

class SuperSynteny:
    def __init__(self, pairwise_syntenies, max_distance, q_genome):
        self.pairwise_syntenies = {synteny.synteny_file: synteny for synteny in pairwise_syntenies}
        self.max_distance = max_distance
        self.q_genome = q_genome

    def generate_global_synteny_graph(self, fasta_out_dir):
        global_synteny_indices = []
        for synteny in self.pairwise_syntenies.values():
            synteny.synteny_structure['name'] = np.array(synteny.synteny_structure.index)#synteny.synteny_structure['q_chr'] + '\t' + synteny.synteny_structure['q_xi'] + '\t' + synteny.synteny_structure['q_xf']#[['q_chr','q_xi','q_xf']]
            global_synteny_indices.extend([(synteny.synteny_file,name) for name in synteny.synteny_structure['name'].as_matrix().tolist()])
        global_synteny_indices = pd.MultiIndex.from_tuples(global_synteny_indices,names=['synteny_file','region'])
        global_synteny_matrix = pd.DataFrame(np.zeros((len(global_synteny_indices),)*2),index=global_synteny_indices,columns=global_synteny_indices)
        global_synteny_matrix = global_synteny_matrix.sort_index(axis=0).sort_index(axis=1)
        referenced_regions = list(global_synteny_matrix.index)
        regions = np.array(referenced_regions)[:,1].tolist()
        regions_bed_files = dict(zip(regions,map(lambda x: BedTool(x,from_string=True),regions)))
        for region_i,region_j in combinations(regions,r=2) + zip(regions,regions):
            if region_i == region_j:
                    global_synteny_matrix.loc[pd.IndexSlice[:,region_i],pd.IndexSlice[:,region_i]] = 1
            else:
                distance = int(str(regions_bed_files[region_i].closest(regions_bed_files[region_j],d=True)).split('\t')[-1])
                if distance >=0 and distance <= self.max_distance:
                    global_synteny_matrix.loc[pd.IndexSlice[:,region_i],pd.IndexSlice[:,region_j]] = 1
                    global_synteny_matrix.loc[pd.IndexSlice[:,region_j],pd.IndexSlice[:,region_i]] = 1
        n_components, connected_components = sps.csgraph.connected_components(global_synteny_matrix.as_matrix())
        referenced_regions = np.array(referenced_regions)
        for i,connected_component_list in enumerate(connected_components):
            out_fasta = []
            q_regions = BedTool('\n'.join(np.vectorize(lambda region: region + '%s.%s'%(self.q_genome.protID,'_'.join(region)))(str(BedTool(map(lambda x: x.split(),referenced_regions[connected_component_list,1].tolist())).sort().merge(d=100000000000)).splitlines())),from_string=True).sequence(fi=self.q_genome.fasta,name=True)#np.vectorize(lambda region: region + '%s %s'%(self.q_genome.protID,' '.join(region)))(referenced_regions[connected_component_list]).tolist()
            out_fasta.append(q_regions)
            referenced_regions_s_species = referenced_regions[connected_component_list,:]
            for species in set(referenced_regions_s_species[:,0]):
                s_regions = BedTool('\n'.join(np.vectorize(lambda region: region + '%s.%s'%(self.pairwise_syntenies[species].s_genome.protID,'_'.join(region)))(str(BedTool(map(lambda x: x.split(),np.vectorize(lambda x: '\t'.join(self.pairwise_syntenies[species].synteny_structure.loc[x,['s_chr','s_xi','s_xf']].as_matrix().tolist()))(referenced_regions_s_species[referenced_regions_s_species[:,0]==species,1]).tolist())).sort().merge(d=self.max_distance)).splitlines())),from_string=True).sequence(fi=self.pairwise_syntenies[species].s_genome.fasta,name=True)
                out_fasta.append(s_regions)
            with open(fasta_out_dir+'/'+'fasta_output%d.fasta'%i,'w') as f:
                f.write(reduce(lambda x,y: x+'\n'+y,out_fasta))

class PairwiseSynteny:
    def __init__(self, q_genome, s_genome, synteny_file='', loci_threshold = 4):
        self.synteny_file = synteny_file
        self.q_genome = q_genome
        self.s_genome = s_genome
        self.loci_threshold = loci_threshold

    def generate_synteny_structure(self,synteny_path):
        """Take anchor file or synteny file and searches for starting and ending genes for each syntenic block"""
        if self.synteny_file.endswith('.unout'):
            self.synteny_structure = self.unout2structure(self.q_genome, self.s_genome)
        elif self.synteny_file.endswith('.anchors'):
            self.synteny_structure = self.anchor2structure(self.q_genome, self.s_genome)
        elif self.synteny_file.endswith('.bed'):
            self.import_synteny_structure()
        else:
            self.run_synteny(self.q_genome,self.s_genome,synteny_path)
            self.synteny_structure = self.anchor2structure(self.q_genome, self.s_genome)

    def import_synteny_structure(self):
        self.synteny_structure = pd.read_table(self.synteny_file,header=None,names=['q_chr','q_xi','q_xf','s_chr','s_xi','s_xf'])
        self.synteny_structure = self.synteny_structure.rename(dict(enumerate((self.synteny_structure['q_chr']+'\t'+self.synteny_structure['q_xi']+'\t'+self.synteny_structure['q_xf']).as_matrix().tolist())))

    def unout2structure(self,q_genome,s_genome):
        with open(self.synteny_file,'r') as f:
            lines = np.array(f.read().splitlines())
        anchors = np.array_split(lines,np.where(np.vectorize(lambda line: line.startswith('\t') == 0)(lines))[0])
        synteny_structure = []
        for anchor in anchors:
            anchor = np.array(map(lambda line: line.split('\t')[1].split(','),anchor[1:].tolist()))
            if len(anchor) >= self.loci_threshold and anchor.tolist():
                q_genes, s_genes = anchor[:,2], anchor[:,5]
                q_coords, s_coords = q_genome.df[q_genes,:], s_genome.df[s_genes,:]
                synteny_structure.append([q_coords.iloc[0,0],q_coords[['xi','xf']].min(),q_coords[['xi','xf']].max(),s_coords.iloc[0,0],s_coords.loc[['xi','xf']].min(),s_coords.loc[['xi','xf']].max()])
        self.synteny_structure = pd.DataFrame(synteny_structure,columns=['q_chr','q_xi','q_xf','s_chr','s_xi','s_xf'],index = np.vectorize(lambda x: '\t'.join(map(str,x[:3])))(synteny_structure))

    def run_synteny(self,genome1,genome2, synteny_path):
        pwd = os.getcwd()
        os.chdir(synteny_path)
        for abs_path, link_name in zip([genome1.bed_file,genome2.bed_file,genome1.CDS_file,genome2.CDS_file],[genome1.protID+'.bed',genome2.protID+'.bed',genome1.protID+'.cds',genome2.protID+'.cds']):
            subprocess.call('ln -s %s %s'%(abs_path,link_name),shell=True)
        subprocess.call('python -m jcvi.compara.catalog ortholog %s %s'%(genome1.short_name,genome2.short_name),shell=True)
        if genome1.short_name != genome1.protID and genome2.short_name != genome2.protID:
            subprocess.call('mv %s.%s.lifted.anchors %s.%s.lifted.anchors'%(genome1.short_name,genome2.short_name,genome1.protID,genome2.protID),shell=True)
        self.synteny_file = os.path.abspath('%s.%s.lifted.anchors'%(genome1.protID,genome2.protID))
        os.chdir(pwd)

    def anchor2structure(self, q_genome, s_genome):
        anchor_file = self.synteny_file
        with open(anchor_file,'r') as f:
            anchors = f.read().split('###')[1:]
        synteny_structure = []
        for anchor in anchors:
            #print anchor
            if anchor:
                genes = np.array([line.split()[:2] for line in anchor.splitlines() if line])
                if genes.shape[0] >= self.loci_threshold:
                    q_genes, s_genes = genes[:,0] , genes[:,1]
                    q_coords, s_coords = q_genome.df.loc[q_genes,:], s_genome.df.loc[s_genes,:]
                    #print q_coords[['xi','xf']]
                    synteny_structure.append([q_coords.iloc[0,0],q_coords[['xi','xf']].values.min(),q_coords[['xi','xf']].values.max(),s_coords.iloc[0,0],s_coords[['xi','xf']].values.min(),s_coords[['xi','xf']].values.max()])
        self.synteny_structure = pd.DataFrame(synteny_structure,columns=['q_chr','q_xi','q_xf','s_chr','s_xi','s_xf'])

    def synteny_structure_2_bed(self,filename):
        self.synteny_structure.to_csv(filename,sep='\t',index=False,header=None)

    def synteny_structure_2_link(self, filename):
        df = self.synteny_structure
        df['q_chr'] = np.vectorize(lambda x: self.q_genome.protID+'-'+x)(df['q_chr'])
        df['s_chr'] = np.vectorize(lambda x: self.s_genome.protID+'-'+x)(df['s_chr'])
        df.to_csv(filename,sep=' ',index=False,header=None)
        self.link = os.path.abspath(filename)

    def export_karyotypes(self,circos_input):
        self.s_genome.export_karyotype(circos_input+'/'+self.s_genome.protID+'.karyotype.txt')
        self.q_genome.export_karyotype(circos_input+'/'+self.q_genome.protID+'.karyotype.txt')

######################
#### GENOME CLASS ####

class Genome:
    def __init__(self, fasta_file, bed_file, protID, gff_file = '', gene_info = 'Name'):
        self.fasta_file = fasta_file
        self.fasta = Fasta(fasta_file)
        self.gene_info = gene_info
        self.bed_file = os.path.abspath(bed_file)
        self.short_name = self.bed_file.split('/')[-1].replace('.bed3','').replace('.bed','')
        self.protID = protID
        if gff_file and os.path.exists(self.bed_file) == 0 or (os.path.exists(self.bed_file) and os.stat(self.bed_file).st_size == 0):
            self.gff_file = gff_file
            #click.echo('python -m jcvi.formats.gff bed --type=mRNA --key=%s %s > %s'%(self.gene_info,self.gff_file,self.bed_file))
            subprocess.call('python -m jcvi.formats.gff bed --type=mRNA --key=%s %s > %s'%(self.gene_info,self.gff_file,self.bed_file),shell=True)
            """
            with open(gff_file,'r') as f:
                for header_line,line in enumerate(f):
                    if line.startswith('#') == 0:
                        break
            df = pd.read_table(gff_file, skiprows=header_line, header=None,names=['chr','rm_1','feature','xi','xf','rm_3','rm_4','rm_5','Annotation'])
            df = df[df['feature'] == 'mRNA'].drop([feat for feat in list(df) if 'rm' in feat],axis=1)
            df = df[np.vectorize(lambda line: 'longest=1' in line)(df['Annotation']).astype(bool)]
            df['Gene'] = np.vectorize(lambda line: line.split(';')[1].replace('Name=',''))(df['Annotation'])
            df['xi'] -= 1
            df = df.drop(['feature','Annotation'],axis=1).reindex(columns=['chr','xi','xf','Gene'])

            self.df = df
            """
        self.df = pd.read_table(self.bed_file,header=None,names=['chr','xi','xf','Gene'],dtype={'chr':str,'xi':np.int,'xf':np.int,'Gene':str},usecols=[0,1,2,3])
        self.df = self.df.set_index('Gene')

    def export_bed(self,filename):
        df = self.df.reset_index().rename(dict(index='Gene'),axis='columns').reindex(columns=['chr','xi','xf','Gene'])
        df.to_csv(filename,sep='\t',index=False,header=None)

    def extract_CDS(self):
        self.CDS_file = self.bed_file.replace('.bed3','.cds').replace('.bed','.cds')
        subprocess.call('python -m jcvi.formats.gff load %s %s --parents=mRNA --children=CDS --id_attribute=%s -o %s'%(self.gff_file,self.fasta_file,self.gene_info,self.CDS_file),shell=True)

    def export_karyotype(self, filename, n_chromosomes=25, shorten_chr=False):
        df = pd.read_table(self.fasta_file+'.fai', header=None,names=['chr','length'],usecols=[0,1],dtype=dict(zip(['chr','length'],[str,np.int])))
        df = df.sort_values(['length'],ascending=False)
        if n_chromosomes < df.shape[0]:
            df = df.iloc[:n_chromosomes,:].reset_index(drop=True)
        out_txt = []
        for i in range(df.shape[0]):
            chrom = df.loc[i,'chr']
            chr_name = chrom if not shorten_chr else chrom[0] + chrom.split('_')[-1]
            if i >= 25:
                out_txt.append('chr - %s-%s %s 0 %d %d,%d,%d\n'%(self.protID,chrom,chr_name,df.loc[i,'length'],randint(1,255),randint(1,255),randint(1,255)))
            else:
                out_txt.append('chr - %s-%s %s 0 %d chr%d\n'%(self.protID,chrom,chr_name,df.loc[i,'length'],i+1))
        with open(filename,'w') as f:
            f.writelines(out_txt)
        self.karyotype = os.path.abspath(filename)

######################
#### CIRCOS CLASS ####

class Circos:
    def __init__(self,PairwiseSynteny):
        self.synteny = PairwiseSynteny

    def write_ideogram_config(self, filename='txideogram.conf'):
        with open(filename,'w') as f:
            f.write("""<ideogram>
                show = yes
                <spacing>
                default = 0.005r
                </spacing>
                radius    = 0.9r
                thickness = 40p
                fill      = yes
                show_label = yes
                label_font = default
                label_radius = 1.08r
                label_size = 40
                label_parallel = yes
                show_bands = yes
                fill_bands = yes
                </ideogram>""")
        self.ideogram = filename
        return filename

    def write_ticks_config(self, filename='txticks.conf'):
        with open(filename,'w') as f:
            f.write("""show_ticks = yes
                show_tick_labels = yes
                <ticks>
                radius = 1.01r
                color = black
                thickness = 2p
                multiplier = 1e-6
                format = %d
                <tick>
                spacing = 1u
                size = 5p
                </tick>
                <tick>
                spacing = 5u
                size = 10p
                show_label = yes
                label_size = 20p
                label_offset = 10p
                format = %d
                </tick>
                </ticks>""")
        self.ticks = filename
        return filename

    def generate_config(self, ticks = 'txticks.conf', ideogram = 'txideogram.conf', links_and_rules = 'linksAndrules.conf', config='circos.conf'):
        colors = pd.read_table(self.synteny.s_genome.karyotype,header=None,usecols=[2,6]).as_matrix()
        self.links_and_rules = links_and_rules
        self.config = config
        if hasattr(self, 'ticks'):
            self.write_ticks_config(ticks)
        if hasattr(self, 'ideogram'):
            self.write_ideogram_config(ideogram)
        with open(self.config,'w') as f:
            f.write("""# circos.conf
                karyotype = %s, %s
                chromosomes_units = 1000000
                chromosomes_display_default = yes
                <<include %s>>
                <<include %s>>
                <<include %s>>
                <image><<include etc/image.conf>></image>
                <<include etc/colors_fonts_patterns.conf>>
                <<include etc/housekeeping.conf>>
                """%(self.synteny.q_genome.karyotype,self.synteny.s_genome.karyotype,self.ideogram,self.ticks,self.links_and_rules))
        with open(self.links_and_rules,'w') as f:
            f.write("""
                <links>
                <link>
                file = %s
                radius = 0.99r
                bezier_radius = 0r
                ribbon = yes
                color = black_a4
                <rules>
                <rule>
                condition = var(intrachr)
                show = no
                </rule>\n"""%(self.synteny.link) + '\n'.join(['<rule>\ncondition = to(%s)\ncolor = %s\n</rule>'%(chrom,color) for chrom,color in colors.itertuples()]) + '\n</rules>\n</link>\n</links>')

    def run_circos(self, output_dir='./', pdf=False):
        subprocess.call('circos -conf %s -outputfile %s-%s -outputdir %s'%(self.config,self.synteny.q_genome.protID,self.synteny.s_genome.protID,output_dir),shell=True)
        if pdf:
            subprocess.call('convert %s/%s-%s.png %s/%s-%s.pdf'%(os.path.abspath(output_dir),self.synteny.q_genome.protID,self.synteny.s_genome.protID,os.path.abspath(output_dir),self.synteny.q_genome.protID,self.synteny.s_genome.protID),shell=True)

##########################
#### CACTUS RUN CLASS ####

class CactusRun:
    def __init__(self,fasta_output_path,cactus_run_directory,cactus_softlink, cactus_env_softlink, nickname_file = '', fasta_path = ''):
        self.fasta_output_path = fasta_output_path +'/'
        self.cactus_run_directory = cactus_run_directory +'/'
        self.cactus_output = self.cactus_run_directory+'output/'
        self.hal_path = self.cactus_output+'hal/'
        self.cactus_softlink = os.path.abspath(cactus_softlink)
        if not nickname_file:
            self.nickname_file = self.cactus_run_directory + 'prot_dict'
            self.protIDs = [fasta.split('_')[-2] for fasta in glob.glob(fasta_path+'/*.fa')+glob.glob(fasta_path+'/*.fasta')]
            with open(self.nickname_file,'w') as f:
                f.write('\n'.join(['\t'.join((protID,)*2) for protID in self.protIDs]))
        self.nickname_file = nickname_file
        subprocess.call('source %s'%cactus_env_softlink,shell=True)

    def write_fastas_seqfile(self):
        with open(self.nickname_file,'r') as f:
            self.nicknames = dict([tuple(line.split()) for line in f.read().splitlines() if line])
        fastas = glob.glob(self.fasta_output_path+'/*.fa')+glob.glob(self.fasta_output_path+'/*.fasta')
        self.run_files = []
        for fasta in fastas:
            self.nickname2file_dict = {}
            self.fasta_run_dir = self.cactus_output+fasta[fasta.rfind('/')+1:].split('.')[0]+'/'
            try:
                os.mkdir(self.fasta_run_dir)
            except:
                pass
            subprocess.call('rm %s/*.fa %s/*.fasta -r'%(self.fasta_run_dir,self.fasta_run_dir),shell=True)
            with open(fasta) as f1:
                records = SeqIO.parse(f1,'fasta')
                for record in records:
                    protID = record.id.split('.')[0]
                    nickname = self.nicknames[protID]
                    record.id = record.id.replace(protID,nickname)
                    record.description = record.id
                    SeqIO.write(record, open(self.fasta_run_dir+nickname+'.fa','a'),'fasta')
                    self.nickname2file_dict[nickname] = os.path.abspath(self.fasta_run_dir+nickname+'.fa')
            self.generate_trees()
            self.seqfile = self.fasta_run_dir+'seqfile'
            with open(self.fasta_run_dir+'output_tree.nh','r') as f1, open(self.seqfile,'w') as f2:
                f2.write(f1.read()+'\n'+'\n'.join(['%s %s'%(nickname, nickname_file) for nickname, nickname_file in self.nickname2file_dict.items()]))
            run_file = os.path.abspath(self.seqfile+'.sh')
            self.run_files.append(run_file)
            with open(run_file,'w') as f:
                f.write('#!/bin/bash\nexport _JAVA_OPTIONS="-Xmx155g"\n%s --maxThreads 16 %s %s %s >& %s\nscp %s %s'%(os.path.abspath(self.cactus_output),self.seqfile,self.fasta_run_dir,fasta.replace('.fa','.hal'),self.seqfile+'.sh.stdout',fasta.replace('.fa','.hal'),self.hal_path+fasta.split('/')[-1].replace('.fa','.hal')))

    def run_cactus(self):
        for run_file in self.run_files:
            subprocess.call('nohup sh %s &'%(run_file),shell=True)
    # fixme, make sure to activate progressive cactus environment; sourcec environment in cactus folder, add hal2maf, send maf 2 cns analysis, parallelize and pipeline all scripts, maybe call nextflow from within python for each of the jobs for pipeline submission

    def generate_trees(self,scaled = 10000, kmer_length = 23, multi_fasta = False):
        # fixme use mafft to generate new guide trees for each set of input fasta files
        subprocess.call('sourmash compute -f --scaled %d %s -o %s.sig -k %d %s'%(scaled,self.fasta_run_dir,self.fasta_run_dir+'tree',kmer_length,'--singleton' if multi_fasta else ''),shell=True)
        subprocess.call('sourmash compare %s.sig --csv %s.cmp.csv'%(self.fasta_run_dir+'tree',self.fasta_run_dir+'tree'),shell=True)
        df = pd.read_csv('%s.cmp.csv'%self.fasta_run_dir+'tree',index_col = None)
        samples = [fasta.split('/')[-1].replace('.fa','') for fasta in list(df)]
        distance_matrix = df.as_matrix()
        constructor = DistanceTreeConstructor()
        dm = _DistanceMatrix(names=samples,matrix=[list(distance_matrix[i,0:i+1]) for i in range(len(samples))])
        tree = constructor.nj(dm)
        Phylo.write(tree,self.fasta_run_dir+'output_tree.nh','newick')

    def hal2maf(self, n_cpus, hal2maf_softlink):
        self.hal2maf_softlink = hal2maf_softlink
        self.maf_path = self.hal_path.replace('hal','maf')
        run_hal = lambda hal: self.run_hal2maf(hal)
        for hal in glob.glob(self.hal_path+'/*.hal'):
            proc = mp.Process(target=run_hal, args=(hal,))
            proc.daemon = True
            proc.start()
            while len(mp.active_childern()) > n_cpus:
                sleep(1)
        while len(mp.active_childern()) > 0:
            sleep(1)

    def run_hal2maf(self,hal):
        os.system('%s %s %s'%(self.hal2maf_softlink, hal, hal.replace('hal','maf')))

###################
#### MAF CLASS #### fixme try to inheret from circos class, maf class will allow you to extract CNS or VCF, visualize CS and CNS

###################
#### SNP CLASS #### fixme vcf or tab file, as well as df, can local pca, local tree, final tree (run iqtree), visualize tree, produce enw vcf files and interact with other snp objects
# fixme maybe add nextflow class?
# fixme add test classes

###################
#### CNS CLASS ####

#######################
#### RUN CLI GROUP ####

CONTEXT_SETTINGS = dict(help_option_names=['-h','--help'], max_content_width=90)

@click.group(context_settings= CONTEXT_SETTINGS)
@click.version_option(version='0.2')
def joshuatree():
    pass

##################
#### COMMANDS ####

#####################
#### RUN SYNTENY ####

@joshuatree.command() #fixme can use a combo of unout and anchor files, or generate new anchor files after reach end of synteny file list, but there should be more
@click.option('-q', '--query_protID', default = 'xxx', show_default=True, help='Three number or letter proteome identifier of strain being compared against.')
@click.option('-fi', '--fasta_path', default = './fasta_path/', show_default=True, help='Fasta path containing all of the input genomes. Genome naming must conform to xxx_[protID]_xxx.[fa/fasta].', type=click.Path(exists=False))
@click.option('-s', '--synteny_path', default = './synteny_path/', show_default=True, help='Path containing synteny files, .unout or .anchors files. *.unout must conform to following pattern: [PAC4GC/PAC2_0].[q_protID]-[PAC4GC/PAC2_0].[s_protID]_5.unout; *.anchors must conform to: [q_protID].[s_protID].[*].anchors. Not neccessary to add files to this path, synteny will be generated if no specification.', type=click.Path(exists=False))
@click.option('-gff', '--gff_path', default = './gff_path/', show_default=True, help='Gff path containing all of the gff/gff3 files. Gff naming must conform to: xxx.[protID].[gff/gff3].', type=click.Path(exists=False))
@click.option('-bed', '--bed_path', default = './bed_path/', show_default=True, help='Bed path containing all of the bed files.', type=click.Path(exists=False))
@click.option('-info', '--gene_info', default = 'Name', show_default=True, help='Naming convention for gff file\'s gene name field.', type=click.Choice(['Name', 'gene_name']))
@click.option('-fo', '--fasta_out_dir', default = './fasta_output/', show_default=True, help='Path containing all syntenic aligned regions, organized into fasta files for multiple sequence alignment via Cactus.', type=click.Path(exists=False))
@click.option('-bp', '--bps_threshold', default= 100000, show_default=True, help='Maximum distance from which to merge nearby syntenic regions of a particular genome.')
@click.option('-l', '--loci_threshold', default= 4, show_default=True, help='Minimum number of genes in a syntenic block in order to include the block.')
@click.option('-c', '--circos', is_flag=True, help='If selected, visualize each pairwise alignment using Circos.')
@click.option('-ci', '--circos_inputs', default = './circos_inputs/', show_default=True, help='Path containing all of circos inputs and configuration files.', type=click.Path(exists=False))
@click.option('-co', '--circos_outputs', default = './circos_outputs/', show_default=True, help='Path containing all of circos output images.', type=click.Path(exists=False))
def run_synteny_pipeline(query_protID,fasta_path,synteny_path,gff_path, bed_path, gene_info, fasta_out_dir, bps_threshold, loci_threshold, circos, circos_inputs, circos_outputs):
    """Stitch together many pairwise syntenic blocks into multiple species syntenic blocks and output as fasta files for a multiple sequence alignment. If synteny files are not supplied, conduct pairwise synteny between all included strains."""
    fasta_files = {fasta.split('_')[-2] : fasta for fasta in glob.glob(fasta_path+'/*.fa')+glob.glob(fasta_path+'/*.fasta')}
    gff_files = {gff.split('.')[-2] : gff for gff in glob.glob(gff_path+'/*.gff')+glob.glob(gff_path+'/*.gff3') }
    intersect_keys = set(fasta_files.keys()) & set(gff_files.keys())
    fasta_files = {protID:fasta for protID,fasta in fasta_files.items() if protID in intersect_keys}
    gff_files = {protID:gff for protID,gff in gff_files.items() if protID in intersect_keys}
    genomes = {}
    pairwise_syntenies = []
    for protID in gff_files:
        genomes[protID] = Genome(fasta_files[protID],bed_path+'/'+protID+'.bed',protID,gff_files[protID],gene_info)
        if circos:
            genomes[protID].export_karyotype(circos_inputs+'/'+protID+'.karyotype.txt')
    synteny_files = glob.glob(synteny_path+'/*.unout')+glob.glob(synteny_path+'/*.anchors')
    if synteny_files:
        for synteny_file in os.listdir(synteny_path):
            if synteny_file.endswith('.unout'):
                coords = reduce(lambda x,y: x+y, sorted([[m.start(0),m.end(0)] for m in re.finditer('PAC2_0|PAC4GC',synteny_file)]))[1::2]
                q_protID, s_prot_ID = map(lambda x: synteny_file[x+1:x+4],coords)
            else:
                q_protID, s_prot_ID = tuple(synteny_file.split('.')[:2])
            pairwise_synteny = PairwiseSynteny(genomes[q_protID],genomes[s_prot_ID],synteny_file,loci_threshold=loci_threshold)
            pairwise_synteny.generate_synteny_structure(synteny_path)
            pairwise_syntenies.append(pairwise_synteny)
    else:
        for protID in gff_files:
            genomes[protID].extract_CDS()
        for s_prot_ID in set(genomes.keys()) - {query_protID}:
            pairwise_synteny = PairwiseSynteny(genomes[query_protID],genomes[s_prot_ID],loci_threshold=loci_threshold)
            pairwise_synteny.generate_synteny_structure(synteny_path)
            pairwise_syntenies.append(pairwise_synteny)

    if circos:
        for pairwise_synteny in pairwise_syntenies:
            pairwise_synteny.synteny_structure_2_link(circos_inputs+'/%s.%s.link.txt'%(pairwise_synteny.q_genome.protID,pairwise_synteny.s_genome.protID))
            circos_obj = Circos(pairwise_synteny)
            circos_obj.generate_config(ticks = circos_inputs+'/txticks.conf', ideogram = circos_inputs+'/txideogram.conf', links_and_rules = circos_inputs+'/linksAndrules.conf', config=circos_inputs+'/circos.conf')
            circos_obj.run_circos(circos_outputs+'/')

    super_synteny = SuperSynteny(pairwise_syntenies, bps_threshold, genomes[query_protID])
    super_synteny.generate_global_synteny_graph(fasta_out_dir)

####################
#### RUN CIRCOS ####

@joshuatree.command()
@click.option('-fi', '--fasta_path', default = './fasta_path/', show_default=True, help='Fasta path containing all of the input genomes. Genome naming must conform to xxx_[protID]_xxx.[fa/fasta].', type=click.Path(exists=False))
@click.option('-gff', '--gff_path', default = './gff_path/', show_default=True, help='Gff path containing all of the gff/gff3 files. Gff naming must conform to: xxx.[protID].[gff/gff3].', type=click.Path(exists=False))
@click.option('-s', '--synteny_path', default = './synteny_path/', show_default=True, help='Path containing synteny files, .unout or .anchors files. *.unout must conform to following pattern: [PAC4GC/PAC2_0].[q_protID]-[PAC4GC/PAC2_0].[s_protID]_5.unout; *.anchors must conform to: [q_protID].[s_protID].[*].anchors. Not neccessary to add files to this path, synteny will be generated if no specification.', type=click.Path(exists=False))
@click.option('-bed', '--bed_path', default = './bed_path/', show_default=True, help='Bed path containing all of the bed files.', type=click.Path(exists=False))
@click.option('-ci', '--circos_inputs', default = './circos_inputs/', show_default=True, help='Path containing all of circos inputs and configuration files.', type=click.Path(exists=False))
@click.option('-co', '--circos_outputs', default = './circos_outputs/', show_default=True, help='Path containing all of circos output images.', type=click.Path(exists=False))
@click.option('-l', '--loci_threshold', default= 4, show_default=True, help='Minimum number of genes in a syntenic block in order to include the block.')
@click.option('-info', '--gene_info', default = 'Name', show_default=True, help='Naming convention for gff file\'s gene name field.', type=click.Choice(['Name', 'gene_name']))
@click.option('-n', '--n_cpus', default = 16, show_default=True, help='Number of cpus used to convert hal 2 maf files.')
def circos_dropper(fasta_path, gff_path, synteny_path, bed_path, circos_inputs, circos_outputs, loci_threshold,gene_info, n_cpus):
    """Visualize many pairwise synteny results. If synteny files are not supplied, conduct pairwise synteny between all included strains."""
    fasta_files = {fasta.split('_')[-2] : fasta for fasta in glob.glob(fasta_path+'/*.fa')+glob.glob(fasta_path+'/*.fasta')}
    gff_files = {gff.split('.')[-2] : gff for gff in glob.glob(gff_path+'/*.gff')+glob.glob(gff_path+'/*.gff3') }
    intersect_keys = set(fasta_files.keys()) & set(gff_files.keys())
    fasta_files = {protID:fasta for protID,fasta in fasta_files.items() if protID in intersect_keys}
    gff_files = {protID:gff for protID,gff in gff_files.items() if protID in intersect_keys}
    genomes = {}
    pairwise_syntenies = []
    print gff_files, fasta_files
    for protID in gff_files:
        genomes[protID] = Genome(fasta_files[protID],bed_path+'/'+protID+'.bed',protID,gff_files[protID],gene_info)
        genomes[protID].export_karyotype(circos_inputs+'/'+protID+'.karyotype.txt')
    print genomes
    synteny_files = glob.glob(synteny_path+'/*.unout')+glob.glob(synteny_path+'/*.anchors')
    if synteny_files:
        for synteny_file in os.listdir(synteny_path):
            if synteny_file.endswith('.unout'):
                coords = reduce(lambda x,y: x+y, sorted([[m.start(0),m.end(0)] for m in re.finditer('PAC2_0|PAC4GC',synteny_file)]))[1::2]
                q_protID, s_prot_ID = map(lambda x: synteny_file[x+1:x+4],coords)
            else:
                q_protID, s_prot_ID = tuple(synteny_file.split('.')[:2])
            pairwise_synteny = PairwiseSynteny(genomes[q_protID],genomes[s_prot_ID],synteny_file,loci_threshold=loci_threshold)
            pairwise_synteny.generate_synteny_structure(synteny_path)
            pairwise_syntenies.append(pairwise_synteny)
    else:
        for protID in gff_files:
            proc = mp.Process(target=lambda: genomes[protID].extract_CDS, args=None)
            proc.daemon = True
            proc.start()
            while len(mp.active_childern()) > n_cpus:
                sleep(1)
        while len(mp.active_childern()) > 0:
            sleep(1)

        pairwise_syntenies = []

        def mycallback(x):
            pairwise_syntenies.extend(x)

        def p_synteny(protIDs):
            q_protID, s_prot_ID = protIDs
            pairwise_synteny = PairwiseSynteny(genomes[q_protID],genomes[s_prot_ID],loci_threshold=loci_threshold)
            pairwise_synteny.generate_synteny_structure(synteny_path)
            return pairwise_synteny

        """
        for q_protID, s_prot_ID in combinations(genomes.keys(),r=2):
            pairwise_synteny = PairwiseSynteny(genomes[q_protID],genomes[s_prot_ID],loci_threshold=loci_threshold)
            pairwise_synteny.generate_synteny_structure(synteny_path)
            pairwise_syntenies.append(pairwise_synteny)"""
        r = mp.Pool(n_cpus).map_async(p_synteny,combinations(genomes.keys(),r=2),callback=mycallback)
        r.wait()


    for pairwise_synteny in pairwise_syntenies:
        pairwise_synteny.synteny_structure_2_link(circos_inputs+'/%s.%s.link.txt'%(pairwise_synteny.q_genome.protID,pairwise_synteny.s_genome.protID))
        circos_obj = Circos(pairwise_synteny)
        circos_obj.generate_config(ticks = circos_inputs+'/txticks.conf', ideogram = circos_inputs+'/txideogram.conf', links_and_rules = circos_inputs+'/linksAndrules.conf', config=circos_inputs+'/circos.conf')
        circos_obj.run_circos(circos_outputs+'/')

####################
#### RUN CACTUS ####

@joshuatree.command()
@click.option('-fo', '--fasta_output_path', default = './fasta_output/', show_default=True, help='Path containing all syntenic aligned regions, organized into fasta files for multiple sequence alignment via Cactus.', type=click.Path(exists=False))
@click.option('-c', '--cactus_run_directory', default = './cactus_run/', show_default=True, help='Path containing all syntenic aligned regions, organized into fasta files for multiple sequence alignment via Cactus.', type=click.Path(exists=False))
@click.option('-cac', '--cactus_softlink', default = './runProgressiveCactus.sh', show_default=True, help='Name of softlinked Progressive Cactus batch script file.', type=click.Path(exists=False))
@click.option('-env', '--cactus_env_softlink', default = './cactus_environment', show_default=True, help='Name of softlinked Progressive Cactus virtual environment file.', type=click.Path(exists=False))
@click.option('-n', '--n_cpus', default = 16, show_default=True, help='Number of cpus used to convert hal 2 maf files.')
@click.option('-h2m', '--hal2maf_softlink', default = './hal2maf', show_default=True, help='Name of softlinked Progressive Cactus hal2maf program.', type=click.Path(exists=False))
@click.option('-h2m', '--nickname_file', default = '', show_default=True, help='File containing protID nickname in each line for all protIDs, can omit this file by leaving it blank.', type=click.Path(exists=False))
@click.option('-fi', '--fasta_path', default = './fasta_path/', show_default=True, help='Fasta path containing all of the input genomes. Genome naming must conform to xxx_[protID]_xxx.[fa/fasta].', type=click.Path(exists=False))
def run_cactus(fasta_output_path,cactus_run_directory,cactus_softlink, cactus_env_softlink, n_cpus, hal2maf_softlink, nickname_file = '', fasta_path = ''): #fixme get rid of '' and add to command line tool
    """Run multiple sequence alignment via Progressive Cactus on multiple species synteny blocks and export as maf files. Try to run softlink_cactus beforehand, else use official cactus paths instead of softlinks."""
    cactus_run_obj = CactusRun(fasta_output_path,cactus_run_directory,cactus_softlink, cactus_env_softlink, nickname_file, fasta_path)
    cactus_run_obj.write_fastas_seqfile()
    cactus_run_obj.run_cactus()
    cactus_run_obj.hal2maf(n_cpus, hal2maf_softlink)

@joshuatree.command()
@click.option('-cd', '--cactus_distribution_dir', default = './progressiveCactus/', show_default=True, help='Path containing installed Progressive Cactus Distribution.', type=click.Path(exists=False))
@click.option('-cac', '--softlink_cactus_name', default = './runProgressiveCactus.sh', show_default=True, help='Name of softlinked Progressive Cactus batch script file.', type=click.Path(exists=False))
@click.option('-env', '--softlink_env_name', default = './cactus_environment', show_default=True, help='Name of softlinked Progressive Cactus virtual environment file.', type=click.Path(exists=False))
@click.option('-h2m', '--softlink_hal2maf_name', default = './hal2maf', show_default=True, help='Name of softlinked Progressive Cactus hal2maf program.', type=click.Path(exists=False))
def softlink_cactus(cactus_distribution_dir,softlink_cactus_name, softlink_env_name, softlink_hal2maf_name):
    """Softlink cactus distribution's cactus bash script, virtual environment, and hal2maf program. Useful if installed Cactus to particular directory and want to save time in referencing that directory when running cactus."""
    subprocess.call('ln -s %s %s'%(os.path.abspath(cactus_distribution_dir+'/bin/runProgressiveCactus.sh'),softlink_cactus_name),shell=True)
    subprocess.call('ln -s %s %s'%(os.path.abspath(cactus_distribution_dir+'/environment'),softlink_env_name),shell=True)
    subprocess.call('ln -s %s %s'%(os.path.abspath(cactus_distribution_dir+'/submodules/hal/bin/hal2mafMP.py'),softlink_hal2maf_name),shell=True)

@joshuatree.command()
@click.option('-i','--install_path', help='Install Path.')
def install_cactus(install_path):
    """Install the Cactus distribution and hal tools.""" # fixme, make sure hal tools works
    os.chdir(install_path)
    conda_env = os.popen('echo $CONDA_PREFIX').read().split('/')[-1]
    subprocess.call('git config --global --add http.sslVersion tlsv1.2\ngit clone git://github.com/glennhickey/progressiveCactus.git\ncd progressiveCactus\ngit pull\ngit submodule update --init\nmake',shell=True)
    # source deactivate\n \nsource activate %s\n'%conda_env

# fixme, functions to add, nextflow progressiveCactus executor, move maf to new folder, visualize CS, breakdown into CNS elements, extract vcf, local PCA, local Trees, run iqtree, plottree, visualize tree with ete3, merge/intersect vcf, mafstrander, maffilter, vcf to snp matrix, plotPositions and all associated tools,
# fixme anything with vcf or SNPs should be one class; convert mafobject to vcf object, vcf2tab inside vcf object
# fixme main ideas extract vcf, CNS, visualize snps and cns, chrom spatial and PCA, compute phylogeny via SNPs, local tree topology, vcf operations,

####################
#### fixme add commands with snps and maf ####





#### RUN CLI ####

if __name__ == '__main__':
    joshuatree()

"""'PAC2_0.316-PAC4GC.655_5.unout'
>>> sorted([(m.start(0),m.end(0)) for m in re.finditer('\PAC2_0|PAC4GC',unout)])
[(0, 6), (11, 17)]"""