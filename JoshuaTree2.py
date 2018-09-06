#################
#### IMPORTS ####

import click, os, pandas as pd, numpy as np, subprocess, sys, shutil
from pybedtools import BedTool
from pyfaidx import Fasta
from itertools import combinations
import scipy.sparse as sps
import glob, re
from random import randint
from Bio import SeqIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
#import multiprocessing as mp
from time import sleep
from pathos import multiprocessing as mp
from pathos.pp_map import ppmap
from Bio import AlignIO
from ete3 import PhyloTree, Tree,TreeStyle,NodeStyle,EvolTree
import networkx as nx
import matplotlib.pyplot as plt
import pickle

#################
#### CLASSES ####


#########################
#### SYNTENY CLASSES ####

class SuperSynteny:
    def __init__(self, pairwise_syntenies, max_distance, q_genome = ''):
        self.pairwise_syntenies = {synteny.synteny_file: synteny for synteny in pairwise_syntenies}
        self.max_distance = max_distance
        self.q_genome = q_genome

    def generate_global_synteny_graph(self, fasta_out_dir, memory=150):
        synteny_graph = nx.Graph()
        for synteny in self.pairwise_syntenies.values():
            synteny_graph.add_edges_from(zip(synteny.q_genome.protID+'.'+synteny.synteny_structure['q_chr']+'\t'+synteny.synteny_structure['q_xi']+'\t'+synteny.synteny_structure['q_xf'],
                                             synteny.s_genome.protID+'.'+synteny.synteny_structure['s_chr']+'\t'+synteny.synteny_structure['s_xi']+'\t'+synteny.synteny_structure['s_xf']))
        genomes = {genome.protID: genome for genome in set(reduce(lambda x,y: x+y,[[synteny.q_genome, synteny.s_genome] for synteny in self.pairwise_syntenies.values()]))}
        nodes_all = np.array(synteny_graph.nodes())
        if self.q_genome:
            nodes = np.array(filter(lambda x: x.startswith(self.q_genome.protID), nodes_all))
        else:
            nodes = nodes_all
        clustered_regions = BedTool('\n'.join(nodes),from_string=True).sort().cluster(d=self.max_distance).to_dataframe().astype(str)
        nodes = nodes_all
        clustered_regions = clustered_regions.rename(dict(enumerate(clustered_regions['chrom']+'\t'+clustered_regions['start']+'\t'+clustered_regions['end'])))
        for cluster in clustered_regions['name'].unique():
            synteny_graph.add_edges_from(list(combinations(list(clustered_regions[clustered_regions['name'] == cluster].index),r=2)))
        synteny_graph.add_edges_from(zip(nodes,nodes))
        connected_components = enumerate(list(nx.connected_components(synteny_graph)))
        #species_colors = {species:i for i,species in enumerate(genomes.keys())}
        #node_labels = {}
        #for i, regions in connected_components:
        #    for region in regions:
        #        node_labels[region] = species_colors[region.split('.')[0]]#i
        #plt.figure()
        #nx.draw_spectral(synteny_graph,node_color=np.vectorize(lambda x: node_labels[x])(nodes),arrows=False,node_size=10)
        #plt.savefig('./'+'spectral_layout.png',dpi=300)
        for i, regions in connected_components:#list(nx.connected_component_subgraphs(synteny_graph)):
            regions = BedTool('\n'.join(regions),from_string=True).sort().merge(d=self.max_distance).to_dataframe().astype(dict(chrom=str,start=int,end=int))
            regions['species'] = regions['chrom'].map(lambda x: x.split('.')[0])
            regions['chrom'] = regions['chrom'].map(lambda x: x.split('.',1)[1])
            print(regions)

            with open(fasta_out_dir+'fasta_output.fa','w') as f:
                for species in regions['species'].unique():
                    f.write('\n'.join(['>%s.%s_%d_%d\n%s'%(species, chrom, start, end, str(genomes[species].fasta[chrom][start:end])) for chrom, start, end, species in map(tuple,regions[regions['species'] == species].values)]) + '\n')
            subprocess.call('reformat.sh overwrite=true in=%s out=%s fastawrap=60 Xmx=%dG && rm %s'%(fasta_out_dir+'fasta_output.fa',fasta_out_dir+'FastaOut%d.fasta'%i,memory,fasta_out_dir+'fasta_output.fa') , shell=True)
        #global_synteny_matrix = nx.to_scipy_sparse_matrix(synteny_graph,nodelist=nodes)
        #n_components, connected_components = sps.csgraph.connected_components(global_synteny_matrix,directed=False)
        #for connected_component in connected_components:
        #    regions = np.
        self.synteny_graph = synteny_graph
        self.genomes = genomes

    def visualize_synteny_graph(self, work_dir = './'):
        work_dir += '/'
        nodes = np.array(self.synteny_graph.nodes())
        connected_components = enumerate(list(nx.connected_components(self.synteny_graph)))
        species_colors = {species:i for i,species in enumerate(self.genomes.keys())}
        node_labels_species = {}
        node_labels_component = {}
        for i, regions in connected_components:
            for region in regions:
                node_labels_species[region] = species_colors[region.split('.')[0]]#i
                node_labels_component[region] = i
        plt.figure()
        nx.draw_spectral(self.synteny_graph,node_color=np.vectorize(lambda x: node_labels_species[x])(nodes),arrows=False,node_size=10)
        plt.savefig(work_dir+'spectral_layout_species.png',dpi=300)
        plt.figure()
        nx.draw_spectral(self.synteny_graph,node_color=np.vectorize(lambda x: node_labels_component[x])(nodes),arrows=False,node_size=10)
        plt.savefig(work_dir+'spectral_layout_component.png',dpi=300)


    def generate_global_synteny_graph_reference(self, fasta_out_dir): # fixme add global synteny graph with no references... all sequences if syntenic Aij = 1, if within merge distance Aij = 1, just merge and output all connected components, probably easier and more straight forward, add global option for this, maybe separate method and circos_dropper like synteny computations, networkx plot final nodes
        global_synteny_indices = []
        for synteny in self.pairwise_syntenies.values():
            synteny.synteny_structure['name'] = np.array(synteny.synteny_structure.index)#synteny.synteny_structure['q_chr'] + '\t' + synteny.synteny_structure['q_xi'] + '\t' + synteny.synteny_structure['q_xf']#[['q_chr','q_xi','q_xf']]
            global_synteny_indices.extend([(synteny.synteny_file,name) for name in synteny.synteny_structure['name'].as_matrix().tolist()])
        global_synteny_indices = pd.MultiIndex.from_tuples(global_synteny_indices,names=['synteny_file','region'])
        global_synteny_matrix = pd.DataFrame(np.zeros((len(global_synteny_indices),)*2),index=global_synteny_indices,columns=global_synteny_indices)
        global_synteny_matrix = global_synteny_matrix.sort_index(axis=0).sort_index(axis=1)
        referenced_regions = list(global_synteny_matrix.index)
        clustered_regions = BedTool(np.array(referenced_regions)[:,1].tolist()).sort().cluster(d=self.max_distance).to_dataframe().astype(str)
        #print(clustered_regions)
        clustered_regions = clustered_regions.rename(dict(enumerate(clustered_regions['chrom']+'\t'+clustered_regions['start']+'\t'+clustered_regions['end'])))

        for cluster in clustered_regions['name'].unique():
            list_regions = list(clustered_regions[clustered_regions['name'] == cluster].index)
            for region_i,region_j in list(combinations(list_regions,r=2)):
                global_synteny_matrix.loc[pd.IndexSlice[:,region_i],pd.IndexSlice[:,region_j]] = 1
                global_synteny_matrix.loc[pd.IndexSlice[:,region_j],pd.IndexSlice[:,region_i]] = 1
            for region_i in list_regions:
                global_synteny_matrix.loc[pd.IndexSlice[:,region_i],pd.IndexSlice[:,region_i]] = 1 # fixme add global mode

        """
        regions_bed_files = dict(zip(regions,map(lambda x: BedTool(x,from_string=True),regions)))
        for region_i,region_j in list(combinations(regions,r=2)) + zip(regions,regions):
            if region_i == region_j:
                    global_synteny_matrix.loc[pd.IndexSlice[:,region_i],pd.IndexSlice[:,region_i]] = 1
            else:
                distance = int(str(regions_bed_files[region_i].closest(regions_bed_files[region_j],d=True)).split('\t')[-1])
                if distance >=0 and distance <= self.max_distance:
                    global_synteny_matrix.loc[pd.IndexSlice[:,region_i],pd.IndexSlice[:,region_j]] = 1
                    global_synteny_matrix.loc[pd.IndexSlice[:,region_j],pd.IndexSlice[:,region_i]] = 1"""
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
            with open(fasta_out_dir+'/'+'FastaOut%d.fa'%i,'w') as f:
                f.write(reduce(lambda x,y: x+'\n'+y,out_fasta))

    def test(self,region1, region2):
        syntenic_sequences = ['future_testing_in_development_global_mode']
        region1 = region1.split()
        region2 = region2.split()
        return 1 if region1 == region2 or (region1[0] == region2[0] and (int(region1[2]) + self.max_distance >= int(region2[1]) or int(region2[2]) + self.max_distance >= int(region1[1]))) or (region1,region2) in syntenic_sequences else 0

class PairwiseSynteny:
    def __init__(self, q_genome, s_genome, synteny_file='', loci_threshold = 4):
        self.synteny_file = synteny_file
        self.q_genome = q_genome
        self.s_genome = s_genome
        self.loci_threshold = loci_threshold
        self.chrom_colors = {}

    def generate_synteny_structure(self,synteny_path):
        """Take anchor file or synteny file and searches for starting and ending genes for each syntenic block"""
        if self.synteny_file.endswith('.unout'):
            self.unout2structure(self.q_genome, self.s_genome)
        elif self.synteny_file.endswith('.lifted.anchors'):
            self.anchor2structure(self.q_genome, self.s_genome)
        elif self.synteny_file.endswith('.bed'):
            self.import_synteny_structure()
        else:
            self.run_synteny(self.q_genome,self.s_genome,synteny_path)
            self.anchor2structure(self.q_genome, self.s_genome)
        self.synteny_structure_index()

    def synteny_structure_index(self):
        self.synteny_structure = self.synteny_structure.rename(dict(zip(range(self.synteny_structure.shape[0]),self.synteny_structure['q_chr']+'\t'+self.synteny_structure['q_xi']+'\t'+self.synteny_structure['q_xf'])))

    def import_synteny_structure(self):
        self.synteny_structure = pd.read_table(self.synteny_file,header=None,names=['q_chr','q_xi','q_xf','s_chr','s_xi','s_xf'])
        self.synteny_structure_index()
        #self.synteny_structure = self.synteny_structure.rename(dict(enumerate((self.synteny_structure['q_chr']+'\t'+self.synteny_structure['q_xi']+'\t'+self.synteny_structure['q_xf']).as_matrix().tolist())))

    def unout2structure(self,q_genome,s_genome):
        with open(self.synteny_file,'r') as f:
            lines = np.array(f.read().splitlines())
        anchors = np.array_split(lines,np.where(np.vectorize(lambda line: line.startswith('\t') == 0)(lines))[0])
        synteny_structure = []
        for anchor in anchors:
            anchor = np.array(map(lambda line: line.split('\t')[1].split(','),anchor[1:].tolist()))
            if len(anchor) >= self.loci_threshold and anchor.tolist():
                q_genes, s_genes = anchor[:,2], anchor[:,5]
                q_coords, s_coords = q_genome.df.loc[q_genes,:], s_genome.df.loc[s_genes,:]
                synteny_structure.append([q_coords.iloc[0,0],q_coords[['xi','xf']].values.min(),q_coords[['xi','xf']].values.max(),s_coords.iloc[0,0],s_coords[['xi','xf']].values.min(),s_coords[['xi','xf']].values.max()])
        self.synteny_structure = pd.DataFrame(synteny_structure,columns=['q_chr','q_xi','q_xf','s_chr','s_xi','s_xf']).astype(str)#,index = np.vectorize(lambda x: '\t'.join(map(str,x[:3])))(synteny_structure))

    def run_synteny(self,genome1,genome2, synteny_path):
        pwd = os.getcwd()
        os.chdir(synteny_path)
        #subprocess.call('rm {0}/*.bck {0}/*.prj {0}/*.sds {0}/*.ssp {0}/*.suf {0}/*.tis {0}/*.des {0}/*.bed {0}/*.cds'.format(synteny_path),shell=True)
        for abs_path, link_name in zip([genome1.bed_file,genome2.bed_file,genome1.CDS_file,genome2.CDS_file],[genome1.protID+'.bed',genome2.protID+'.bed',genome1.protID+'.cds',genome2.protID+'.cds']):
            subprocess.call('ln -s %s %s'%(abs_path,link_name),shell=True)
        try:
            subprocess.call('python -m jcvi.compara.catalog ortholog --no_strip_names %s %s'%(genome1.short_name,genome2.short_name),shell=True)
        except:
            subprocess.call('python -m jcvi.compara.catalog ortholog --no_strip_names %s %s'%(genome1.short_name,genome2.short_name),shell=True)
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
            if anchor:
                genes = np.array([line.split()[:2] for line in anchor.splitlines() if line])
                if genes.shape[0] >= self.loci_threshold:
                    q_genes, s_genes = genes[:,0] , genes[:,1]
                    q_coords, s_coords = q_genome.df.loc[q_genes,:], s_genome.df.loc[s_genes,:]
                    #print q_coords[['xi','xf']]
                    synteny_structure.append([q_coords.iloc[0,0],q_coords[['xi','xf']].values.min(),q_coords[['xi','xf']].values.max(),s_coords.iloc[0,0],s_coords[['xi','xf']].values.min(),s_coords[['xi','xf']].values.max()])
        self.synteny_structure = pd.DataFrame(synteny_structure,columns=['q_chr','q_xi','q_xf','s_chr','s_xi','s_xf']).astype(str)

    def synteny_structure_2_bed(self,filename):
        df = self.synteny_structure
        df['q_chr'] = np.vectorize(lambda x: self.q_genome.protID+'-'+x)(df['q_chr'])
        df['s_chr'] = np.vectorize(lambda x: self.s_genome.protID+'-'+x)(df['s_chr'])
        df.to_csv(filename,sep='\t',index=False,header=None)

    def synteny_structure_2_link(self, filename, bundle_links = False, link_gap = 10000):
        self.chrom_colors = dict(self.q_genome.chrom_colors,**self.s_genome.chrom_colors)
        click.echo(str(self.chrom_colors))
        df = self.synteny_structure
        df['q_chr'] = np.vectorize(lambda x: self.q_genome.protID+'-'+x)(df['q_chr'])
        df['s_chr'] = np.vectorize(lambda x: self.s_genome.protID+'-'+x)(df['s_chr'])
        if 0:
            df = df[['s_chr','s_xi','s_xf','q_chr','q_xi','q_xf']]
        df.to_csv(filename,sep=' ',index=False,header=None)
        if bundle_links:
            click.echo("./helper_scripts/circos-tools-0.22/tools/bundlelinks/bin/bundlelinks -max_gap {0} -links {1} > {1}.temp && cut -f 1-6 -d " " {1}.temp > {1} && rm {1}.temp".format(str(link_gap), filename))
            subprocess.call("./helper_scripts/circos-tools-0.22/tools/bundlelinks/bin/bundlelinks -max_gap {0} -links {1} > {1}.temp && cut -f 1-6 -d \" \" {1}.temp > {1} && rm {1}.temp".format(str(link_gap), filename), shell=True)
        if 0:
            df = pd.read_table(filename ,sep = ' ', header = None, names = ['q_chr','q_xi','q_xf','s_chr','s_xi','s_xf'])
            df['color'] = np.vectorize(lambda x: 'color=%s'%(self.chrom_colors[x]))(df['q_chr'])
            print(df)
            df.to_csv(filename, sep=' ', index=False, header=None)
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
        self.bed_file = bed_file#os.path.abspath(bed_file)
        self.short_name = self.bed_file.split('/')[-1].replace('.bed3','').replace('.bed','')
        self.protID = protID
        self.gff_file = gff_file
        self.chrom_colors = {}
        if self.gff_file and (os.path.exists(self.bed_file) == 0 or (os.path.exists(self.bed_file) and os.stat(self.bed_file).st_size == 0)):
            #click.echo('python -m jcvi.formats.gff bed --type=mRNA --key=%s %s > %s'%(self.gene_info,self.gff_file,self.bed_file))
            #print('python -m jcvi.formats.gff bed %s --type=mRNA --key=%s -o %s'%(self.gff_file,self.gene_info,self.bed_file))
            subprocess.call('python -m jcvi.formats.gff bed %s --type=mRNA --key=%s -o %s'%(self.gff_file,self.gene_info,self.bed_file),shell=True)#FIXME
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
        self.bed_file = os.path.abspath(self.bed_file)
        self.CDS_file = self.bed_file.replace('.bed3','.cds').replace('.bed','.cds')
        self.df = pd.read_table(self.bed_file,header=None,names=['chr','xi','xf','Gene'],dtype={'chr':str,'xi':np.int,'xf':np.int,'Gene':str},usecols=[0,1,2,3])
        self.df = self.df.set_index('Gene')

    def export_bed(self,filename):
        df = self.df.reset_index().rename(dict(index='Gene'),axis='columns').reindex(columns=['chr','xi','xf','Gene'])
        df.to_csv(filename,sep='\t',index=False,header=None)

    def extract_CDS(self): # python -m jcvi.formats.gff uniq t.PAC2_0.316.gff3 -o uniq.gff3
        if os.path.exists(self.CDS_file) == 0 or (os.path.exists(self.CDS_file) and os.stat(self.CDS_file).st_size == 0):
            subprocess.call('python -m jcvi.formats.gff load %s %s --parents=mRNA --children=CDS --id_attribute=%s -o %s'%(self.gff_file,self.fasta_file,self.gene_info,self.CDS_file),shell=True)

    def export_karyotype(self, filename, n_chromosomes=25, shorten_chr=False, chrom_file = ''):
        df = pd.read_table(self.fasta_file+'.fai', header=None,names=['chr','length'],usecols=[0,1],dtype=dict(zip(['chr','length'],[str,np.int])))
        df = df.sort_values(['length'],ascending=False)
        chromosomes_names_dict = {}
        if chrom_file:
            with open(chrom_file) as f:
                chromosomes = f.read().splitlines()
            if len(chromosomes[0].split()) == 2:
                chr_tuples = map(lambda l: tuple(l.split()), chromosomes)
                chromosomes_names_dict = dict(chr_tuples)
                chromosomes = np.array(chr_tuples)[:,0]
            else:
                chromosomes = np.array(chromosomes)
            #print(df[np.isin(df['chr'].values,chromosomes)].set_index('chr'))

            df = df[np.isin(df['chr'].values,chromosomes)].set_index('chr').reindex(chromosomes).reset_index().rename({'index':'chr'})[['chr','length']]
        if n_chromosomes < df.shape[0]:
            df = df.iloc[:n_chromosomes,:].reset_index(drop=True)
        out_txt = []
        for i in range(df.shape[0]):
            chrom = df.loc[i,'chr']
            chr_name = (chromosomes_names_dict[chrom] if chromosomes_names_dict else chrom) if not shorten_chr else chrom[0] + chrom.split('_')[-1]
            if i >= 25:
                color = '%d,%d,%d'%(randint(1,255),randint(1,255),randint(1,255))
            else:
                color = 'chr%d'%(i+1)
            out_txt.append('chr - %s-%s %s 0 %d %s\n'%(self.protID,chrom,chr_name,df.loc[i,'length'],color))
            self.chrom_colors['%s-%s'%(self.protID,chrom)] = color
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

    def generate_config(self, ticks = 'txticks.conf', ideogram = 'txideogram.conf', links_and_rules = 'linksAndrules.conf', config='circos.conf', variable_thickness = False, thickness_factor=1000, switch_lines = False):
        colors = pd.read_table(self.synteny.s_genome.karyotype if not switch_lines else self.synteny.q_genome.karyotype,header=None,usecols=[2,6],sep=' ').as_matrix()
        self.links_and_rules = links_and_rules
        self.config = config
        self.ideogram,self.ticks = ideogram, ticks
        if hasattr(self, 'ticks'):
            self.write_ticks_config(self.ticks)
        if hasattr(self, 'ideogram'):
            self.write_ideogram_config(self.ideogram)
        with open(self.config,'w') as f:
            f.write("""# circos.conf
                karyotype = %s, %s
                chromosomes_units = 1000000
                chromosomes_display_default = yes
                <<include %s>>
                <<include %s>>
                <<include %s>>
                <image>
                <<include etc/image.conf>>
                </image>
                <<include etc/colors_fonts_patterns.conf>>
                <<include etc/housekeeping.conf>>
                """%(self.synteny.q_genome.karyotype , self.synteny.s_genome.karyotype, self.ideogram,self.ticks,self.links_and_rules))
        with open(self.links_and_rules,'w') as f:
            f.write("""
                <links>
                <link>
                file = %s
                radius = 0.99r
                bezier_radius = 0r
                %s
                ribbon = yes
                color = black
                <rules>
                <rule>
                condition = var(intrachr)
                show = no
                </rule>\n"""%(self.synteny.link, 'thickness = eval(max(1,round(var(size1)/%d)))'%thickness_factor if variable_thickness else '') + '\n'.join(['<rule>\ncondition = %s(%s)\ncolor = %s\n</rule>'%('from' if switch_lines else 'to',chrom,color) for chrom,color in (self.synteny.q_genome.chrom_colors.items() if switch_lines else self.synteny.s_genome.chrom_colors.items())]) + '\n</rules>\n</link>\n</links>') # map(tuple,colors)

    def run_circos(self, output_dir='./', pdf=False):
        subprocess.call('circos -conf %s -outputfile %s-%s -outputdir %s'%(self.config,self.synteny.q_genome.protID,self.synteny.s_genome.protID,output_dir),shell=True)
        if pdf:
            subprocess.call('convert %s/%s-%s.png %s/%s-%s.pdf'%(os.path.abspath(output_dir),self.synteny.q_genome.protID,self.synteny.s_genome.protID,os.path.abspath(output_dir),self.synteny.q_genome.protID,self.synteny.s_genome.protID),shell=True)

##########################
#### CACTUS RUN CLASS ####

class CactusRun:
    def __init__(self,fasta_output_path,cactus_run_directory,cactus_softlink, nickname_file = '', fasta_path = ''):
        self.fasta_output_path = fasta_output_path +'/'
        self.cactus_run_directory = cactus_run_directory +'/'
        self.cactus_output = self.cactus_run_directory+'output/'
        self.hal_path = self.cactus_output+'hal/'
        self.cactus_softlink = os.path.abspath(cactus_softlink+'/bin/runProgressiveCactus.sh')
        if not nickname_file:
            self.nickname_file = self.cactus_run_directory + 'prot_dict'
            self.protIDs = [fasta.split('_')[-2] for fasta in glob.glob(fasta_path+'/*.fa')+glob.glob(fasta_path+'/*.fasta')]
            with open(self.nickname_file,'w') as f:
                f.write('\n'.join(['\t'.join((protID,)*2) for protID in self.protIDs]))
        else:
            self.nickname_file = nickname_file
        self.cactus_env_softlink = os.path.abspath(cactus_softlink+'/environment')
        #subprocess.call('source %s'%cactus_env_softlink,shell=True)

    def fasta2seq(self, fasta):
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
            f.write('#!/bin/bash\nexport _JAVA_OPTIONS="-Xmx155g"\n%s --maxThreads 16 %s %s %s >& %s\nscp %s %s'%(os.path.abspath(self.cactus_softlink),self.seqfile,self.fasta_run_dir,fasta.replace('.fasta','.hal'),self.seqfile+'.sh.stdout',fasta.replace('.fasta','.hal'),self.hal_path+fasta.split('/')[-1].replace('.fasta','.hal')))
        return run_file

    def write_fastas_seqfile(self):
        with open(self.nickname_file,'r') as f:
            self.nicknames = dict([tuple(line.split()) for line in f.read().splitlines() if line])
        fastas = glob.glob(self.fasta_output_path+'/FastaOut*.fa')+glob.glob(self.fasta_output_path+'/FastaOut*.fasta')
        self.run_files = []
        fasta2seq = lambda fasta: self.fasta2seq(fasta)
        p = mp.ProcessingPool()
        r = p.amap(fasta2seq,fastas)
        r.wait()
        self.run_files = r.get()
        p.close()

    def run_cactus(self, submission = 'local', shifter = False):
        with open('nextflow.config','w') as f:
            f.write('\n'.join(["process.%s = '%s'"%(i,j) for i,j in zip(['executor','memory', 'clusterOptions'],[submission,'155G', '' if submission != 'sge' else '-P plant-analysis.p -cwd -l h_rt=24:00:00 -pe pe_slots 16 -e OutputFile.txt'])]+(["shifter.enabled = true"] if shifter else ["docker.enabled = true"])+["process.container = 'lepbase/progressive-cactus:latest'"]))
        subprocess.call("export SHIFTER_RUNTIME='' && nextflow cactus_run.nf --work_dir %s --cactus_run_files %s --environment %s"%(os.getcwd(),','.join(os.path.abspath(run_file) for run_file in self.run_files),os.path.abspath(self.cactus_env_softlink)),shell=True)
        #for run_file in self.run_files:
        #    subprocess.call('nohup sh %s &'%(run_file),shell=True)
    # fixme, make sure to activate progressive cactus environment; sourcec environment in cactus folder, add hal2maf, send maf 2 cns analysis, parallelize and pipeline all scripts, maybe call nextflow from within python for each of the jobs for pipeline submission

    def generate_trees(self,scaled = 10000, kmer_length = 23, multi_fasta = False):
        # fixme use mafft to generate new guide trees for each set of input fasta files
        subprocess.call('sourmash compute -f --scaled %d %s/*.fa -o %s.sig -k %d %s'%(scaled,self.fasta_run_dir,self.fasta_run_dir+'tree',kmer_length,'--singleton' if multi_fasta else ''),shell=True)
        subprocess.call('sourmash compare %s.sig --csv %s.cmp.csv'%(self.fasta_run_dir+'tree',self.fasta_run_dir+'tree'),shell=True)
        df = pd.read_csv('%s.cmp.csv'%(self.fasta_run_dir+'tree'),index_col = None)
        samples = [fasta.split('/')[-1].replace('.fa','') for fasta in list(df)]
        distance_matrix = df.as_matrix()
        constructor = DistanceTreeConstructor()
        dm = _DistanceMatrix(names=samples,matrix=[list(distance_matrix[i,0:i+1]) for i in range(len(samples))])
        tree = constructor.nj(dm)
        Phylo.write(tree,self.fasta_run_dir+'output_tree.nh','newick') # fixme bug, tree being generated has negative branch lengths, this is why cactus is failing

    def hal2maf(self, n_cpus, hal2maf_softlink):
        self.hal2maf_softlink = hal2maf_softlink
        self.maf_path = self.hal_path.replace('hal','maf')
        run_hal = lambda hal: self.run_hal2maf(hal)
        p = mp.ProcessingPool()
        r = p.amap(run_hal, glob.glob(self.hal_path+'/*.hal'))
        r.wait()
        p.close()
        """
        for hal in glob.glob(self.hal_path+'/*.hal'):
            proc = mp.Process(target=run_hal, args=(hal,))
            proc.daemon = True
            proc.start()
            while len(mp.active_childern()) > n_cpus:
                sleep(1)
        while len(mp.active_childern()) > 0:
            sleep(1)"""

    def run_hal2maf(self,hal):
        os.system('%s %s %s'%(self.hal2maf_softlink, hal, hal.replace('hal','maf')))

###################
#### MAF CLASS #### fixme try to inheret from circos class, maf class will allow you to extract CNS or VCF, visualize CS and CNS

class MAF_filter_config:
    def __init__(self, config_file='maf_filter_config.bpp', input_file='merged.maf',species = 'list_species.txt', reference_species = '', log_file = 'log.out', out_all_species = True):
        self.config_file = config_file
        self.config_txt = ''
        self.input_file = input_file
        self.input_format = 'Maf'
        self.log_file = log_file
        self.all_species = species.split(',') if not species.endswith('.txt') else open(species,'r').read().splitlines()
        all_species_but_one = set(self.all_species) - {reference_species}
        if out_all_species:
            self.species = self.all_species
        else:
            self.species = all_species_but_one
        self.reference_species = reference_species
        self.endtxt = []

    def export(self):
        with open(self.config_file,'w') as f:
            f.write(self.config_txt)

    def add_subset_text(self, species = '', keep=True):
        self.endtxt.append("""Subset(\\
                strict=yes,\\
                keep=%s,\\
                species=(%s),\\
                remove_duplicates=yes),\\"""%('yes' if keep else 'no', species if species else ','.join(self.all_species)))


    def add_vcf_text(self, vcf_file='vcfs/merged.vcf'):
        self.endtxt.append("""VcfOutput(\\
                file=%s,\\
                genotypes=(%s),\\
                all=no,\\
                reference=%s),\\"""%(vcf_file,','.join(self.species),self.reference_species))

    def create_txt(self):
        self.config_txt += """input.file=%s
    input.format=%s
    output.log=%s
    maf.filter=\\\n"""%(self.input_file,self.input_format,self.log_file) +'\n'.join(self.endtxt[:-3])
        with open(self.config_file,'w') as f:
            f.write(self.config_txt)

    def run_maffilter(self):
        subprocess.call('./helper_scripts/maffilter param=%s'%self.config_file,shell=True)


class MAF:
    def __init__(self, maf_file, special_format = True):
        self.maf_file = maf_file
        self.special_format = special_format

    def merge(self, new_maf_file='merged.maf'):
        if '*' in new_maf_file:
            maf_files = glob.glob(self.maf_file)
        else:
            maf_files = self.maf_file.split(',')
        subprocess.call('rm %s'%new_maf_file,shell=True)
        subprocess.call("(echo '##maf version=1'; ( cat %s | sed -e '/Anc/d;/#/d' ) ;) > %s"%(' '.join(maf_files),new_maf_file),shell=True)
        self.maf_file = new_maf_file

    def index(self):
        idxs = []
        with open(self.maf_file,'r') as f:
            offset = 0
            for line in iter(f.readline, ''):
                if line.startswith('a'):
                    idxs.append(offset)
                offset = f.tell()
            idxs.append(f.tell())
        idxs = sorted(set(idxs))
        self.idx = {idxs[i]:idxs[i+1] for i in range(len(idxs)-1)}
        self.idxs = sorted(self.idx.keys())

    def change_coordinates(self, reference_species,reference_species_chromosomes, changed_coordinates_file):
        def change_segment(segment, ref_species, ref_species_chr):
            aln_lines = segment.splitlines()
            for i,line in enumerate(aln_lines):
                if line.startswith('s'):
                    lineList = line.split()
                    orientation = lineList[4]
                    lineList2 = lineList[1].split('.')
                    lineList3 = lineList2[-1].split('_')[-2:]
                    lineList2[-1] = lineList2[-1].replace('_'+'_'.join(lineList3),'')
                    lineList[1] = '.'.join(lineList2[1:])
                    if lineList2[0] == ref_species:
                        chrom = lineList2[-1]
                        lineList[5] = str(ref_species_chr[chrom])
                        if orientation == '-':
                            lineList[2] = str(ref_species_chr[chrom]-int(lineList3[-1])+int(lineList[2]))#-int(lineList[3]))
                        else:
                            lineList[2] = str(int(lineList3[-2]) + int(lineList[2]))
                        position = int(lineList[2])
                    else:
                        lineList[2] = str(int(lineList3[-2]) + int(lineList[2]))
                    aln_lines[i] = '\t'.join(lineList)
            try:
                return chrom,position,'\n'.join(sorted(filter(None,aln_lines)))+'\n\n'
            except:
                return '', '', ''

        chunks = [self.idxs[i:i+50000] for i in range(0,len(self.idxs),50000)]
        with open(self.maf_file,'r') as f, open(changed_coordinates_file,'w') as f2:
            for chunk in chunks:
                out_segments = []
                for idx in chunk:
                    f.seek(idx)
                    chrom, position, segment = change_segment(f.read(self.idx[idx] - idx),reference_species,reference_species_chromosomes)
                    if chrom:
                        out_segments.append(segment)
                f2.write(''.join(out_segments))
        self.maf_file_new_coords = changed_coordinates_file

    def strand(self, reference_species):
        subprocess.call('./helper_scripts/mafStrander -m %s --seq %s > temp.maf'%(self.maf_file,reference_species),shell=True)
        subprocess.call('mv temp.maf %s'%self.maf_file)

    def maf2vcf(self, maf_filter_config, species, reference_species, reference_species_fai, vcf_out, change_coordinates = True):
        """Run on a merged maf file first by using merger."""
        reference_species_chromosomes = dict(zip(os.popen("awk '{print $1}' %s"%reference_species_fai).read().splitlines(),map(int,os.popen("awk '{print $2}' %s"%reference_species_fai).read().splitlines())))
        try:
            os.mkdir('vcfs')
        except:
            pass

        self.strand(reference_species)

        self.index()

        if change_coordinates:
            self.change_coordinates(reference_species,reference_species_chromosomes, self.maf_file.replace('.maf','.new_coords.maf'))
        else:
            self.maf_file_new_coords = self.maf_file

        self.config = MAF_filter_config(maf_filter_config, self.maf_file_new_coords, species, reference_species)
        self.config.add_subset_text(keep=False)
        self.config.add_subset_text(species=reference_species)
        self.config.add_vcf_text(vcf_out)
        self.config.run_maffilter()

        vcf_obj = SNP(vcf_out)
        vcf_obj.concat_vcf(vcf_out)
        vcf_obj.generate_new_header(vcf_out)


###################
#### SNP CLASS #### fixme vcf or tab file, as well as df, can local pca, local tree, final tree (run iqtree), visualize tree, produce enw vcf files and interact with other snp objects
# fixme maybe add nextflow class?
# fixme add test classes

class SNP:
    def __init__(self, snp_file, snp_format = 'vcf'):
        self.snp_file = snp_file
        self.format = snp_format
        if snp_format == 'tab':
            self.tab = self.snp_file
            print("Warning: Do not use vcf manipulation methods.")

    def concat_vcf(self, vcf_out):
        if self.format == 'vcf':
            list_vcfs = self.snp_file.split(',')
            master_df = pd.DataFrame()
            header_lines = []
            print list_vcfs
            for vcf_in in list_vcfs:
                with os.popen(('z' if vcf_in.endswith('.gz') else '')+'cat %s'%vcf_in) as f:
                    line_count = 0
                    for line in f:
                        header_lines.append(line)
                        if line.startswith('#CHROM'):
                            line_info = line.strip('/n').split() # FIXME can grab header line number here
                            break
                        line_count += 1
                if vcf_in.endswith('.gz'):
                    master_df = master_df.append(pd.DataFrame(np.hstack([np.array(os.popen(('z' if vcf_in.endswith('.gz') else '')+"cat %s | grep -v ^# | awk '{ print $%d }'"%(vcf_in,i+1)).read().splitlines())[:,None] for i in range(len(line_info))]),columns = line_info))
                else:
                    master_df = master_df.append(pd.read_table(vcf_in,header=line_count))
            header_lines = set(header_lines)
            master_df['POS'] = np.vectorize(int)(master_df['POS'])
            master_df = master_df.sort_values(['#CHROM','POS'])
            master_df.to_csv(vcf_out,sep='\t',index=False, na_rep = '.')
            with open(vcf_out.replace('.vcf','.headers.vcf'),'w') as f, open(vcf_out,'r') as f2:
                for line in [line2 for line2 in header_lines if '#CHROM' not in line2]:
                    f.write(line)
                f.write(f2.read())
            subprocess.call('mv %s %s'%(vcf_out.replace('.vcf','.headers.vcf'),vcf_out),shell=True)
            self.snp_file = vcf_out
        else:
            print('File(s) must be in vcf/vcf.gz format before proceeding.')

    def generate_new_header(self, vcf_out):
        if self.format == 'vcf':
            sort_vcf_in = self.snp_file
            header_line = '\n'+os.popen("grep '^#CHROM' %s"%sort_vcf_in).read().strip('\n')
            chrms = set(os.popen("awk '{ print $1}' %s | grep -v ^#"%sort_vcf_in).read().splitlines())
            new_lines = """##fileformat=VCFv4.1\n"""+'\n'.join(sorted(['##contig=<ID=' + chrm + ',length=' + os.popen('grep %s %s | tail -n 1'%(chrm,sort_vcf_in)).read().strip('/n').split()[1] + '>' for chrm in chrms]))+'\n'+'\n'.join(['##FILTER=<ID=gap,Description="At least one sequence contains a gap">','##FILTER=<ID=unk,Description="At least one sequence contains an unresolved character">','##FILTER=<ID=gunk,Description="At least one sequence contains an unresolved character and gap.">','##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'])+header_line
            with open('new_header.txt','w') as f:
                f.write(new_lines+'\n')
            subprocess.call("""(cat new_header.txt;
                                sed 's/gap,unk/gunk/g' %s | grep -v ^#;) \
                                > %s"""%(sort_vcf_in,sort_vcf_in.replace('.vcf','.new_head.vcf')),shell=True)
            df = pd.read_table(sort_vcf_in.replace('.vcf','.new_head.vcf'),header=new_lines.count('\n')).fillna('.')
            df['INFO'] = '.'
            df['FORMAT'] = 'GT'
            df.to_csv(sort_vcf_in.replace('.vcf','.new_head.vcf'),sep='\t',index=False, na_rep = '.')
            subprocess.call("""(cat new_header.txt;
                                cat %s | grep -v ^#;) \
                                > %s"""%(sort_vcf_in.replace('.vcf','.new_head.vcf'),vcf_out),shell=True)# sort_vcf_in.replace('.vcf','.new_head_final.vcf')
            self.snp_file = vcf_out
        else:
            print('File must be in vcf format before proceeding.')

    def merge_vcf(self, vcf_out, excluded_species, reference_species, *list_vcf_obj):
        #list_vcfs = list_vcfs.split(',')
        list_vcfs = [vcf_obj.snp_file for vcf_obj in list_vcf_obj]
        excluded_species = ','.join(excluded_species.split(',')+['%d:%s'%(i+2,reference_species) for i in range(len(list_vcfs) - 1)])
        list_vcfs2 = []
        for vcf in list_vcfs:
            subprocess.call("(cat %s | grep ^# ; bedtools intersect -a %s -b %s ;) > %s"%(vcf,vcf,' '.join([v for v in list_vcfs if v != vcf]),vcf.replace('.vcf','.intersect.vcf')),shell=True)
            vcf = vcf.replace('.vcf','.intersect.vcf')
            subprocess.call('bcftools view %s -O z -o %s'%(vcf, vcf+'.gz'),shell=True)
            subprocess.call('bcftools index %s.gz'%vcf,shell=True)
            list_vcfs2.append(vcf)
        subprocess.call('bcftools merge %s -O v --force-samples -o %s.gz'%(' '.join([vcf+'.gz' for vcf in list_vcfs2]), vcf_out),shell=True)
        subprocess.call('bcftools view %s -O v -s ^%s -o %s'%(vcf_out+'.gz',excluded_species,vcf_out), shell=True)
        return SNP(vcf_out)

    def vcf2tab(self,tab_out):
        subprocess.call("bcftools query -Hf '%%CHROM\\t%%POS\\t%%REF[\\t%%TGT]\\n' -o %s %s"%(tab_out,self.snp_file),shell=True)
        self.tab = tab_out

    def tab2fasta(self, fasta_out, sample = 0, tab_in = '', df = ''):
        if tab_in:
            self.tab = tab_in
        if not list(df):
            sample = sample
            df = pd.read_table(self.tab,header=0)
            if sample:
                df = df.sample(n=sample)
        with open(fasta_out,'w') as f:
            for col in list(df)[3:]:
                species = col[col.find(']')+1:col.rfind(':')]
                f.write('>%s\n%s\n'%(species,''.join((df[col].as_matrix())).replace('.','-').replace('*','-')))
        return fasta_aln(fasta_out)

    def intersect_vcf(self, bed_regions, vcf_out):
        subprocess.call('bedtools intersect -wa -a %s -b %s > %s'%(self.snp_file, bed_regions, vcf_out),shell=True)
        return SNP(vcf_out)

    def visualize_distribution(self, work_dir):
        import seaborn as sns, matplotlib.pyplot as plt
        work_dir += '/'
        df = pd.read_table(self.snp_file,header=next(i for i,line in enumerate(open(self.snp_file,'r')) if line.startswith('#CHROM')))
        for chrom in df['#CHROM'].unique():
            dff = df[df['#CHROM'] == chrom]
            plt.figure()
            sns.distplot(dff['POS'].as_matrix())
            plt.title(chrom + ' SNP distribution')
            plt.savefig(work_dir+chrom+'_snp_distribution.png',dpi=300)

    def erase_indels(self, vcf_out):
        subprocess.call("vcftools --vcf %s --remove-indels --recode --recode-INFO-all --out %s"%(self.snp_file,vcf_out),shell=True)
        return SNP(vcf_out)

    def run_local_trees(self,snps_interval, phylogeny, work_dir):
        vcf_file = os.path.abspath(self.snp_file)
        tab_file = vcf_file.replace('.vcf','.tab')
        work_dir = os.path.abspath(work_dir)+'/'
        subprocess.call("export OPENBLAS_NUM_THREADS=1 && nextflow run local_trees.nf --vcf_file %s --tab_file %s --snps_interval %d --phylogeny %s --work_dir %s"%(vcf_file,tab_file,snps_interval,phylogeny,work_dir),shell=True)

    def tab2chunks(self, tab_in, snp_intervals, write_file):
        tab_df = pd.read_table(tab_in,header=0)
        tab_df['[2]POS'] = tab_df['[2]POS'].as_matrix().astype(np.int)
        tab_df = tab_df.sort_values(['# [1]CHROM','[2]POS'])
        subprocess.call('rm %s'%write_file,shell=True)
        for index, df in tab_df.groupby(np.arange(len(tab_df))//snp_intervals):
            chroms = set(df['# [1]CHROM'])
            if len(chroms) == 1:
                interval_data = '_'.join([list(chroms)[0],str(np.min(df['[2]POS'])),str(np.max(df['[2]POS']))]) # fixme maybe try to do in parallel
                self.tab2fasta('./%s.fa'%interval_data,0,'', df = df)
                with open(write_file,'a') as f:
                    f.write('\t'.join([interval_data,os.path.abspath('./%s.fa'%interval_data)])+'\n')

class fasta_aln:
    def __init__(self, fasta):
        self.fasta = fasta

    def to_snp(self, vcf_out ,fasta_out ,monomorphic):
        subprocess.call('snp-sites -v%s -o %s %s'%('b' if int(monomorphic) else '',vcf_out,self.fasta),shell=True)
        subprocess.call('snp-sites -m%s -o %s %s'%('b' if int(monomorphic) else '',fasta_out,self.fasta),shell=True)
        return SNP(vcf_out), fasta_aln(fasta_out)

    def generate_tree(self, phylogeny, tree_out = './out.treefile', model='MF',bootstrap=1, n_threads = 'AUTO'):
        return_tree = 1
        if phylogeny == 'iqtree':
            iqtree_line = next(path for path in sys.path if 'conda/' in path and '/lib/' in path).split('/lib/')[0]+'/bin/ete3_apps/bin/iqtree'
            subprocess.call('rm %s.ckp.gz'%self.fasta,shell=True)
            subprocess.call(iqtree_line + ' -s %s -m %s -nt %s %s'%(self.fasta,model,n_threads,'-b %d'%int(bootstrap) if int(bootstrap) > 1 else ''), shell=True) #GTR
            shutil.copy(self.fasta+'.treefile',tree_out)
        elif phylogeny == 'phyml':
            phylip_in = self.fasta
            if not phylip_in.endswith('.phylip'):#phylip_in.endswith('.fasta') or phylip_in.endswith('.fa'): #FIXME can add biopython conversion
                #subprocess.call(['perl', fasta2phylip, phylip_in, phylip_in.replace('fasta','phylip').replace('fa','phylip')])
                phylip_in = self.fasta.replace('fasta','phylip').replace('fa','phylip')
                AlignIO.convert(self.fasta,'fasta',phylip_in,'phylip-relaxed')
            phymlLine = next(path for path in sys.path if 'conda/' in path and '/lib/' in path).split('/lib/')[0]+'/bin/ete3_apps/bin/phyml'
            subprocess.call([phymlLine, '-i', phylip_in, '-s', 'BEST', '-q', '-b', bootstrap, '-m', 'GTR'])
            if tree_out:
                shutil.copy(phylip_in+'_phyml_tree.txt',tree_out)
        elif phylogeny == 'fasttree':
            fasttreeLine = next(path for path in sys.path if 'conda/' in path and '/lib/' in path).split('/lib/')[0]+'/bin/ete3_apps/bin/FastTree'
            subprocess.call(fasttreeLine + ' -gtr -nt < ' + self.fasta + ' > ' + tree_out, shell=True)
        else:
            print('Please select different phylogenetic analysis tool [iqtree|phyml|fasttree].')
            return_tree = 0
        if return_tree:
            return TreeObj(tree_out)

class TreeObj:
    def __init__(self, treefile):
        self.treefile = treefile

    def reroot(self, root_species, tree_out):
        t = PhyloTree(open(self.treefile,'r').read())
        t.set_outgroup(root_species)
        t.write(tree_out)
        self.treefile = tree_out

    def tree2matrix(self, distance_matrix = 'distance_matrix.csv'):
        tree = Phylo.read(self.treefile,'newick')
        allclades = list(tree.find_clades(order='level'))
        species_names = [clade.name for clade in allclades if clade.name]
        df = pd.DataFrame(np.nan, index=species_names, columns=species_names)
        for i,j in combinations(species_names,r=2):
            if i == j:
                df.set_value(i,j,0)
            if i != j:
                distance = tree.distance(i,j)
                df.set_value(i,j,distance)
                df.set_value(j,i,distance)
        df_keys = sorted(list(df))
        df = df.reindex(index=df_keys,columns=df_keys)
        df.to_csv(distance_matrix)
        return df

    def output_tree_image(self, fasta_obj, output_image):
        if fasta_obj.fasta.endswith('.fasta') or fasta_obj.fasta.endswith('.fa'):
            subprocess.call("awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' %s > aln.fasta"%fasta_obj.fasta,shell=True)
            t = PhyloTree(self.treefile,alignment='aln.fasta',alg_format='fasta')
        else:
            t = Tree(self.treefile)
        ts = TreeStyle()
        ns = NodeStyle()
        ns['size']=0
        ts.show_leaf_name = True
        ts.show_branch_length = False
        ts.show_branch_support = True
        for n in t.traverse():
            n.set_style(ns)
        t.render(output_image,tree_style = ts, dpi=300)

    def write_trees_intervals(self, interval, out_file):
        with open(self.treefile,'r') as f1, open(out_file,'w') as f2:
            f2.write('\t'.join([interval,f1.read().strip('\n')]))

    def local_trees2final_output(self,work_dir):
        from collections import defaultdict
        import seaborn as sns, matplotlib.pyplot as plt
        from sklearn.manifold import MDS
        import plotly.graph_objs as go
        import plotly.offline as py
        import dendropy
        from dendropy.calculate import treecompare # ete3 calculates this distance as well
        RFDistance = 1
        work_dir += '/'
        cluster_data = defaultdict(list)
        tns = dendropy.TaxonNamespace()
        for interval_data, tree in pd.read_table(self.treefile,header=None,dtype=str).as_matrix().tolist():
                cluster_data[interval_data] = dendropy.Tree.get(data=tree,schema='newick',taxon_namespace=tns)
                cluster_data[interval_data].encode_bipartitions()
        cluster_keys = sorted(cluster_data.keys())
        df = pd.DataFrame(np.nan, index=cluster_keys, columns=cluster_keys)
        for i,j in list(combinations(cluster_keys,r=2)) + zip(cluster_keys,cluster_keys):
            if i == j:
                df.set_value(i,j,0)
            if i != j:
                if RFDistance:
                    dissimilarity = treecompare.weighted_robinson_foulds_distance(cluster_data[i], cluster_data[j])
                else:
                    dissimilarity = np.linalg.norm(cluster_data[i]-cluster_data[j],None)
                df.set_value(i,j,dissimilarity)
                df.set_value(j,i,dissimilarity)
        keys_df = pd.DataFrame(np.array([key.rsplit('_',2) for key in cluster_keys]))
        keys_df[1] = keys_df[1].as_matrix().astype(np.int)
        keys_df = keys_df.sort_values([0,1])
        keys_df[1] = keys_df[1].as_matrix().astype(str)
        new_keys = np.array(['_'.join(x) for x in keys_df.as_matrix()])
        # FIXME sort by chromosome and position, integer... find a way, maybe feed to new dataframe and sort that way, break labels by _ and sort by [0,1] and not [2]
        df = df.reindex(index=new_keys,columns=new_keys)
        df.to_csv(work_dir+'dissimilarity_matrix_local_pca.csv')
        if 0:
            plt.figure()
            sns.heatmap(df)
            plt.savefig(work_dir+'dissimilarity_matrix_local_pca.png',dpi=300)
        local_pca_dissimilarity = df.as_matrix()
        local_pca_dissimilarity = np.nan_to_num(local_pca_dissimilarity)
        mds = MDS(n_components=3,dissimilarity='precomputed')
        transformed_data = mds.fit_transform(local_pca_dissimilarity)
        np.save(work_dir+'local_pca_MDS_transform.npy',transformed_data)
        pickle.dump(list(df.index.values),open(work_dir+'local_pca_window_names.p','wb'))
        plots = []
        plots.append(go.Scatter3d(x=transformed_data[:,0],y=transformed_data[:,1],z=transformed_data[:,2],text=list(df.index.values), mode='markers',marker=dict(color='blue', size=5),name='Regions'))
        py.plot(go.Figure(data=plots),filename=work_dir+'Local_Topology_Differences.html',auto_open=False)




###################
#### CNS CLASS #### # fixme interface circos plots with CNS analysis, plot density of elements

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
@click.option('-q', '--query_prot_id', default = 'all', show_default=True, help='Three number or letter proteome identifier of strain being compared against. If all selected, then synteny graph will be created by running pairwise comparisons between all species.')
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
def run_synteny_pipeline(query_prot_id,fasta_path,synteny_path,gff_path, bed_path, gene_info, fasta_out_dir, bps_threshold, loci_threshold, circos, circos_inputs, circos_outputs):
    """Stitch together many pairwise syntenic blocks into multiple species syntenic blocks and output as fasta files for a multiple sequence alignment. If synteny files are not supplied, conduct pairwise synteny between all included strains."""
    query_protID = query_prot_id
    fasta_files = {fasta.split('_')[-2] : fasta for fasta in glob.glob(fasta_path+'/*.fa')+glob.glob(fasta_path+'/*.fasta')}
    gff_files = {gff.split('.')[-2] : gff for gff in glob.glob(gff_path+'/*.gff')+glob.glob(gff_path+'/*.gff3') }
    intersect_keys = set(fasta_files.keys()) & set(gff_files.keys())
    fasta_files = {protID:fasta for protID,fasta in fasta_files.items() if protID in intersect_keys}
    gff_files = {protID:gff for protID,gff in gff_files.items() if protID in intersect_keys}
    genomes = {}
    pairwise_syntenies = []
    for protID in intersect_keys:
        genomes[protID] = Genome(fasta_files[protID],bed_path+'/'+protID+'.bed',protID,gff_files[protID],gene_info)
        if circos:
            genomes[protID].export_karyotype(circos_inputs+'/'+protID+'.karyotype.txt')
    synteny_files = glob.glob(synteny_path+'/*.unout')+glob.glob(synteny_path+'/*.lifted.anchors')
    synteny_protIDs, synteny_protIDs2 = [], []
    remaining_protIDs, remaining_synteny = [], []
    if synteny_files:
        for synteny_file in synteny_files:
            if synteny_file.endswith('.unout'):
                coords = reduce(lambda x,y: x+y, sorted([[m.start(0),m.end(0)] for m in re.finditer('PAC2_0|PAC4GC',synteny_file)]))[1::2]
                q_protID, s_prot_ID = map(lambda x: synteny_file[x+1:x+4],coords)
            else:
                q_protID, s_prot_ID = tuple(synteny_file[synteny_file.rfind('/')+1:].split('.')[:2])
            if q_protID == query_protID or query_protID == 'all': # fixme for now... in future, implement -global option so all sequences can be included
                synteny_protIDs.append((q_protID, s_prot_ID))
                synteny_protIDs2.extend([(q_protID, s_prot_ID),(s_prot_ID, q_protID)])
                pairwise_synteny = PairwiseSynteny(genomes[q_protID],genomes[s_prot_ID],synteny_file,loci_threshold=loci_threshold)
                pairwise_synteny.generate_synteny_structure(synteny_path)
                pairwise_syntenies.append(pairwise_synteny)

        if len(pairwise_syntenies) < len(intersect_keys)-1 and query_protID != 'all':
            synteny_protIDs = set(np.array(synteny_protIDs)[:,1]).union({query_protID})
            remaining_protIDs = set(intersect_keys) - synteny_protIDs
        elif query_protID == 'all':
            remaining_synteny = set(list(combinations(intersect_keys,r=2))) - set(synteny_protIDs)
        else:
            remaining_protIDs = set()
    else:
        if query_protID == 'all':
            remaining_synteny = list(combinations(intersect_keys,r=2))
        else:
            remaining_protIDs = set(intersect_keys) - {query_protID}


    if list(remaining_protIDs) or remaining_synteny:
        def generate_CDS(protID):
            print(protID)
            genomes[protID].extract_CDS()
            return protID
        p = mp.ProcessingPool()
        r = p.amap(generate_CDS,list(set(reduce(lambda x,y: list(x)+list(y),remaining_synteny))) if query_protID == 'all' else remaining_protIDs.union({query_protID}))
        r.wait()
        protIDs = r.get()
        #for protID in remaining_protIDs.union({query_protID}):
        #    genomes[protID].extract_CDS()

        def p_synteny(protIDs):
            q_protID, s_prot_ID = protIDs
            pairwise_synteny = PairwiseSynteny(genomes[q_protID],genomes[s_prot_ID],loci_threshold=loci_threshold)
            pairwise_synteny.generate_synteny_structure(synteny_path)
            return pairwise_synteny

        r = p.amap(p_synteny,remaining_synteny if query_protID == 'all' else [(query_protID,s_prot_ID) for s_prot_ID in remaining_protIDs])#,callback=mycallback) # _async
        r.wait()
        pairwise_syntenies.extend(r.get())
        p.close()

        """
        for s_prot_ID in remaining_protIDs:
            pairwise_synteny = PairwiseSynteny(genomes[query_protID],genomes[s_prot_ID],loci_threshold=loci_threshold)
            pairwise_synteny.generate_synteny_structure(synteny_path)
            pairwise_syntenies.append(pairwise_synteny)"""

    if circos:
        for pairwise_synteny in pairwise_syntenies:
            pairwise_synteny.synteny_structure_2_link(circos_inputs+'/%s.%s.link.txt'%(pairwise_synteny.q_genome.protID,pairwise_synteny.s_genome.protID))
            circos_obj = Circos(pairwise_synteny)
            circos_obj.generate_config(ticks = circos_inputs+'/txticks.conf', ideogram = circos_inputs+'/txideogram.conf', links_and_rules = circos_inputs+'/linksAndrules.conf', config=circos_inputs+'/circos.conf')
            circos_obj.run_circos(circos_outputs+'/')

    super_synteny = SuperSynteny(pairwise_syntenies, bps_threshold, genomes[query_protID])
    super_synteny.generate_global_synteny_graph(fasta_out_dir)

@joshuatree.command()
@click.option('-f1', '--fasta_1', default = '1.fasta', show_default=True, help='Fasta file 1.', type=click.Path(exists=False))
@click.option('-f2', '--fasta_2', default = '2.fasta', show_default=True, help='Fasta file 2.', type=click.Path(exists=False))
@click.option('-g1', '--gff_1', default = '1.gff', show_default=True, help='GFF file 1.', type=click.Path(exists=False))
@click.option('-g2', '--gff_2', default = '2.gff', show_default=True, help='GFF file 2.', type=click.Path(exists=False))
@click.option('-link', '--link_file', default = '', show_default=True, help='Link file, either lifted anchors or unout.', type=click.Path(exists=False))
@click.option('-info', '--gene_info', default = 'Name', show_default=True, help='Naming convention for gff file\'s gene name field.', type=click.Choice(['Name', 'gene_name']))
@click.option('-l', '--loci_threshold', default= 4, show_default=True, help='Minimum number of genes in a syntenic block in order to include the block.')
@click.option('-w', '--work_dir', default = './', show_default=True, help='Working Directory.')
def extract_syntenic_blocks(fasta_1, fasta_2, gff_1, gff_2, link_file, gene_info, loci_threshold, work_dir):
    """Run pairwise circos in local directory."""
    work_dir += '/'
    genome1 = Genome(fasta_file=fasta_1, bed_file=work_dir+gff_1.split('.')[-2]+'.bed', protID=gff_1.split('.')[-2], gff_file=gff_1, gene_info=gene_info)
    genome2 = Genome(fasta_file=fasta_2, bed_file=work_dir+gff_2.split('.')[-2]+'.bed', protID=gff_2.split('.')[-2], gff_file=gff_2, gene_info=gene_info)
    if not link_file:
        genome1.extract_CDS()
        genome2.extract_CDS()
    pairwise_synteny = PairwiseSynteny(genome1,genome2,link_file,loci_threshold=loci_threshold)
    pairwise_synteny.generate_synteny_structure('./')
    pairwise_synteny.synteny_structure_2_bed(work_dir+'/%s.%s.synteny.bed'%(pairwise_synteny.q_genome.protID,pairwise_synteny.s_genome.protID))

####################
#### RUN CIRCOS ####

@joshuatree.command()
@click.option('-f1', '--fasta_1', default = '1.fasta', show_default=True, help='Fasta file 1.', type=click.Path(exists=False))
@click.option('-f2', '--fasta_2', default = '2.fasta', show_default=True, help='Fasta file 2.', type=click.Path(exists=False))
@click.option('-g1', '--gff_1', default = '1.gff', show_default=True, help='GFF file 1.', type=click.Path(exists=False))
@click.option('-g2', '--gff_2', default = '2.gff', show_default=True, help='GFF file 2.', type=click.Path(exists=False))
@click.option('-link', '--link_file', default = '', show_default=True, help='Link file, either lifted anchors or unout.', type=click.Path(exists=False))
@click.option('-chr1', '--chrom_file1', default = '', show_default=True, help='File listing chromosomes in new order for species 1, can change names of all chromosomes via space delimiting in each line: old_chr_name new_chr_name.', type=click.Path(exists=False))
@click.option('-chr2', '--chrom_file2', default = '', show_default=True, help='File listing chromosomes in new order for species 2, can change names of all chromosomes via space delimiting in each line: old_chr_name new_chr_name.', type=click.Path(exists=False))
@click.option('-info', '--gene_info', default = 'Name', show_default=True, help='Naming convention for gff file\'s gene name field.', type=click.Choice(['Name', 'gene_name']))
@click.option('-l', '--loci_threshold', default= 4, show_default=True, help='Minimum number of genes in a syntenic block in order to include the block.')
@click.option('-n', '--n_chromosomes', default= 25, show_default=True, help='Number of chromosomes in synteny.')
@click.option('-w', '--work_dir', default = './', show_default=True, help='Working Directory.')
@click.option('-v', '--variable_thickness', is_flag=True, help="Variable thickness for the links.")
@click.option('-t', '--thickness_factor', default=1000, show_default=True, help="If variable, thickness of link is length of link divided by factor.")
@click.option('-b', '--bundle_links', is_flag=True, help="Bundle closely spaced links.")
@click.option('-g', '--link_gap', default=10000, show_default=True, help="Gap between closely spaced links.")
@click.option('-s', '--switch_lines', is_flag=True, help="Switch reference and query lines for circos production.")
def pairwise_circos(fasta_1, fasta_2, gff_1, gff_2, link_file, chrom_file1, chrom_file2, gene_info, loci_threshold, n_chromosomes, work_dir, variable_thickness, thickness_factor, bundle_links, link_gap, switch_lines):
    """Run pairwise circos in local directory."""
    work_dir += '/'
    genome1 = Genome(fasta_file=fasta_1, bed_file=work_dir+gff_1.split('.')[-2]+'.bed', protID=gff_1.split('.')[-2], gff_file=gff_1, gene_info=gene_info)
    genome1.export_karyotype(work_dir+fasta_1[:fasta_1.rfind('.')]+'.karyotype.txt',n_chromosomes, chrom_file=chrom_file1)
    genome2 = Genome(fasta_file=fasta_2, bed_file=work_dir+gff_2.split('.')[-2]+'.bed', protID=gff_2.split('.')[-2], gff_file=gff_2, gene_info=gene_info)
    genome2.export_karyotype(work_dir+fasta_2[:fasta_2.rfind('.')]+'.karyotype.txt',n_chromosomes, chrom_file=chrom_file2)
    if not link_file:
        genome1.extract_CDS()
        genome2.extract_CDS()
    pairwise_synteny = PairwiseSynteny(genome1,genome2,link_file,loci_threshold=loci_threshold)
    pairwise_synteny.generate_synteny_structure('./')
    pairwise_synteny.synteny_structure_2_link(work_dir+'/%s.%s.link.txt'%(pairwise_synteny.q_genome.protID,pairwise_synteny.s_genome.protID), bundle_links = bundle_links, link_gap = link_gap)
    circos_obj = Circos(pairwise_synteny)
    circos_obj.generate_config(ticks = work_dir+'./txticks.conf', ideogram = work_dir+'/txideogram.conf', links_and_rules = work_dir+'/linksAndrules.conf', config=work_dir+'/circos.conf', variable_thickness=variable_thickness, thickness_factor=thickness_factor, switch_lines=switch_lines)
    circos_obj.run_circos(work_dir)

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
    for protID in intersect_keys:
        genomes[protID] = Genome(fasta_files[protID],bed_path+'/'+protID+'.bed',protID,gff_files[protID],gene_info)
        genomes[protID].export_karyotype(circos_inputs+'/'+protID+'.karyotype.txt')
    print genomes
    synteny_files = glob.glob(synteny_path+'/*.unout')+glob.glob(synteny_path+'/*.lifted.anchors')
    synteny_protIDs = []
    if synteny_files:
        for synteny_file in synteny_files:
            if synteny_file.endswith('.unout'):
                coords = reduce(lambda x,y: x+y, sorted([[m.start(0),m.end(0)] for m in re.finditer('PAC2_0|PAC4GC',synteny_file)]))[1::2]
                q_protID, s_prot_ID = map(lambda x: synteny_file[x+1:x+4],coords)
            else:
                q_protID, s_prot_ID = tuple(synteny_file[synteny_file.rfind('/')+1:].split('.')[:2])
            synteny_protIDs.extend([(q_protID, s_prot_ID),(s_prot_ID, q_protID)])
            pairwise_synteny = PairwiseSynteny(genomes[q_protID],genomes[s_prot_ID],synteny_file,loci_threshold=loci_threshold)
            pairwise_synteny.generate_synteny_structure(synteny_path)
            pairwise_syntenies.append(pairwise_synteny)
        remaining_synteny = set(list(combinations(intersect_keys,r=2))) - set(synteny_protIDs)
    else:
        remaining_synteny = list(combinations(intersect_keys,r=2))
    if remaining_synteny:
        #print remaining_synteny

        def generate_CDS(protID):
            print(protID)
            genomes[protID].extract_CDS()
            return protID
        p = mp.ProcessingPool(n_cpus)
        r = p.amap(generate_CDS,list(set(reduce(lambda x,y: list(x)+list(y),remaining_synteny))))
        r.wait()
        protIDs = r.get()
        print(protIDs)
        #p = mp.ProcessingPool(ncpus=n_cpus)
        #p.daemon = True
        #r = p.amap(generate_CDS,list(set(reduce(lambda x,y: list(x)+list(y),remaining_synteny)))) # get pathos multiprocessing https://github.com/uqfoundation/pathos
        #while not r.ready():
        #    sleep(5)
        #r.wait()
        """
        for protID in set(reduce(lambda x,y: list(x)+list(y),remaining_synteny)):
            proc = mp.Process(target=lambda: genomes[protID].extract_CDS, args=None)
            proc.daemon = True
            proc.start()
            while len(mp.active_childern()) > n_cpus:
                sleep(1)
        while len(mp.active_childern()) > 0:
            sleep(1)"""
        # fixme add remaining prot ID feauture

        #def mycallback(x):
        #    pairwise_syntenies.extend(x)

        def p_synteny(protIDs):
            print(protIDs)
            q_protID, s_prot_ID = protIDs
            pairwise_synteny = PairwiseSynteny(genomes[q_protID],genomes[s_prot_ID],loci_threshold=loci_threshold)
            pairwise_synteny.generate_synteny_structure(synteny_path)
            return pairwise_synteny

        """
        for q_protID, s_prot_ID in combinations(genomes.keys(),r=2):
            pairwise_synteny = PairwiseSynteny(genomes[q_protID],genomes[s_prot_ID],loci_threshold=loci_threshold)
            pairwise_synteny.generate_synteny_structure(synteny_path)
            pairwise_syntenies.append(pairwise_synteny)"""
        r = p.amap(p_synteny,remaining_synteny)#,callback=mycallback) # _async
        r.wait()
        pairwise_syntenies.extend(r.get())
        p.close()


    for pairwise_synteny in pairwise_syntenies:
        pairwise_synteny.synteny_structure_2_link(circos_inputs+'/%s.%s.link.txt'%(pairwise_synteny.q_genome.protID,pairwise_synteny.s_genome.protID))
        circos_obj = Circos(pairwise_synteny)
        circos_obj.generate_config(ticks = circos_inputs+'/txticks.conf', ideogram = circos_inputs+'/txideogram.conf', links_and_rules = circos_inputs+'/linksAndrules.conf', config=circos_inputs+'/circos.conf')
        circos_obj.run_circos(circos_outputs+'/')

#####################################
#### PAIRWISE SEQUENCE ALIGNMENT ####

@joshuatree.command()
@click.option('-f1', default = './genome1.fa', show_default=True, help='Path containing fasta file one.', type=click.Path(exists=False))
@click.option('-f2', default = './genome2.fa', show_default=True, help='Path containing fasta file two.', type=click.Path(exists=False))
@click.option('-maf', '--out_file', default = './output.maf', show_default=True, help='Output in maf format.', type=click.Path(exists=False))
def pairwise_alignment(f1,f2, out_file):
    """Compute pairwise alignment between two genomes."""
    subprocess.call("samtools faidx %s && samtools faidx %s && lastz --format=maf %s %s > %s"%(f1,f2,f1,f2,out_file),shell=True)

####################
#### RUN CACTUS ####

@joshuatree.command()
@click.option('-fo', '--fasta_output_path', default = './fasta_output/', show_default=True, help='Path containing all syntenic aligned regions, organized into fasta files for multiple sequence alignment via Cactus.', type=click.Path(exists=False))
@click.option('-c', '--cactus_run_directory', default = './cactus_run/', show_default=True, help='Directory containing cactus run information.', type=click.Path(exists=False))
#@click.option('-cac', '--cactus_softlink', default = './runProgressiveCactus.sh', show_default=True, help='Name of softlinked Progressive Cactus batch script file.', type=click.Path(exists=False))
@click.option('-cac', '--cactus_softlink', default = './runProgressiveCactus/', show_default=True, help='Name of softlinked Progressive Cactus distribution.', type=click.Path(exists=False))
#@click.option('-env', '--cactus_env_softlink', default = './cactus_environment', show_default=True, help='Name of softlinked Progressive Cactus virtual environment file.', type=click.Path(exists=False))
@click.option('-n', '--n_cpus', default = 16, show_default=True, help='Number of cpus used to convert hal 2 maf files.')
@click.option('-h2m', '--hal2maf_softlink', default = './hal2maf', show_default=True, help='Name of softlinked Progressive Cactus hal2maf program.', type=click.Path(exists=False))
@click.option('-h2m', '--nickname_file', default = '', show_default=True, help='File containing protID nickname in each line for all protIDs, can omit this file by leaving it blank.', type=click.Path(exists=False))
@click.option('-fi', '--fasta_path', default = './fasta_path/', show_default=True, help='Fasta path containing all of the input genomes. Genome naming must conform to xxx_[protID]_xxx.[fa/fasta].', type=click.Path(exists=False))
@click.option('-s', '--submission_system', default = 'local', show_default=True, help='Different nextflow submission system to use.', type=click.Choice(['local','sge','slurm']))
@click.option('-shift', '--shifter', is_flag = True, help='Use shifter instead of docker.')
def run_cactus(fasta_output_path,cactus_run_directory,cactus_softlink, n_cpus, hal2maf_softlink, nickname_file, fasta_path, submission_system, shifter): #fixme get rid of '' and add to command line tool
    """Run multiple sequence alignment via Progressive Cactus on multiple species synteny blocks and export as maf files. Try to run softlink_cactus beforehand, else use official cactus paths instead of softlinks."""
    cactus_run_obj = CactusRun(fasta_output_path,cactus_run_directory,cactus_softlink, nickname_file, fasta_path)
    cactus_run_obj.write_fastas_seqfile()
    cactus_run_obj.run_cactus(submission_system, shifter=shifter)
    cactus_run_obj.hal2maf(n_cpus, hal2maf_softlink)

@joshuatree.command()
@click.option('-cd', '--cactus_distribution_dir', default = './progressiveCactus/', show_default=True, help='Path containing installed Progressive Cactus Distribution.', type=click.Path(exists=False))
#@click.option('-cac', '--softlink_cactus_name', default = './runProgressiveCactus.sh', show_default=True, help='Name of softlinked Progressive Cactus batch script file.', type=click.Path(exists=False))
@click.option('-cac', '--softlink_cactus_name', default = './runProgressiveCactus', show_default=True, help='Name of softlinked Progressive Cactus distribution.', type=click.Path(exists=False))
#@click.option('-env', '--softlink_env_name', default = './cactus_environment', show_default=True, help='Name of softlinked Progressive Cactus virtual environment file.', type=click.Path(exists=False))
@click.option('-h2m', '--softlink_hal2maf_name', default = './hal2maf', show_default=True, help='Name of softlinked Progressive Cactus hal2maf program.', type=click.Path(exists=False))
def softlink_cactus(cactus_distribution_dir,softlink_cactus_name, softlink_hal2maf_name):
    """Softlink cactus distribution's cactus bash script, virtual environment, and hal2maf program. Useful if installed Cactus to particular directory and want to save time in referencing that directory when running cactus."""
    #subprocess.call('ln -s %s %s'%(os.path.abspath(cactus_distribution_dir+'/bin/runProgressiveCactus.sh'),softlink_cactus_name),shell=True)
    #subprocess.call('ln -s %s %s'%(os.path.abspath(cactus_distribution_dir+'/environment'),softlink_env_name),shell=True)
    subprocess.call('ln -s %s %s'%(os.path.abspath(cactus_distribution_dir),softlink_cactus_name),shell=True)
    subprocess.call('ln -s %s %s'%(os.path.abspath(cactus_distribution_dir+'/submodules/hal/bin/hal2mafMP.py'),softlink_hal2maf_name),shell=True)

@joshuatree.command()
@click.option('-i','--install_path', help='Install Path.', type=click.Path(exists=False))
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

####################
#### MAF Commands

####################
#### VCF Commands

####################
#### Tree Commands

####################
#### Local Tree Topology + Local PCA

@joshuatree.command()
@click.option('-vcf','--vcf_file', help='Input vcf file.', type=click.Path(exists=False))
@click.option('-i','--snps_interval', default = 4000, show_default=True, help='How many snps to use per localized window.', type=click.Path(exists=False))
@click.option('-phy','--phylogeny', default='iqtree', show_default=True, help='Phylogenetic analysis to use.', type=click.Choice(['iqtree','phyml','fasttree']))
@click.option('-w','--work_dir', default = './', show_default = True, help='Work directory for local tree analysis.', type=click.Path(exists=False))
def run_local_trees(vcf_file,snps_interval, phylogeny, work_dir):
    snp = SNP(vcf_file)
    snp.run_local_trees(snps_interval, phylogeny, work_dir)

@joshuatree.command()
@click.option('-vcf','--vcf_file', help='Input vcf file.', type=click.Path(exists=False))
@click.option('-tab','--tab_file', default='out.tab', show_default=True, help='Output tab file.', type=click.Path(exists=False))
def vcf2tab(vcf_file, tab_file):
    snp = SNP(vcf_file)
    snp.vcf2tab(tab_file)

@joshuatree.command()
@click.option('-tab','--tab_file', default='in.tab', show_default=True, help='Input tab file.', type=click.Path(exists=False))
@click.option('-i','--snp_interval', default = 4000, show_default=True, help='How many snps to use per localized window.', type=click.Path(exists=False))
@click.option('-o', '--write_file', default = './output.trees', show_default=True, help='File to write information on local fasta file chunks from local SNPs.',type=click.Path(exists=False))
def tab2chunks(tab_file, snp_interval, write_file):
    snp = SNP(tab_file,'tab')
    snp.tab2chunks(tab_file,snp_interval, write_file)

@joshuatree.command()
@click.option('-f','--fasta_file', default='in.fasta', show_default=True, help='Input fasta file.', type=click.Path(exists=False))
@click.option('-p','--phylogeny', default='iqtree', show_default=True, help='Phylogenetic analysis to use.', type=click.Choice(['iqtree','phyml','fasttree']))
@click.option('-t','--tree_file', default = './out.treefile', show_default=True, help='Output tree file, newick format.', type=click.Path(exists=False))
def generate_phylogeny(fasta_file, phylogeny, tree_file):
    fasta = fasta_aln(fasta_file)
    fasta.generate_tree(phylogeny, tree_out = tree_file, model='GTR',bootstrap=1, n_threads = '1')

@joshuatree.command()
@click.option('-t','--tree_file', default = './in.treefile', show_default=True, help='Input tree file, newick format.', type=click.Path(exists=False))
@click.option('-i','--interval', default = '', show_default=True, help='Bed interval of tree, delimited by underscore.', type=click.Path(exists=False))
@click.option('-o','--out_file', default = 'out.interval', show_default=True, help='Output file containing interval and tree', type=click.Path(exists=False))
def write_trees_intervals(tree_file,interval, out_file):
    TreeObj(tree_file).write_trees_intervals(interval, out_file)

@joshuatree.command()
@click.option('-t','--trees_file', default = './in.treesfile', show_default=True, help='Input trees file, containing intervals and trees.', type=click.Path(exists=False))
@click.option('-w','--work_dir', default = './', show_default = True, help='Work directory for local tree analysis.', type=click.Path(exists=False))
def local_trees2final_output(trees_file,work_dir):
    TreeObj(trees_file).local_trees2final_output(work_dir)



#### RUN CLI ####

if __name__ == '__main__':
    joshuatree()
