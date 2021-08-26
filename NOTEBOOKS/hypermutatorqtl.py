## Bring in mods and write ftns for analysis
import pandas as pd, numpy as np, os, glob, scipy.stats as ss
from Bio import SeqIO
from Bio.Seq import Seq

## Bring in and disable warnings
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

## Set centromere dataframe and mat locus
## This was taken from Vikas Yadav's paper in PNAS 2018
## Chromosome 3 was edited here to include a 1 at beginning.
Centromeres = pd.DataFrame([ 
    [970169,835384,1370568,708804,1559983,780649,525714,
    451162,801830,199434,868824,139633,579772,441845],
    [1006931,889427,1409632,752337,1587231,821756,584338,
    512653,839446,243741,933658,171048,632362,477986]]).T

## MAT is on chromosome 5 in H99 reference genome
MAT = np.array([185700,275730])

## Define colormap
colors = ['#E69F00','#009E73','#CC79A7',
              '#56B4E9','#0072B2','#D55E00','#0072B2']

## Set parnetal lables
## For H99
h99_label = 'H99'+'$\mathrm{\u03B1}$'

## Bt65 label
bt_label =  r'Bt65$\bf{a}$'

## Define ftn for conducted a kruskal test
def crypto_kruskal(site,pheno):
    """
    Conducts a bi-allelic comparison at SNP loci coded in SITE on phenotypes PHENO.
    SITE: An array or list coded with either zero or one for a given SNP site.
    PHENO: An array of list with phenotypic values for sample progeny.
    This function assumes that the order of SITE and PHENO (left to right) are paired 
    for each segregant.
    """
    refxpheno = np.array(pheno,dtype=float)[np.array(site)==0]
    altxpheno = np.array(pheno,dtype=float)[np.array(site)==1]
    return -np.log10(ss.kruskal(refxpheno,altxpheno)[1])

## Define ftn for calculating CI
def qtl_boot(segs,gt,df,pheno,pos='Pos',bs=1000,alpha=5):
    """
    Calculate 95 CI for QTL using a Kruskal-W H-test.
    SEGS: labels of segregant names
    GT: Genotype dataframe with rows represent SNP sites along a chromosome and columns represent progeny values.
    DF: Phenotype dataframe with index of progeny names and columns containing phenotype values.
    PHENO: Name of phenotype column in DF used in analysis.
    POS: The positional column in the genotype dataframe containing the positon of a variant on a chromosome.
    BS: Number of iterations used int he bootstraps.
    ALPHA: The significance threshold used in calculating CI. Default is 5 representing an alpha of 0.05.
    """
    
    boot_straped = []

    while len(boot_straped)<bs:
    
        ransegs = np.random.choice(segs,len(segs))
    
        test = gt[ransegs].drop_duplicates().copy()
    
        test['Pval'] = test.apply(crypto_kruskal,
                        args=[df.T[ransegs].T[pheno].values],axis=1)
    
        test = test.T.drop_duplicates().T
        res = gt.merge(test)
        boot_straped.append(res[(res.Pval==res.Pval.max())][pos].median())
        
    return (np.percentile(boot_straped,alpha/2),np.percentile(boot_straped,100 -(alpha/2)))

## Set features of interest
foi = ['gene','three_prime_UTR','five_prime_UTR','CDS']

def addparent(gff,parentsplit='Parent=',idsplit='ID='):
    """
    Adds a parent transcript to gff from attribute field. 
    """
    return [a.split(parentsplit)[-1].split(';')[0].split(idsplit)[-1] 
                      for a in gff.Attribute]

def strandint(gff,strandcol='Strand',svalues = ['-','+'],nsvalues=[-1,1]):
    """
    Converts a strand column in gff from + and - to 1 and -1. 
    """
    return gff[strandcol].replace(dict(zip(svalues,nsvalues)))

## Define ftns for use in QTL mapping
def loadgff(path,
            sep='\t',
            header=None,
            comment='#',
            genecol = 'Gene',
            strandcol = 'Strand',
            parentcol='Parent',
            genesplit = '-t26',
            parentsplit = 'Parent=',idsplit = 'ID=',
            names = ["Seqid", "Source", "Type", "Start",
                     "End", "Score","Strand", "Phase", "Attribute"],
            dtype = ["str","str","str","int","int",
                     "str","str","str","str"],
           svalues = ['-','+'],nsvalues=[-1,1]):
    """
    Loads in and returns a gene features file (GFF) given a specified PATH.
    The GFF is parsed using the seperated specified in SEP.
    A default header of NONE is set and COMMENTs are assumed to start with 
    an # symbol. 
    
    New columns add to the GFF by this function are named in
    GENECOL, STRANDCOL, and PARENTCOL which name the 
    genes, strand, and parent transcripts per gene, respectively. 
    
    The strings set in GENESPLIT, PARENTSPLIT, and IDSPLIT are used
    for splitting the attributes field within the GFF file. 
    
    NAMES: the set of column names of the GFF file.
    DTYPE: the data types set per column in the GFF.
    """
    gff = pd.read_csv(path,comment=comment,sep=sep,header=header,
                  names=names,dtype=dict(zip(names,dtype)))

    gff[strandcol] = strandint(gff,strandcol,svalues,nsvalues)
    gff[parentcol] = addparent(gff,parentsplit,idsplit)
    gff[genecol] = [a.split(genesplit)[0] for a in gff.Parent.tolist()]
    
    return gff


def makeorf(df,fill='.'):
    ref = Seq(''.join(df[(df.Ref!='.')]['Ref'].tolist()))
    alt = Seq(''.join(df[(df.Alt!='.')]['Alt'].tolist()))
    
    if df.Strand.min() < 0:
        ref = ref.reverse_complement()
        alt = alt.reverse_complement()
        
    return ref,alt

def loadref(path,fasta='fasta'):
    return [s for s in SeqIO.parse(path,format=fasta)]


def dotfill(r,a,fill='.'):
    return ''.join([fill for i in range(abs(len(a)-len(r)))])


def scorealign(r,a,fill='.',verbose=False):
    score = 0
    ns = []
    for i in range(len(r))[1:-1]:
        temp = r[:i] + dotfill(r,a) + r[i:]
        c = 0
        for j,n in enumerate(a):
            if n == temp[j]:
                c = c + 1
        if c > score:
            ns = temp
            score = c
            if verbose:
                print(ns,c)
    if len(ns) == 0:
        ns = r[0]+dotfill(r,a)+r[-1]
    return ns


def calcpos(seq,start,end,pad=0,fill='.'):
    c = 0;
    pos = np.arange(start,end+pad)
    poss= []
    if fill in seq:
        for n in seq:
            if n != fill:
                poss.append(pos[c])
                c = c + 1
            else:
                poss.append(pos[c])
    else:
        poss = pos
    return poss


def realign(ref,alt,verbose=False):
    if len(alt) > len(ref):
        newref = scorealign(ref,alt,verbose=verbose)
        newalt = alt
    elif len(ref) > len(alt):
        newalt = scorealign(alt,ref,verbose=verbose)
        newref = ref
    else:
        assert len(alt) == len(ref)
        newref = ref
        newalt = alt
    return newref,newalt

rdict = dict(zip(['three_prime_UTR',
                  'five_prime_UTR',
                  'CDS'],
                 [3,5,0]))


def gvtodf(gv,refseqdf,gffdf,rs='REF',fill='.',rdict=rdict,verbose=False):

    refs = []
    alts = []
    poss = []
    for i,j in gv.iterrows():
    
        ref = j.Ref
        alt = j.Alt.split('.')[int(j.Sample)]
        start = j.Pos
        end = start + len(ref)
    
        newref,newalt = realign(ref,alt)
        newpos = calcpos(newref,start,end)
    
        refs.append(newref)
        alts.append(newalt)
        poss.append(newpos)
    
    pos = np.concatenate(poss)
    ref = [n for n in ''.join(refs)]
    alt = [n for n in ''.join(alts)]

    df = pd.DataFrame([pos,ref,alt],
                      index=['Pos','Ref','Alt']).T
    df['Isvar'] = 1

    check = refseqdf[(refseqdf.Pos.isin(df.Pos))]
    failed = []
    for i,j in check.iterrows():
        
        tocheck = [s for s in 
                   df[(df.Pos==j.Pos)].Ref.tolist() 
                   if (s != fill)]
        if len(tocheck) > 1: ## Should be one base
            failed.append(j.Pos)
        elif not tocheck[0] == j[rs]:
            failed.append(j.Pos)
            
    if (len(failed) > 0) & verbose:
        print('Reference alignment failed')
    
    dfadd = refseqdf[~(refseqdf.Pos.isin(df.Pos))].copy()
    dfadd['Ref'] = dfadd[rs]
    dfadd['Alt'] = dfadd[rs]
    dfadd['Isvar'] = 0
    dfadd.drop(rs,axis=1,inplace=True)

    res = pd.concat([dfadd,df]).sort_values('Pos')
    res['Strand'] = gffdf.Strand.min()
    res['Type'] = -1
    res['Phase'] = -1

    for i,j in gffdf[(gffdf.Type!='gene')].iterrows():
    
        nb = np.arange(j.Start,j.End+1)
        res.loc[(res.Pos.isin(nb)),'Type'] = j.Type
    
        if j.Type == 'CDS':
            res.loc[(res.Pos.isin(nb)),'Phase'] = int(j.Phase)
    
    assert -1 not in res[(res.Type=='CDS')].Phase.tolist()
    
    res.Type.replace(rdict,inplace=True)

    return res,failed

def genesummary(nv,nt,fa,savepath,
                fname = 'No_snpeffect_results.csv',
                reasons = ['Invariant','NoAnnotation',
                           'FailedAlignment'],
                column_names = ['Parent','Reason']):
    """
    Writes a report of variant prediction resutls from GVTODF function
    """
    genes = [nv,nt,fa] ## Group gene names
    report = [] ## initilize report

    for i,j in enumerate(genes): ## Iterate over each gene

        ## Zip and append gene names with reason for no effect
        nz = list(zip(j,np.repeat(reasons[i],len(j)))) 
        report.append(nz) 
    #print(report)
    ## Save gene names and reasons for report
    finrep = np.concatenate([r for r in report if len(r)>0])
    pd.DataFrame(finrep, 
                 columns = column_names
                 ).to_csv(savepath+fname,index=False)
    
    return savepath+fname

def snpeffect(samplename,genos,savepath, gff, REF,
                zerobase=False, skip=False, verbos=False, 
                info_cols = None, parents = None,
                Start_column = 'Start', 
                Sample_column = 'Sample', 
                Gene_column = 'Gene',
                variant_column = 'Isvar',
                position_column = 'Pos',
                reference_column = 'REF',
                fasta='fasta',
                varindicator=1):
    """
    Imputs effects of genetic variants on annotated sequence of genes. 
    """
    ## initilize lists
    novars = [] 
    notype = []
    failalign = []

    ## Gather infocols
    if (info_cols is None):
        info_cols = [c for c in genos.columns if c != samplename]

    ## Make path if not there
    if not os.path.exists(savepath):
        os.mkdir(savepath)

    ## Load GFF
    if (type(gff) == str):
        gff = loadgff(gff)
        gff = gff[(gff.Type.isin(foi))]

    elif (type(gff) == pd.core.frame.DataFrame) and ('Parent' in gff.columns.tolist()):
        pass
    else:
        print('Please specify a GFF or a proper path to a GFF with a column denoting the parent transcript per gene.\n')

    ## Load reference
    if (type(REF) == str):
        REF = loadref(REF,fasta)
    elif (type(REF) == list) or (type(REF) == np.array):
        pass
    else:
        print('Please specify a reference file or path to file\n')

    if (parents is None):
        parents = np.unique(np.concatenate([g[(g.Type!='gene') & (g.Type.isin(foi))].Parent.unique() for i,g in gff.groupby('Gene')]))
    
        if verbos:
            print(len(parents))

    ## Conduct analysis per gene
    for genename in parents: ## iterate over the gene parent names
        genesave = savepath+genename+'.csv.gz' ## Set save path
        if os.path.exists(genesave) and skip: ## check if df exists
            continue

    ## Take gene info from gff file
        genedf = gff[(gff.Type.isin(foi)) & (gff.Parent==genename)].sort_values(Start_column).copy()
        genedf[Start_column] = genedf[Start_column] + (1 if zerobase else 0)
        
        if not genedf.Start.min() == gff[(gff.Parent==genename)].Start.min() + (1 if zerobase else 0):
            print('start error in %s.\n%s'%(genename,genedf.Start.min()-gff[(gff.Parent==genename)].Start.min()))
        
    ## Check work
        assert len(genedf.Strand.unique()) <= 1, "More than one strand on gene %s."%genename
        if 'CDS' not in genedf.Type.unique():#(len(genedf.Type.unique()) >= 3):
            notype.append(genename)
            continue
            
    ## For a sample of interest bring in data
        sample = genos.loc[genos[samplename].dropna().index,info_cols+[samplename]]
        sample[Sample_column] = sample[samplename]
    
    ## Locate genetic variants within a gene
        gv = sample[(sample.Seqid==genedf.Seqid.min()) & (sample.Pos>=genedf.Start.min()) & (sample.Pos<=genedf.End.max())]
        
        if (gv.shape[0] == 0):
            novars.append(genename) ## append the list genes with zero variants
            if verbos:
                print('No variants within %s for strain %s\n'%(genename,samplename))
            continue
    
    ## From the referecne take the sequence of interest
        gene_seq = [s.seq[genedf.Start.min()-(0 if zerobase else 1):genedf.End.max()] 
                                    for s in REF if s.id == genedf.Seqid.min()][0]
    
    ## Make a dataframe from these sequences
        refseqdf = pd.DataFrame([np.arange(genedf.Start.min()+(1 if zerobase else 0),
                                genedf.End.max()+1), [n for n in str(gene_seq)]],
                                index = [position_column,reference_column]).T;
        #print(refseqdf.shape,genedf.End.max()-genedf.Start.min())

    ## Make gene variant dataframe
        res,failed  = gvtodf(gv,refseqdf,genedf,verbose=verbos)
        res[Gene_column] = genename
        res[Sample_column] = samplename

        if res.Pos.min() != genedf.Start.min() and verbos:
            print("Position mismatch at start: %s %s %s\n"%(genename,res.Pos.min(),genedf.Start.min()))
        if res.Pos.max() != genedf.End.max() and verbos:
            print("Position mismatch at end: %s %s %s \n"%(genename,res.Pos.max(),genedf.End.max()))

    ## Append failed alignments
        if (len(failed) > 0):
            failalign.append(genename)
    
    ## Recode variant index
        res.reset_index(drop=True, inplace=True)
        res[variant_column] = 0 
        revar_ix = res[(res.Ref != res.Alt)].index.values
        res.loc[revar_ix,variant_column] = varindicator

    ## Recheck variants
        if (res[(res[variant_column]==varindicator)].shape[0]==0):
            novars.append(genename)
            continue

    ## Save
        if verbos:
            print(genesave)

        res.to_csv(genesave, index=False)
        if verbos:
            print('%s\n-------------\n%s\n%s\n-------------\n%s\n'%(
                genename,res[position_column].unique().shape[0],
                res.shape[0],res[(res.Isvar==varindicator)].shape[0]))
        else:
            pass
    
    #print(np.sum([len(novars),len(notype),len(failalign)]))
    ## Save summary results of genes we failed to analyze
    if (np.sum([len(novars),len(notype),len(failalign)]) == 0):
        respath = []
    else:
        respath = genesummary(novars,notype,failalign,savepath)

    return respath


def snpresults(samplepath,
                file_end = '.csv.gz',
                parent_end = '-t26',
                column_names = ['Strain',
                                'Gene','Parent','Expected',
                                'Ref','Alt','Refstop','Altstop','Nonsyn',
                                'Nvars','Utr5','Utr3',
                                'Exon','Intron','Lengthdif','Strand'],
                tostop = True,
                stopchr = '*'):
    """
    Summarize results of genes containing genetic variants
    """
    assert samplepath[-1]=='/', "Incomplete sample path be sure to add an / to end."

    genepath = sorted(glob.glob(samplepath+'*'+file_end))

    assert len(genepath) > 0, "No gene files detected, check path."

    todf = []
    for s in genepath:

        genep = s.split('/')[-1].split(file_end)[0]
        gene_name = s.split('/')[-1].split(parent_end)[0]

        temp = pd.read_csv(s)
        gene = temp.Gene.min()
        sample = temp.Sample.min()
        strand = temp.Strand.min()

        #assert sample == temp.Sample.min()
        #assert gene == genep
        
        cds = temp[(temp.Type==0)]
        
        ref,alt = makeorf(cds)
        
        el = (cds.Pos.unique().shape[0]/3) - 1 ## Removes stop-codon
        
        ra = ref.translate(to_stop=tostop)
        aa = alt.translate(to_stop=tostop)
        
        rl = len(ra)
        al = len(aa)
        
        sr = ref.translate().count(stopchr)
        sa = alt.translate().count(stopchr)
        
        r = ref.translate()
        a = alt.translate()
        res = np.unique(list(r)+list(a))

        ## Count the number on non-synonmous changes
        ns = np.sum(np.abs(np.array([r.count(rs) for rs in res]) - np.array([a.count(rs) for rs in res])))
        
        nvars = temp[(temp.Isvar==1)].shape[0]
        utr3 = temp[(temp.Type==3) & (temp.Isvar==1)].shape[0]
        utr5 = temp[(temp.Type==5) & (temp.Isvar==1)].shape[0]
        inexon = temp[(temp.Type==0) & (temp.Isvar==1)].shape[0]
        inintron = temp[(temp.Type==-1) & (temp.Isvar==1)].shape[0]
        
        todf.append((sample,gene_name,gene,
                     el,rl,al,sr,sa,ns,
                     nvars,utr5,utr3,inexon,inintron,np.abs(rl-al)/rl,strand))
    
    return pd.DataFrame(todf, columns=column_names)