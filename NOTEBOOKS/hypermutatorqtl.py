## Load in needed mods
import pandas as pd, numpy as np, scipy.stats as ss

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
colors = ['olive','salmon','tab:blue','tab:green',
          'tab:red','tab:grey','tab:brown',
          'lightgrey','tan','gold','lightblue']

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