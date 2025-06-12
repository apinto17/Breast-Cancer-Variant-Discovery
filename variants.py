import matplotlib.pyplot as plt
from cyvcf2 import VCF
from collections import defaultdict
import pandas as pd

gtf_file = 'gencode.gtf.gz'


def plot_variant_density():
    bin_size = 1_000 
    vcf = VCF("output.vcf.gz") # returns a generator

    # find the frequencies for every variant bin of 1000 bases
    freq = defaultdict(int)

    for variant in vcf():
        # only check for passing genes
        if variant.FILTER not in (None, 'PASS'):
            continue
        # NOTE: comment back in to check for more than one allele
        # check for more than one alternative allele - optional
        # if len(variant.ALT) <= 1:
        #     continue

        bin_index = variant.POS // bin_size
        key = (variant.CHROM, bin_index)
        freq[key] += 1

    # get the top 10 regions
    top_10 = sorted(freq.items(), key=lambda x: x[1], reverse=True)[:10]

    for (chrom, bin_index), count in top_10:
        start = bin_index * bin_size
        end = start + bin_size
        print(f"Chromosome: {chrom}, Region: {start}-{end}, Variants: {count}")

    # plot the top 10 regions
    labels = [f"{i * bin_size}-{(i + 1) * bin_size}" for (_, i), _ in top_10]
    counts = [count for _, count in top_10]

    plt.figure(figsize=(10, 6))
    plt.bar(labels, counts, color='skyblue')
    plt.xticks(rotation=45, ha='right')
    plt.xlabel("Genomic Region (bp)")
    plt.ylabel("Variant Count")
    plt.title("Top 10 Variant-Dense Genomic Regions")
    plt.tight_layout()
    plt.savefig("genomic_variant_regions.png")


def parse_gtf(gtf_path):
    # load only relevant columns, use tab delimited seperator
    gtf = pd.read_csv(gtf_path, sep='\t', comment='#', header=None,
                      names=['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    # filter for "gene" feature to get gene-level regions
    genes = gtf[gtf['feature'] == 'gene'].copy()

    # extract gene_name from attribute field
    genes['gene_name'] = genes['attribute'].str.extract('gene_name "([^"]+)"')
    return genes[['chrom', 'start', 'end', 'gene_name']]


def find_overlapping_genes(chrom, pos, genes_by_chrom):
    if chrom not in genes_by_chrom:
        return []
    genes_chr = genes_by_chrom[chrom]
    overlapping = genes_chr[(genes_chr['start'] <= pos) & (genes_chr['end'] >= pos)]
    return overlapping['gene_name'].tolist()


def get_variants_per_gene():
    genes_df = parse_gtf(gtf_file)
    # get all the genes grouped by chromosome
    genes_by_chrom = {}
    for chrom, group in genes_df.groupby('chrom'):
        genes_by_chrom[chrom] = group.sort_values('start') # sort by start regions

    # count variants per gene
    variant_counts = {}
    vcf = VCF("output.vcf.gz")

    for variant in vcf:
        chrom = variant.CHROM
        pos = variant.POS
        # find overlapping genes per chromosome and position
        overlapping_genes = find_overlapping_genes(chrom, pos, genes_by_chrom)
        for gene in overlapping_genes:
            variant_counts[gene] = variant_counts.get(gene, 0) + 1

    print(f"Total genes with variants: {len(variant_counts)}")
    return variant_counts

def plot_variants_per_gene(variant_counts):
    top_genes = sorted(variant_counts.items(), key=lambda x: x[1], reverse=True)[:20]
    genes, counts = zip(*top_genes) # separate genes and counts

    # plot variants per gene
    plt.figure(figsize=(10,6))
    plt.barh(genes[::-1], counts[::-1])  # reverse for descending order top->bottom
    plt.xlabel('Variant Count')
    plt.title('Top 20 Genes by Variant Count')
    plt.tight_layout()
    plt.savefig("variants_per_gene.png")


if __name__ == "__main__":
    plot_variant_density()
    variant_counts = get_variants_per_gene()
    plot_variants_per_gene(variant_counts)
