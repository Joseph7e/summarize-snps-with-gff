#!/usr/bin/python3


import argparse, sys
from collections import OrderedDict
from Bio.Seq import Seq


def parse_args():
    """ Parse command line arguments """

    # Parse arguments
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-vcf", dest="vcf", type=str, required=True, help="input VCF file produced from GATK mutect2")
    arg_parser.add_argument("-gff", dest="gff", type=str, required=True, help="input gff from RefSeq")
    # arg_parser.add_argument("-cds", dest="cds", type=str, required=True, help="CDS sequences, relies heavily on locus tag")
    arg_parser.add_argument("-report", dest="report", type=str, required=False, help="Assembly report to link chromosome names from ensemble")
    # arg_parser.add_argument("-code", dest="code", type=str, default=1, required=False,
    #                         help="Genetic code, default 1")
    arg_parser.add_argument("-out", dest="out", type=str, default='vcf_data.tsv', required=False,
                            help="Output name to hold data")

    arg_parser.add_argument("-interest", dest="interest", type=str, required=False,
                            help="list of positions of interest")

    return arg_parser.parse_args()

def map_names(report_file):
    """

    :param report_file:
    :return: dictioanry to link vcf names to gff names
    """
    chromosome_dict = {}
    for line in open(report_file):
        if line[0] != '#':
            elements = line.rstrip().split('\t')
            genbank_acc = elements[4].split('.')[0] # these are found in the variant calls
            refseq_acc = elements[6] # these are found in the gff file
            chromosome_dict[genbank_acc] = refseq_acc
    return chromosome_dict


def parse_gff(gff_file):
    """
    :param gff_file:
    :return: gff_lookup, chromosome:{ }
    """
    gff_lookup = {}
    catch = 'protein_id' # make 'locus_tag' for others
    for line in open(gff_file):
        if line[0] != '#':  # skip the headers
            elements = line.rstrip().split('\t')
            chr = elements[0]; start = int(elements[3]); stop = int(elements[4])
            info = elements[8]; direction = elements[6],
            elements_info = info.split(';')
            # for element in elements_info:
            #     if element.startswith(catch):
            #         locus_tag = element.split('=')[-1]
            #locus_tag = info.split('locus_tag=')[-1].split(';')[0]
            gff_lookup.setdefault(chr, {})
            if elements[2] != 'exon' and elements[2] != 'region':
                gff_lookup[chr][str(start)+'_' + str(stop)] = [info, direction, elements[2]]

    return gff_lookup


def parse_vcf(input_file, gff_lookup, output_name, chromosome_map, interest_list):
    cur_count = 0
    with open(output_name, 'w') as out_handle:
        out_handle.writelines('Count\tVCF_contig\tGFF_contig\tPosition\tREF\tALT\tIn_Gene\tStart_Stop\tGene\tDbxref\tProduct\tGene_left\tDbxref\tProduct_left\tLeft_Start_Stop\tDistance\tGene_right\tDbxref\tProduct_right\tRight_Start_Stop\tDistance\n')
        for line in open(input_file):
            if line[0] != '#':  # skip the headers
                cur_count += 1
                if cur_count % 5000 == 0:
                    print ('Completed snp', cur_count)
                elements = line.rstrip().split('\t')
                og_chr = elements[0]; location = int(elements[1]); ref = elements[3]; alt = elements[4]
                if chromosome_map:
                    chr = chromosome_map[og_chr]

                if ',' in ref:
                    ref = ref.split(',')[0]
                if ',' in alt:
                    alt = alt.split(',')[0]

                if alt == '*':
                    alt = ''
                if ref == '*':
                    ref = ''

                # check if it falls in a gene
                in_gene = 'False'
                gene = '-'
                dbx = '-'
                product = '-'
                start_stop = '-'

                closest_gene_left  = '-'
                closest_product_left = '-'
                distance_left = '0'
                left_start_stop = '-'
                closest_gene_right  = '-'
                closest_product_right = '-'
                right_start_stop = '-'
                distance_right = '0'
                right_dbx = '-'
                left_dbx = '-'



                if chr in gff_lookup.keys():
                    for range, gene_info in gff_lookup[chr].items():
                        gff_start, gff_stop = [int(x) for x in range.split('_')]

                        if gff_stop >= location >= gff_start:
                            in_gene = 'True'
                            gene = 'unknown'
                            product = 'unknown'
                            start_stop = (str(gff_start) + ':' + str(gff_stop))
                            for data in gene_info[0].split(';'):
                                key, i = data.split('=')
                                if key == 'gene':
                                    gene = i  + '::' +gene_info[-1]
                                if key == 'product':
                                    product = i
                                if key == 'Dbxref':
                                    dbx = i

                        else: # check to see what the cloesest features are to the left and right
                            if location > gff_stop and location - gff_stop < int(distance_left) or location > gff_stop and distance_left == '0':
                                distance_left = location - gff_start
                                for data in gene_info[0].split(';'):
                                    key, i = data.split('=')
                                    if key == 'gene':
                                       closest_gene_left = i +  '::' + gene_info[-1]
                                    if key == 'product':
                                        closest_product_left = i
                                    if key == 'Dbxref':
                                        left_dbx = i

                                left_start_stop = (str(gff_start) + ':' + str(gff_stop))
                            if location < gff_start and gff_start - location < int(distance_right) or location < gff_start and distance_right == '0':
                                distance_right = gff_start - location
                                for data in gene_info[0].split(';'):
                                    key, i = data.split('=')
                                    if key == 'gene':
                                       closest_gene_right = i + '::' +gene_info[-1]
                                    if key == 'product':
                                        closest_product_right = i

                                    if key == 'Dbxref':
                                        right_dbx = i

                                right_start_stop = (str(gff_start) + ':' + str(gff_stop))

                if in_gene == 'True':
                    closest_gene_left = '-'
                    closest_product_left = '-'
                    distance_left = '-'
                    left_start_stop = '-'
                    closest_gene_right = '-'
                    closest_product_right = '-'
                    right_start_stop = '-'
                    distance_right = '-'
                    right_dbx = '-'
                    left_dbx = '-'

                if cur_count in interest_list:
                    out_handle.writelines('\t'.join([str(x) for x in [cur_count, og_chr, chr, location, ref, alt, in_gene, start_stop, gene, dbx, product,closest_gene_left, left_dbx, closest_product_left, left_start_stop, distance_left, closest_gene_right, right_dbx, closest_product_right, right_start_stop, distance_right]]) + '\n')
                    print ('Im interested')
                else:
                    out_handle.writelines('#' + '\t'.join([str(x) for x in
                                                     [cur_count, og_chr, chr, location, ref, alt, in_gene, start_stop, gene, dbx,
                                                      product, closest_gene_left, left_dbx, closest_product_left, left_start_stop, distance_left, closest_gene_right, right_dbx, closest_product_right, right_start_stop, distance_right]]) + '\n')


if __name__ == '__main__':
    # Parse arguments
    args = parse_args()
    vcf = args.vcf
    gff = args.gff
    output = args.out

    gff_lookup = parse_gff(gff)
    chromosome_lookup = map_names(args.report)

    if args.interest:
        interest_list = []
        for line in open(args.interest):
            interest_list.append(int(line.rstrip()))
        parse_vcf(vcf, gff_lookup, output, chromosome_lookup, interest_list)
    else:
        parse_vcf(vcf, gff_lookup, output, chromosome_lookup, [])



