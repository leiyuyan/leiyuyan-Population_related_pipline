#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import sys
import time
import warnings
import argparse
import os
import shutil
import tarfile
import gzip
import re
import subprocess
import copy
# from distutils.spawn import find_executable
from shutil import which

def get_gff_dict(gff,up,down):
    gff_dict = {}
    with open(gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip()
            if not line:# not empty line
                continue
            fields = line.split('\t')
            if fields[2] == 'gene':
                chr = fields[0]
                gene_id = fields[8].split(';')[0].split('=')[1]
                start = fields[3]
                end = fields[4]

                if fields[6] == '+':# because vcf 1-base and bcftools 1-base
                    start = int(start) - up if (int(start) - up) >= 0 else 1
                    end = int(end) + down
                elif fields[6] == '-':
                    start = int(start) - down if (int(start) - down) >= 0 else 1
                    end = int(end) + up
                else:
                    print("Error: the 6th column of gff file is not + or -")
                    sys.exit(1)

                gff_dict[gene_id] = [chr, start, end] # gff是1-based,这里是已经扩展过的1-base区域，bcftools获取区间的信息是，使用的1-base 而且左闭右闭
    return gff_dict

def get_base_dict(ref_information,alt_information):
    base_dict = {}
    base_dict['0'] = ref_information.upper()
    for index,base_type in enumerate(alt_information.split(',')):
        base_dict[str(index + 1)] = base_type.upper()
    base_dict['.'] = 'N'
    return base_dict


def pare_vcf(std_bcftools,query_gene_id,outdir):

    sample_haplotype_dict = {}
    chr_pos_list = []
    for line_bcf in std_bcftools.stdout.split('\n'):
        line_bcf = line_bcf.strip()
        if line_bcf:# 最后一个line是''，非''才能操作
            if line_bcf.startswith('##'):
                continue # skip the header line
            elif line_bcf.startswith('#CHROM'):
                header_line = line_bcf # 获取heade with sample information
                sample_list = header_line.split('\t')[9:]
                for sample in sample_list:
                    sample_haplotype_dict[sample] = []
            else:# data line
                fields = line_bcf.split('\t')
                chr_pos_list.append(fields[0]+'_'+fields[1]) # 获取了位置信息
                ref = fields[3]
                alt = fields[4]
                base_dict = get_base_dict(ref,alt)
                for index,information in enumerate(fields[9:]):
                    genotype_list = sorted(re.split(r'[/|]', information.split(':')[0]))# 排过序了，所以0/1是01，1/0是01，0|1是01，1|0是01
                    genotype_base = base_dict[genotype_list[0]] + base_dict[genotype_list[1]]
                    sample_haplotype_dict[sample_list[index]].append(genotype_base)
    haplotype_dict = {}
    for sample in sample_haplotype_dict:
        haplotype = ''.join(sample_haplotype_dict[sample])
        if haplotype not in haplotype_dict:
            haplotype_dict[haplotype] = [[],[]] # 第一个装sample_id，第二个装Haplotype id
            haplotype_dict[haplotype][0].append(sample) # 容易遗漏，这里要加sample_id
        else:
            haplotype_dict[haplotype][0].append(sample)
    id_index = 0
    sum_sample = 0
    for key,value_list in sorted(haplotype_dict.items(),key=lambda x:len(x[1][0]),reverse=True): # sort(dict.items(),key=lambda x:x[1],reverse=True) --->[(key,value),(),]
        id_index += 1
        haplotype_dict[key][1].append('Hap'+str(id_index))
        sum_sample += len(value_list[0])
    print(f'{query_gene_id} total sample num is {sum_sample}')
    # print(haplotype_dict)
    sum_sample = 0
    # 写入文件
    with open(f'{outdir}/{query_gene_id}_haplotype.xls', 'w') as out_f:
        print('sample_id' + '\t' + 'Haplotype' + '\t'  + 'Haplotype_id' + '\t' + 'sample_num' + '\t' + '\t'.join(chr_pos_list), file=out_f)
        for sample in sample_haplotype_dict:
            haplotype = ''.join(sample_haplotype_dict[sample])
            print(sample + '\t' + haplotype + '\t'  + haplotype_dict[haplotype][1][0] + '\t' + str(len(haplotype_dict[haplotype][0])) + '\t' + '\t'.join(sample_haplotype_dict[sample]), file=out_f)
            


if __name__ == '__main__':

    # parameter parse
    my_parser = argparse.ArgumentParser(allow_abbrev=False, add_help=True, prefix_chars='-',description='Help information')  #allow_abbrev: allow to abbreviate argv or not
    my_parser.add_argument('--vcf', '-v', type=str, required=True, action='store', help='the vcf file ')
    my_parser.add_argument('--gff', '-g', type=str, required=True, action='store', help='the gff file')
    my_parser.add_argument('--list', '-l', type=str, default='gene.list', required=True, action='store', help='the gene list file')
    my_parser.add_argument('--outdir', '-o', type=str, default='./', required=True, action='store', help='the outdir')
    my_parser.add_argument('--up', '-u', type=int, default=2000, required=True, action='store', help='the upstream flank length')
    my_parser.add_argument('--down', '-d', type=int, default=0, required=True, action='store', help='the downstream flank length')

    # -------------------------------------------------------------------------
    # Search bcftools、bgzip、tabix
    # -------------------------------------------------------------------------

    # 检查依赖
    for tool in ['bcftools', 'bgzip', 'tabix']:
        if not which(tool):
            print(f"Error: Required tool {tool} not found in PATH.")
            sys.exit(1)
    # -------------------------------------------------------------------------
    # Get args
    # -------------------------------------------------------------------------
    args = my_parser.parse_args()
    vcf = os.path.abspath(args.vcf)
    gff = os.path.abspath(args.gff)
    gene_list = os.path.abspath(args.list)
    outdir = os.path.abspath(args.outdir)
    up_bp = int(args.up)
    d_bp = int(args.down)

    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print("Create output directory:"+ outdir)
    else:
        print("Exist output directory:"+ outdir)


    # -------------------------------------------------------------------------
    # Check vcf file 
    # -------------------------------------------------------------------------
    if not vcf.endswith(".gz"):# not gz file
        print('{vcf} not gz file, start bgzip and tabix and suggest to use the gz file with index file next time to save time')
        try:
            subprocess.run(args=f"bgzip {vcf} && tabix -p vcf {vcf}.gz",shell=True,check=True,stdout=subprocess.PIPE,universal_newlines=True)
        except subprocess.CalledProcessError as e:
            print(f"Error compressing or indexing {vcf}: {e}")
            sys.exit(1)
    else:# gz file
        if  not os.path.exists(f"{vcf}.tbi"):
            try:
                subprocess.run(args=f"tabix -p vcf {vcf}",shell=True,check=True,stdout=subprocess.PIPE,universal_newlines=True)
            except subprocess.CalledProcessError as e:
                print(f"Error indexing {vcf}: {e}")
                sys.exit(1)
       
    # -------------------------------------------------------------------------
    # read data 
    # -------------------------------------------------------------------------
    gff_dict = get_gff_dict(gff,up_bp,d_bp)

    with open(gene_list, 'r') as f_gene:
        for line in f_gene:
            line = line.strip()
            if not line:# not empty line
                continue
            query_gene_id = line
            if query_gene_id in gff_dict:
                query_chr, query_start, query_end = gff_dict[query_gene_id]
                try:
                    # 这里使用的是bcftools view -r Chr01:1000-2000 test.vcf.gz 来获取信息,bcftools 使用1-base且是左闭右开的
                    # bcftools中chr的名称如果是错误的不会报错，这里需要检查，有时间再修改
                    std_bcftools= subprocess.run(args=f"bcftools view -r {query_chr}:{query_start}-{query_end} {vcf}",shell=True,check=True,stdout=subprocess.PIPE,universal_newlines=True)# susprocess 错误输出没有管
                    pare_vcf(std_bcftools,query_gene_id,outdir)
                except subprocess.CalledProcessError as e:
                    print(f"Error getting information {vcf}: {e}")
                    sys.exit(1)
            else:
                print(f"Warning: {query_gene_id} not found in {gff}")
                sys.exit(1)
    #文件排下序
    for file in os.listdir(outdir):
        if file.endswith("_haplotype.xls"):
            # (head -n 1 SC20502G13390_haplotype.xls; tail -n +2 SC20502G13390_haplotype.xls | sort -k4,4nr)
            os.system(f"(head -n 1 {outdir}/{file}; tail -n +2 {outdir}/{file} | sort -k4,4nr) > {outdir}/sort_{file}")

    print("Done")