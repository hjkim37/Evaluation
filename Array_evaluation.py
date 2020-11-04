#!/usr/bin python3

import os
import re
import shlex
import shutil
import subprocess
import multiprocessing as mp
import pandas as pd
import click

@click.group()
def cli():
    """ Performance evaluation using two vcf files for HealthSCAN Array pipeline """

@cli.command()
@click.option('-v', '--vcf', required=True, type=click.Path(exists=True, dir_okay=False), help='input vcf') 
def split_vcf(vcf):
    """ VCF splited by chromosomes
    """
    vcf = os.path.abspath(vcf)
    dirname = os.path.dirname(vcf)
    basename = os.path.basename(vcf)
    name = basename.split('.')[0]
    
    chrom_set = set() 
    with open(vcf, 'r') as f:
        for line in f:
            if ('#' not in line) & ('##' not in line):
                chrom, pos, ID, ref, alt, qual, _filter, info, _format, sample = line.split('\t')
                chrom_set.add(chrom)

    for _chr in chrom_set:
        new_vcf = os.path.join(dirname, f'{name}.{_chr}.vcf')
        cmd = (f'grep -w "^{_chr}" '
               f'{vcf} > {new_vcf}')
        print(cmd)
        subprocess.run(cmd, shell = True)

@cli.command()
@click.option('-s', '--sampleid', required=True, type=click.STRING, help='sample ID')
@click.option('-v', '--input_vcf', required=True, type=click.Path(exists=True, dir_okay=False), help='vcf including information you want to subset')
@click.option('-o', '--output', required=True, type=click.Path(exists=False, dir_okay=True), help='output directory')
def subset_vcf(sampleid, input_vcf, output):
    """ Subset vcf from input_vcf by sample ID """
    
    if not os.path.exists(os.path.join(output, sampleid)):
        os.mkdir(os.path.join(output, sampleid))
    out = os.path.join(output, sampleid, f'{sampleid}.ALL.vcf')
    
    #subset vcf by sample ID
    cmd = (f'bcftools view -s {sampleid} {input_vcf} '
           f'> {out}') 
    print(f'[cmd] {cmd}')       
    subprocess.run(cmd, shell = True)


@cli.command()
@click.option('-v', '--vcf', required=True, type=click.Path(exists=True, dir_okay=False), help='input vcf')
@click.option('-t', '--threads', required=False, type=click.INT, default = 1, help='Number of threads to use')
def vcf2bfile(vcf, threads):
    """ vcf converted to bfile """
    vcf = os.path.abspath(vcf)
    dirname = os.path.dirname(vcf)
    basename = os.path.basename(vcf)
    name = basename.split('.')[0]

    new_bfile = os.path.join(dirname, name)
    cmd = (f'plink --vcf {vcf} '
           f'--make-bed --const-fid '
           f'--out {new_bfile} '
           f'--threads {threads}') 
    subprocess.run(shlex.split(cmd))

@cli.command()
@click.option('-b1', '--bfile1', required=True, type=click.Path(exists=False, dir_okay=False), help='input bfile directory and name W/O extension')
@click.option('-b2', '--bfile2', required=True, type=click.Path(exists=False, dir_okay=False), help='input bfile directory and name W/O extension')
@click.option('-o', '--output', required=True, type=click.Path(exists=False, dir_okay=True), help='output directory')
@click.option('-p', '--prefix', required=True, type=click.STRING, help='Output prefix string')
@click.option('-t', '--threads', required=False, type=click.INT, default = 1, help='Number of threads to use')
def calculate_concordance(bfile1, bfile2, output, prefix, threads):
    """ Two bfiles compare and get different genotype by plink merge mode 7\n
        And calculate genotype concordance"""
    out = os.path.join(output, f'{prefix}.concord')
    diff = get_concordance(bfile1, bfile2, out, threads)
    
    ###Error: multi allele murge fail
    if (diff.returncode != 0) and (os.path.exists(f'{out}.missnp')):
        print('[INFO] Error: get concordance')
        missnp= f'{out}.missnp'

        with open (missnp) as file:
            missnp_count = len(file.readlines())
        print(f'Error: {missnp_count} variants with 3+ alleles present')

        bfile1_exclude = f'{bfile1}_excluded'
        bfile2_exclude = f'{bfile2}_excluded'

        exclude_snps(input_file = bfile1, missnp=missnp, output = bfile1_exclude)
        exclude_snps(input_file = bfile2, missnp=missnp, output = bfile2_exclude) 
        
        ###Merge bfile after exclude missnp
        get_concordance(bfile1_exclude, bfile2_exclude, out, threads)

    result = cal_concordance(f'{bfile2}.bim', f'{out}.diff', prefix)
    final = os.path.join(output, f'{prefix}_result_concordance.xlsx')
    result.to_excel(final, index = False)

def get_concordance(input_bfile, bfile, output, threads):
    cmd = (f'plink --bfile {input_bfile} '
           f'--bmerge {bfile}.bed {bfile}.bim {bfile}.fam '
           f'--merge-mode 7 '
           f'--out {output} '
           f'--threads {threads}')        
    process = subprocess.run(shlex.split(cmd))
    return process

def exclude_snps(input_file, missnp, output):
    print('[INFO] Exclude SNPs Start')
    cmd = (f'plink --bfile {input_file} '
           f'--exclude {missnp} '
           f'--make-bed --out {output}')      
    subprocess.run(shlex.split(cmd))
    print('[INFO] Exclude SNPs End')

def cal_concordance(bim, diff, prefix):
    
    _dict = {
    'sampleID' : list(),
    'SNPs' : list(),
    'diff' : list(),
    'concord' : list()
    }
    
    bim = pd.read_csv(bim, sep = '\t', header = None)
    result = pd.read_csv(diff, sep = '\s+')
    
    snp = len(bim)
    diff = len(result)
    concord = ((snp - diff) / (snp))
    
    _dict['sampleID'].append(prefix)
    _dict['SNPs'].append(snp)
    _dict['diff'].append(diff)
    _dict['concord'].append(concord)
            
    concordance = pd.DataFrame(_dict)
    return concordance

@cli.command()
@click.option('-r', '--standard', required=True, type=click.Path(exists=False, dir_okay=False), help='Basename of Referfence(NA12878) vcf W/O header') #/data/HealthSCAN/Array/NA12878/HG001_GRCh37_intersectedISCAN
@click.option('-i', '--sample', required=True, type=click.Path(exists=True, dir_okay=False), help='input sample vcf') #/data/HealthSCAN/Array/NA12878/iscan_data/CD_20_01588_DE_HS_ARR/CD_20_01588_DE_HS_ARR.intersectedNA12878.vcf
@click.option('-n', '--tn', required=True, type=click.INT, help='variants number of standard')
@click.option('-o', '--output', required=True, type=click.Path(exists=False, dir_okay=True), help='output directory')
@click.option('-p', '--prefix', required=True, type=click.STRING, help='Output prefix string')
#@click.option('-g', '--gender', required=False, default = 'F', type=click.STRING, help='Female is default, add M for male(Y chromosome)')
@click.option('-t', '--threads', required=False, type=click.INT, default = 1, help='Number of threads to use')
def compare_all_vcf_chrom(standard, sample, output, tn, prefix, threads):
    """ Performance evaluation using two vcf files\n
        This requires only correct answer vcf by chromosome and sample vcf\n
        Sample vcf is devided by chromoson automatically\n
        And get True positive, False negative and False positive vcf of all variants\n
    """

    #1~22 chromosome
    cpu = threads
    pool = mp.Pool(processes = cpu) 

    for i in range(1, 23):
        subset_chrom(sample, i)
        out = f'{sample}.{i}'
        ref = f'{standard}.{i}.vcf'
        pool.apply_async(calculate_variants_chrom, args=(ref, out))
    
    pool.close()
    pool.join()
    
    #X chromosome
    subset_chrom(sample, 'X')
    out = f'{sample}.X'
    ref = f'{standard}.X.vcf'
    calculate_variants_chrom(ref, out)

    #Y chromosome
    #if gender == 'M':
    #    subset_chrom(sample, 'Y')
    #    out = f'{sample}.Y'
    #    ref = f'{standard}.Y.vcf'
    #    calculate_variants_chrom(ref, out)
    
    #MT chromosome
    #subset_chrom(sample, 'MT')
    #out = f'{sample}.MT'
    #ref = f'{standard}.MT.vcf'
    #WCcalculate_variants_chrom(ref, out)
    
    #calculate
    outdir = os.path.dirname(sample)
    root, _, files = next(os.walk(outdir))

    TP_files = [f for f in files if re.search(r'vcf\.[0-9XY]+\.TP', f)]
    FP_files = [f for f in files if re.search(r'vcf\.[0-9XY]+\.FP', f)]
    FN_files = [f for f in files if re.search(r'vcf\.[0-9XY]+\.FN', f)]
    #uncalled_files = [f for f in files if f.endswith('uncalled')]
    
    TP_all = os.path.join(output, f'{prefix}.tp')
    FP_all = os.path.join(output, f'{prefix}.fp')
    FN_all = os.path.join(output, f'{prefix}.fn')
    #uncalled_all = os.path.join(output, f'{prefix}.uncall')

    #write all chromosome to one file
    tp = pd.DataFrame()
    for f in TP_files:
        file = os.path.join(root, f)
        df = pd.read_csv(file, sep = '\t')
        tp = pd.concat([tp, df]).reset_index(drop=True)
    tp.to_csv(TP_all, sep = '\t' , index = False)  

    fp = pd.DataFrame()
    for f in FP_files:
        file = os.path.join(root, f)
        df = pd.read_csv(file, sep = '\t')
        fp = pd.concat([fp, df]).reset_index(drop=True)
    fp.to_csv(FP_all, sep = '\t' , index = False)

    fn = pd.DataFrame()
    for f in FN_files:
        file = os.path.join(root, f)
        df = pd.read_csv(file, sep = '\t')
        fn = pd.concat([fn, df]).reset_index(drop=True)
    fn.to_csv(FN_all, sep = '\t' , index = False)  

    #delete chromosome
    delete_files(TP_files, root)
    delete_files(FP_files, root)
    delete_files(FN_files, root)
    #delete_files(uncalled_files, root)
     
    print(f'{prefix}: TP = {len(tp)}, FP = {len(fp)}, FN = {len(fn)}')

    TP = len(tp)       
    FP = len(fp)
    FN = len(fn)
    TN = tn - TP - FP - FN
    
    print(f'[INFO] Performance Evaluation Start')
    #####make final result dataframe #########
    _dict = {
    'sampleID' : list(),
    'TruePositive' : list(),
    'FalsePositive' : list(),
    'FalseNegative' : list(),
    'TrueNegative' : list(),
    'Sensitivity' : list(),
    'Specificity' : list(),
    'PPV' : list(),
    'NPV' : list()
    }
    
    Sensitivity = (TP/(TP + FN)*100) 
    Specificity = (TN/(FP + TN)*100)
    ppv = (TP/(TP + FP)*100)
    npv = (TN/(FN + TN)*100)
    
    _dict['sampleID'].append(prefix)
    _dict['TruePositive'].append(TP)
    _dict['FalsePositive'].append(FP)
    _dict['FalseNegative'].append(FN)
    _dict['TrueNegative'].append(TN)
    _dict['Sensitivity'].append(Sensitivity)
    _dict['Specificity'].append(Specificity)
    _dict['PPV'].append(ppv) 
    _dict['NPV'].append(npv)


    result = pd.DataFrame(_dict)
    output = os.path.join(output, f'{prefix}_evaluation.txt')
    result.to_csv(output, sep = '\t', index = False)  
    print(f'[INFO] Performance Evaluation End')  
 

def calculate_variants_chrom(correct_vcf, sample_vcf):
    print(f'[INFO] Compare vcf Start')
    standard = pd.read_csv(correct_vcf, usecols = [0,1,2,3,4,9], dtype = {0: object}, sep = '\t', header = None)
    standard.columns = ['chr', 'pos', 'id', 'ref', 'alt', 'genotype']
    standard['genotype_tmp'] = standard['genotype'].str.slice(start = 0, stop = 3)
    standard['GQ'] = standard['genotype'].str.slice(start = 4, stop = 9).astype('float64')
    
    for i in standard.index:
        if (standard.loc[i, 'genotype_tmp'] == '0/1') or (standard.loc[i, 'genotype_tmp'] == '1/0'):
            standard.loc[i, 'gt'] = ''.join(sorted(standard.loc[i, 'ref'] + standard.loc[i, 'alt']))
        elif (standard.loc[i, 'genotype_tmp'] == '1/1'):
            standard.loc[i, 'gt'] = standard.loc[i, 'alt'] + standard.loc[i, 'alt']
        elif (standard.loc[i, 'genotype_tmp'] == '0/0'):
            standard.loc[i, 'gt'] = standard.loc[i, 'ref'] + standard.loc[i, 'ref']
        elif (standard.loc[i, 'genotype_tmp'] == '0/2'):
            standard.loc[i, 'gt'] = ''.join(sorted(standard.loc[i, 'ref'] + standard.loc[i, 'alt'].split(',')[1]))
        elif (standard.loc[i, 'genotype_tmp'] == '0/3'):
            standard.loc[i, 'gt'] = ''.join(sorted(standard.loc[i, 'ref'] + standard.loc[i, 'alt'].split(',')[2]))
        else:
            standard.loc[i, 'gt'] = 'NaN'


    sample = pd.read_csv(sample_vcf, usecols = [0,1,2,3,4,9], dtype = {0: object}, sep = '\t', header = None)
    sample.columns = ['chr', 'pos', 'id', 'ref', 'alt', 'genotype']
    sample['genotype_tmp'] = sample['genotype'].str.slice(start = 0, stop = 3)
    sample['GQ'] = sample['genotype'].str.slice(start = 4, stop = 9).astype('float64')

    for i in sample.index:
        if (sample.loc[i, 'genotype_tmp'] == '0/1') or (sample.loc[i, 'genotype_tmp'] == '1/0'):
            sample.loc[i, 'gt'] = ''.join(sorted(sample.loc[i, 'ref'] + sample.loc[i, 'alt']))
        elif (sample.loc[i, 'genotype_tmp'] == '1/1'):
            sample.loc[i, 'gt'] = sample.loc[i, 'alt'] + sample.loc[i, 'alt']
        elif (sample.loc[i, 'genotype_tmp'] == '0/0'):
            sample.loc[i, 'gt'] = sample.loc[i, 'ref'] + sample.loc[i, 'ref']
        elif (sample.loc[i, 'genotype_tmp'] == '0/2'):
            sample.loc[i, 'gt'] = ''.join(sorted(sample.loc[i, 'ref'] + sample.loc[i, 'alt'].split(',')[1]))
        elif (sample.loc[i, 'genotype_tmp'] == '0/3'):
            sample.loc[i, 'gt'] = ''.join(sorted(sample.loc[i, 'ref'] + sample.loc[i, 'alt'].split(',')[2]))
        else:
            sample.loc[i, 'gt'] = 'NaN'

    TP = pd.merge(standard, sample, how = 'inner', on = ['pos', 'gt'])

    FP_tmp = pd.concat([sample, TP]).drop_duplicates(subset = ['pos', 'gt'], keep = False)
    FP_tmp.dropna(axis = 1, how = 'all', inplace = True)
    FP_tmp = FP_tmp[(FP_tmp['genotype_tmp'] == '0/1') | (FP_tmp['genotype_tmp'] == '1/0') | (FP_tmp['genotype_tmp'] == '1/1')]
    FP_tmp = pd.merge(FP_tmp, standard, how='left', on = ['pos'])
    FP = FP_tmp[['GQ_x', 'alt_x', 'chr_x', 'genotype_x', 'genotype_tmp_x', 'gt_x', 'id_x', 'pos', 'ref_x', 'genotype_y', 'gt_y']]
    
    FN_tmp = pd.concat([standard, TP]).drop_duplicates(subset = ['pos', 'gt'], keep = False)
    FN_tmp.dropna(axis = 1, how = 'all', inplace = True)
    FN_tmp = FN_tmp[(FN_tmp['genotype_tmp'] == '0/1') | (FN_tmp['genotype_tmp'] == '1/0') | (FN_tmp['genotype_tmp'] == '1/1')]
    FN_tmp = pd.merge(FN_tmp, sample, how='left', on = ['pos'])
    FN_tmp = FN_tmp[['GQ_x', 'alt_x', 'chr_x', 'genotype_x', 'genotype_tmp_x', 'gt_x', 'id_x', 'pos', 'ref_x', 'genotype_y', 'gt_y']]
 
    #FP와 중복되는 SNP삭제 
    FN = pd.merge(FN_tmp, FP, how='inner', on=['pos'])
    FN = pd.concat([FN, FN_tmp]).drop_duplicates(subset = ['pos'], keep = False)
    FN.dropna(axis = 1, how = 'all', inplace = True)
    FN = FN[['GQ_x', 'alt_x', 'chr_x', 'genotype_x', 'genotype_tmp_x', 'gt_x', 'id_x', 'pos', 'ref_x', 'genotype_y', 'gt_y']]

    #uncalled = sample[sample['genotype_tmp'] == './.']

    print(f'{sample_vcf}:{len(sample)} TP = {len(TP)}, FP = {len(FP)}, FN = {len(FN)}')
    
    True_pos = f'{sample_vcf}.TP'
    False_pos = f'{sample_vcf}.FP'
    False_neg = f'{sample_vcf}.FN'

    TP.to_csv(True_pos, sep = '\t' , index = False)
    FP.to_csv(False_pos, sep = '\t' , index = False)
    FN.to_csv(False_neg, sep = '\t' , index = False)
    #uncalled.to_csv(f'{sample_vcf}.uncalled', sep = '\t' , index = False)
    print(f'[INFO] Compare vcf End')

def subset_chrom(vcf, chrom):
    cmd = (f'grep -w "^{chrom}" '
           f'{vcf} > {vcf}.{chrom}')
    subprocess.run(cmd, shell = True)

def write_all(name, files, root):
    with open(name, 'wb') as combined_f:
        for f in files:
            f_path = os.path.join(root, f)
            with open(f_path, 'rb') as f_obj:
                shutil.copyfileobj(f_obj, combined_f)

def delete_files(files, root):
    for file in files:
        os.remove(os.path.join(root, file))

def count_row(file):
    with open(file, 'r') as f:
        length = len(f.readlines()) - 1
    return length

if __name__ == '__main__':
    cli.main()