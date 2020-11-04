#!/data/shared_env/bin python3

import os
import re
import shlex
import shutil
import subprocess
import multiprocessing as mp
import pandas as pd
import click

#NA12878_BED = '/data/users/hjkim/NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed'
#NA12878_VCF = '/data/users/hjkim/NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf'
#우리가 가지고 있는 파일에는 chr문자가 들어가 있어서 ref정보를 수정해 줌 
#NA12878_BED = '/data/users/hjkim/NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_rev.bed'
#NA12878_VCF = '/data/users/hjkim/NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer_rev.vcf'
NA12878_BED = '/data/HealthSCAN/NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_rev.bed'
NA12878_VCF = '/data/HealthSCAN/NA12878/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer_rev.vcf'

@click.group()
def cli():
    """ Performance evaluation using standard matarials """

@cli.command()
@click.option('-b', '--input_bed', required=True, type=click.Path(exists=True, dir_okay=False), help='bed file including region information you want to subset')
@click.option('-v', '--input_vcf', required=True, type=click.Path(exists=True, dir_okay=False), help='vcf file you want to subset')
@click.option('-t', '--tn', required=True, type=click.INT, help='variants number of standard')
@click.option('-o', '--output', required=True, type=click.Path(exists=False, dir_okay=True), help='output directory')
@click.option('-p', '--prefix', required=True,  type=click.STRING, help='Output prefix string. e.g. output file: prefix.intersectedNA12878.bed')
def performance(input_bed, input_vcf, output, prefix  ):
    """ Performance evaluation using standard materials(NA12878) 
        This requires sample bed, vcf files 
        And get True positive, False negative and False positive  
    """
    intersected_bed, intersected_vcf = intersected_bed(input_bed, output, prefix)
    sample_vcf = intersect_vcf_bed(input_vcf, input_bed, output, prefix)
    compare_vcf(intersected_vcf, sample_vcf, output, tn, prefix)


@cli.command()
@click.option('-b', '--input_bed', required=True, type=click.Path(exists=True, dir_okay=False), help='bed file including region information you want to subset')
@click.option('-o', '--output', required=True, type=click.Path(exists=False, dir_okay=True), help='output directory')
@click.option('-p', '--prefix', required=True,  type=click.STRING, help='Output prefix string. e.g. output file: prefix.intersectedNA12878.bed')
def intersect_bed(input_bed, output, prefix):
    """ Intersect sample bed with NA12878 bed of Standard Materials(HG001) 
        And get intersected NA12878 VCF with output bed
    """
    out_bed = os.path.join(output, f'{prefix}.intersectedNA12878.bed')
    cmd = (f'bedtools intersect -a {input_bed} '
           f'-b {NA12878_BED} '
           f'> {out_bed}') 
    print(f'[cmd] {cmd}')       
    subprocess.run(cmd, shell = True)

    bed = pd.read_csv(out_bed, sep = '\t', header = None)
    bed.columns = ['chr', 'start', 'end']
    bed['length'] = bed['end'] - bed['start'] + 1
    bed_lengh = bed['length'].sum()

    print(f'bed: {out_bed} \n bed_lengh: {bed_lengh} \n')

    out_vcf = os.path.join(output, f'{prefix}.intersectedNA12878.vcf')
    print(f'Intersect NA12878 VCF with {out_bed}')
    cmd = (f'bedtools intersect -a {NA12878_VCF} '
           f'-b {out_bed} '
           f'> {out_vcf}') 
    print(f'[cmd] {cmd}')       
    subprocess.run(cmd, shell = True)
    
    return out_bed, out_vcf

@cli.command()
@click.option('-v', '--input_vcf', required=True, type=click.Path(exists=True, dir_okay=False), help='vcf file you want to subset')
@click.option('-b', '--input_bed', required=True, type=click.Path(exists=True, dir_okay=False), help='bed file including region information you want to subset')
@click.option('-o', '--output', required=True, type=click.Path(exists=False, dir_okay=True), help='output directory')
@click.option('-p', '--prefix', required=True,  type=click.STRING, help='Output prefix string')
def intersect_vcf_bed(input_vcf, input_bed, output, prefix):
    """ Intersect sample VCF with sample bed """
   
    out_vcf = os.path.join(output, f'{prefix}.intersected.vcf') 
    cmd = (f'bedtools intersect -a {input_vcf} '
           f'-b {input_bed} '
           f'> {out_vcf}') 
    print(f'[cmd] {cmd}')       
    subprocess.run(cmd, shell = True)

    return out_vcf


@cli.command()
@click.option('-r', '--standard', required=True, type=click.Path(exists=True, dir_okay=False), help='Referfence(NA12878) vcf W/O header')
@click.option('-i', '--sample', required=True, type=click.Path(exists=True, dir_okay=False), help='input sample vcf W/O header')
@click.option('-t', '--tn', required=True, type=click.INT, help='variants number of standard')
@click.option('-o', '--output', required=True, type=click.Path(exists=False, dir_okay=True), help='output directory')
@click.option('-p', '--prefix', required=True, type=click.STRING, help='Output prefix string')
def compare_vcf(standard, sample, output, tn, prefix):
    """ Performance evaluation \n
        This requires only correct answer vcf and sample vcf with no header\n
        The VCF files must have whole 9 cloummns\n
        Finally get True positive, False negative and False positive vcf of all variants \n
    """

    print(f'[INFO] Compare vcf Start')
    standard = pd.read_csv(standard, usecols = [0,1,2,3,4,9], dtype = {0: object}, sep = '\t', header = None)
    standard.columns = ['chr', 'pos', 'rsid', 'ref', 'alt', 'genotype']
    standard['genotype_tmp'] = standard['genotype'].str.slice(start = 0, stop = 3)
    for i in standard.index:
        if (standard.loc[i, 'genotype_tmp'] == '0|1') or (standard.loc[i, 'genotype_tmp'] == '1|0') or (standard.loc[i, 'genotype_tmp'] == '0/1') or (standard.loc[i, 'genotype_tmp'] == '1/0'):
            standard.loc[i, 'gt'] = ''.join(sorted(standard.loc[i, 'ref'] + standard.loc[i, 'alt']))
        elif (standard.loc[i, 'genotype_tmp'] == '1|2') or (standard.loc[i, 'genotype_tmp'] == '2|1') :
            standard.loc[i, 'gt'] = ''.join(sorted(standard.loc[i, 'alt'].split(',')[0] + standard.loc[i, 'alt'].split(',')[1]))
        #elif (standard.loc[i, 'genotype_tmp'] == '2|1'):
        #    standard.loc[i, 'gt'] = standard.loc[i, 'alt'].split(',')[1] + standard.loc[i, 'alt'].split(',')[0]
        elif (standard.loc[i, 'genotype_tmp'] == '1|1') or (standard.loc[i, 'genotype_tmp'] == '1/1') :
            standard.loc[i, 'gt'] = standard.loc[i, 'alt'] + standard.loc[i, 'alt']
        elif (standard.loc[i, 'genotype_tmp'] == '0|0') or (standard.loc[i, 'genotype_tmp'] == '0/0') :
            standard.loc[i, 'gt'] = standard.loc[i, 'ref'] + standard.loc[i, 'ref']   
        else:
            standard.loc[i, 'gt'] = 'NaN'


    sample = pd.read_csv(sample, usecols = [0,1,2,3,4,9], dtype = {0: object}, sep = '\t', header = None)
    sample.columns = ['chr', 'pos', 'id', 'ref', 'alt', 'genotype']
    sample['genotype_tmp'] = sample['genotype'].str.slice(start = 0, stop = 3)
    #sample['GQ'] = sample['genotype'].str.slice(start = 4, stop = 9).astype('float64')
    sample['GQ'] = sample.genotype.str.split(':').str[3]

    for i in sample.index:
        if (sample.loc[i, 'genotype_tmp'] == '0/1') or (sample.loc[i, 'genotype_tmp'] == '1/0'):
            sample.loc[i, 'gt'] = ''.join(sorted(sample.loc[i, 'ref'] + sample.loc[i, 'alt']))
        elif (sample.loc[i, 'genotype_tmp'] == '1/1'):
            sample.loc[i, 'gt'] = sample.loc[i, 'alt'] + sample.loc[i, 'alt']
        elif (sample.loc[i, 'genotype_tmp'] == '0/0'):
            sample.loc[i, 'gt'] = sample.loc[i, 'ref'] + sample.loc[i, 'ref']
        elif (sample.loc[i, 'genotype_tmp'] == '1/2') or (sample.loc[i, 'genotype_tmp'] == '2/1'):
            sample.loc[i, 'gt'] = ''.join(sorted(sample.loc[i, 'alt'].split(',')[0] + sample.loc[i, 'alt'].split(',')[1]))
        elif (sample.loc[i, 'genotype_tmp'] == '0/2'):
            sample.loc[i, 'gt'] = ''.join(sorted(sample.loc[i, 'ref'] + sample.loc[i, 'alt'].split(',')[1]))
        elif (sample.loc[i, 'genotype_tmp'] == '0/3'):
            sample.loc[i, 'gt'] = ''.join(sorted(sample.loc[i, 'ref'] + sample.loc[i, 'alt'].split(',')[2]))
        else:
            sample.loc[i, 'gt'] = 'NaN'

    TP = pd.merge(standard, sample, how = 'inner', on = ['chr', 'pos', 'gt'])

    FP_tmp = pd.concat([sample, TP]).drop_duplicates(subset = ['chr', 'pos', 'gt'], keep = False)
    FP_tmp.dropna(axis = 1, how = 'all', inplace = True)
    FP = FP_tmp[(FP_tmp['genotype_tmp'] == '0/1') | (FP_tmp['genotype_tmp'] == '1/0') | (FP_tmp['genotype_tmp'] == '1/1')]

    FN_tmp = pd.concat([standard, TP]).drop_duplicates(subset = ['chr', 'pos', 'gt'], keep = False)
    FN_tmp.dropna(axis = 1, how = 'all', inplace = True)
    #FP와 중복되는 SNP삭제 
    if len(FP) != 0 and len(FN_tmp) != 0:
        tmp = pd.merge(FN_tmp, FP, on = ['chr', 'pos'], how = 'inner')
        FN = pd.concat([FN_tmp, tmp]).drop_duplicates(subset = ['chr', 'pos'], keep = False)
    else:
        FN = FN_tmp

    uncalled = sample[sample['genotype_tmp'] == './.']

    
    True_pos = os.path.join(output, f'{prefix}.TP')
    False_pos = os.path.join(output,f'{prefix}.FP')
    False_neg = os.path.join(output,f'{prefix}.FN')
    Uncalled =  os.path.join(output,f'{prefix}.uncalled')

    
    #결과 저장 

    TP.to_csv(True_pos, sep = '\t' , index = False)
    if len(FP) != 0:
        FP.to_csv(False_pos, sep = '\t' , index = False)
    if len(FN) != 0:   
        FN.to_csv(False_neg, sep = '\t' , index = False)
    if len(uncalled) != 0:
        uncalled.to_csv(Uncalled, sep = '\t' , index = False)

    tp = len(TP)
    fp = len(FP)
    fn = len(FN)
    TN = tn - tp -fp - fn

    print(f'TP = {len(TP)}, FP = {len(FP)}, FN = {len(FN)}')
    print(f'[INFO] Compare vcf End')
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
    
    Sensitivity = (tp/(tp + fn)*100) 
    Specificity = (TN/(fp + TN)*100)
    ppv = (tp/(tp + fp)*100)
    npv = (TN/(fn + TN)*100)
    

    #sampleID =  os.path.basename(sample)
    _dict['sampleID'].append(prefix)
    _dict['TruePositive'].append(tp)
    _dict['FalsePositive'].append(fp)
    _dict['FalseNegative'].append(fn)
    _dict['TrueNegative'].append(TN)
    _dict['Sensitivity'].append(Sensitivity)
    _dict['Specificity'].append(Specificity)
    _dict['PPV'].append(ppv) 
    _dict['NPV'].append(npv)

    result = pd.DataFrame(_dict)
    output = os.path.join(output, f'{prefix}_evaluation.txt')
    result.to_csv(output, sep = '\t', index = False)   
    print(f'[INFO] Performance Evaluation End')     




if __name__ == '__main__':
    cli.main()