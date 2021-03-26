#!/data/shared_env/bin python3

import os
import subprocess
from multiprocessing.pool import ThreadPool as Pool
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
def performance(input_bed, input_vcf, tn, output, prefix):
    """ Performance evaluation using standard materials(NA12878) 
        This requires sample bed, vcf files 
        And get True positive, False negative and False positive  
    """
    bed, vcf = intersect_bed(input_bed, output, prefix)
    sample_vcf = intersect_vcf_bed(input_vcf, input_bed, output, prefix)
    compare(vcf, sample_vcf, output, tn, prefix)

@cli.command()
@click.option('-b', '--input_bed', required=True, type=click.Path(exists=True, dir_okay=False), help='bed file including region information you want to subset')
@click.option('-o', '--output', required=True, type=click.Path(exists=False, dir_okay=True), help='output directory')
@click.option('-p', '--prefix', required=True,  type=click.STRING, help='Output prefix string. e.g. output file: prefix.intersectedNA12878.bed')
def intersect_bed(input_bed, output, prefix):
    """ Intersect sample bed with NA12878 bed of Standard Materials(HG001) 
        And get intersected NA12878 VCF with output bed
    """
    out_bed = os.path.join(output, f'{prefix}.intersectedNA12878.bed')
    cmd = (f'bedtools intersect -a {NA12878_BED} '
           f'-b {input_bed} '
           f'> {out_bed}') 
    print(f'[cmd] {cmd}')       
    subprocess.run(cmd, shell = True)

    bed = pd.read_csv(out_bed, sep = '\t', header = None, usecols=[0,1,2])
    bed.columns = ['chr', 'start', 'end']
    bed['length'] = bed['end'] - bed['start']
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
@click.option('-v', '--input_vcf', required=True, type=click.Path(exists=True, dir_okay=False), help='vcf file you want to subset')
@click.option('-o', '--output', required=True, type=click.Path(exists=False, dir_okay=True), help='output directory')
def subset_chrom(input_vcf, output):
    """Input vcf is devided by chromosome automatically\n"""
    basename = os.path.basename(input_vcf)
    for chrom in range(1, 23):
        print(f"Subset chromosome {chrom}")
        cmd = (f'grep -w "^chr{chrom}" '
               f'{input_vcf} > {output}/{basename}.chr{chrom}')
        subprocess.run(cmd, shell = True)
    print("Subset chromosome X")
    subprocess.run(f'grep -w "^chrX" {input_vcf} > {output}/{basename}.chrX', shell = True)

def compare(standard, sample, output, tn, prefix):
    """ Performance evaluation \n
        This requires only correct answer vcf and sample vcf with no header\n
        The VCF files must have whole 9 cloummns\n
        Finally get True positive, False negative and False positive vcf of all variants \n
    """
    print('[INFO] Compare vcf Start')
    standard = pd.read_csv(standard, usecols = [0,1,2,3,4,9], dtype = {0: object}, sep = '\t', header = None)
    standard = standard.drop_duplicates()
    standard.columns = ['chr', 'pos', 'rsid', 'ref', 'alt', 'genotype']
    standard['genotype_tmp'] = standard['genotype'].str.slice(start = 0, stop = 3)
    for i in standard.index:
        if (standard.loc[i, 'genotype_tmp'] == '0|1') or (standard.loc[i, 'genotype_tmp'] == '1|0') or (standard.loc[i, 'genotype_tmp'] == '0/1') or (standard.loc[i, 'genotype_tmp'] == '1/0'):
            standard.loc[i, 'gt'] = ''.join(sorted(standard.loc[i, 'ref'] + standard.loc[i, 'alt']))
        elif (standard.loc[i, 'genotype_tmp'] == '1|2') or (standard.loc[i, 'genotype_tmp'] == '2|1') or (standard.loc[i, 'genotype_tmp'] == '1/2') or (standard.loc[i, 'genotype_tmp'] == '2/1'):
            standard.loc[i, 'gt'] = ''.join(sorted(standard.loc[i, 'alt'].split(',')[0] + standard.loc[i, 'alt'].split(',')[1]))
        elif (standard.loc[i, 'genotype_tmp'] == '1|1') or (standard.loc[i, 'genotype_tmp'] == '1/1') :
            standard.loc[i, 'gt'] = standard.loc[i, 'alt'] + standard.loc[i, 'alt']
        elif (standard.loc[i, 'genotype_tmp'] == '0|0') or (standard.loc[i, 'genotype_tmp'] == '0/0') :
            standard.loc[i, 'gt'] = standard.loc[i, 'ref'] + standard.loc[i, 'ref']   
        else:
            standard.loc[i, 'gt'] = 'NaN'


    sample = pd.read_csv(sample, usecols = [0,1,2,3,4,9], dtype = {0: object}, sep = '\t', header = None)
    sample = sample.drop_duplicates()
    sample.columns = ['chr', 'pos', 'id', 'ref', 'alt', 'genotype']
    sample['genotype_tmp'] = sample['genotype'].str.slice(start = 0, stop = 3)
    #sample['GQ'] = sample['genotype'].str.slice(start = 4, stop = 9).astype('float64')
    sample['GQ'] = sample.genotype.str.split(':').str[3]
    for i in sample.index:
        if (sample.loc[i, 'genotype_tmp'] == '0/1') or (sample.loc[i, 'genotype_tmp'] == '1/0') or (sample.loc[i, 'genotype_tmp'] == '0|1'):
            sample.loc[i, 'gt'] = ''.join(sorted(sample.loc[i, 'ref'] + sample.loc[i, 'alt']))
        elif (sample.loc[i, 'genotype_tmp'] == '1/1') or (sample.loc[i, 'genotype_tmp'] == '1|1'):
            sample.loc[i, 'gt'] = sample.loc[i, 'alt'] + sample.loc[i, 'alt']
        elif (sample.loc[i, 'genotype_tmp'] == '0/0') or (sample.loc[i, 'genotype_tmp'] == '0|0') :
            sample.loc[i, 'gt'] = sample.loc[i, 'ref'] + sample.loc[i, 'ref']
        elif (sample.loc[i, 'genotype_tmp'] == '1/2') or (sample.loc[i, 'genotype_tmp'] == '2/1'):
            sample.loc[i, 'gt'] = ''.join(sorted(sample.loc[i, 'alt'].split(',')[0] + sample.loc[i, 'alt'].split(',')[1]))
        elif (sample.loc[i, 'genotype_tmp'] == '0/2') or (sample.loc[i, 'genotype_tmp'] == '0|2'):
            sample.loc[i, 'gt'] = ''.join(sorted(sample.loc[i, 'ref'] + sample.loc[i, 'alt'].split(',')[1]))
        elif (sample.loc[i, 'genotype_tmp'] == '2/2') or (sample.loc[i, 'genotype_tmp'] == '2|2'):
            sample.loc[i, 'gt'] = ''.join(sorted(sample.loc[i, 'alt'].split(',')[1] + sample.loc[i, 'alt'].split(',')[1]))
        elif (sample.loc[i, 'genotype_tmp'] == '0/3'):
            sample.loc[i, 'gt'] = ''.join(sorted(sample.loc[i, 'ref'] + sample.loc[i, 'alt'].split(',')[2]))
        elif (sample.loc[i, 'genotype_tmp'] == '1/3'):
            sample.loc[i, 'gt'] = ''.join(sorted(sample.loc[i, 'alt'].split(',')[0] + sample.loc[i, 'alt'].split(',')[2]))
        elif (sample.loc[i, 'genotype_tmp'] == '2/3'):
            sample.loc[i, 'gt'] = ''.join(sorted(sample.loc[i, 'alt'].split(',')[1] + sample.loc[i, 'alt'].split(',')[2]))
        elif (sample.loc[i, 'genotype_tmp'] == '3/3'):
            sample.loc[i, 'gt'] = ''.join(sorted(sample.loc[i, 'alt'].split(',')[2] + sample.loc[i, 'alt'].split(',')[2]))
        else:
            sample.loc[i, 'gt'] = 'NaN'

    TP = pd.merge(standard, sample, how = 'inner', on = ['chr', 'pos', 'gt'], suffixes=("_std", "_sp"))

    FP_tmp = pd.concat([sample, TP]).drop_duplicates(subset = ['chr', 'pos', 'gt'], keep = False)
    # FP가 없는 경우 헤더까지 없어지는것 방지 
    if len(FP_tmp) > 0:
        FP_tmp.dropna(axis = 1, how = 'all', inplace = True)
    #FP = FP_tmp[(FP_tmp['genotype_tmp'] == '0/1') | (FP_tmp['genotype_tmp'] == '1/0') | (FP_tmp['genotype_tmp'] == '1/1')]
    FP = FP_tmp.loc[~FP_tmp['genotype_tmp'].isin(['./.', '0/0'])]
    FP = FP[['chr', 'pos', 'ref', 'alt', 'id', 'genotype', 'genotype_tmp', 'gt', 'GQ']]
    
    FN_tmp = pd.concat([standard, TP]).drop_duplicates(subset = ['chr', 'pos', 'gt'], keep = False)
    # FN가 없는 경우 헤더까지 없어지는것 방지 
    if len(FN_tmp) > 0:
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
        FN.dropna(axis = 1, how = 'all', inplace = True)
        FN = FN[['chr', 'pos', 'ref', 'alt', 'rsid', 'genotype', 'genotype_tmp', 'gt']]
        FN.to_csv(False_neg, sep = '\t' , index = False)
    if len(uncalled) != 0:
        uncalled.to_csv(Uncalled, sep = '\t' , index = False)

    tp = len(TP)
    fp = len(FP)
    fn = len(FN)
    TN = tn - tp -fp - fn

    print(f'TP = {len(TP)}, FP = {len(FP)}, FN = {len(FN)}')
    print('[INFO] Compare vcf End')
    print('[INFO] Performance Evaluation Start')
   
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
    print('[INFO] Performance Evaluation End')     

def delete_files(files, root):
    for file in files:
        os.remove(os.path.join(root, file))

def compare_large(standard, samplefile, output, prefix):
    """ Performance evaluation \n
        This requires only correct answer vcf and sample vcf with no header\n
        The VCF files must have whole 9 cloummns\n
        Finally get True positive, False negative and False positive vcf of all variants \n
    """
    print(f'[INFO] START compare {samplefile}')
    sample = pd.read_csv(samplefile, usecols = [0,1,2,3,4,9], dtype = {0: object}, sep = '\t', header = None)
    sample.columns = ['chr', 'pos', 'id', 'ref', 'alt', 'genotype']
    #sample['genotype_tmp'] = sample['genotype'].str.slice(start = 0, stop = 3)
    #sample['GQ'] = sample.genotype.str.split(':').str[3]
    revised_sample = pd.DataFrame()
    for i in range(0, len(sample), 5000):
        start = i
        end = i+5000
        df = sample.iloc[start:end, :]
        df['genotype_tmp'] = df['genotype'].str.slice(start = 0, stop = 3)
        df['GQ'] = df.genotype.str.split(':').str[3]
        for i in df.index:
            if (df.loc[i, 'genotype_tmp'] == '0/1') or (df.loc[i, 'genotype_tmp'] == '1/0') or (df.loc[i, 'genotype_tmp'] == '0|1'):
                df.loc[i, 'gt'] = ''.join(sorted(df.loc[i, 'ref'] + df.loc[i, 'alt']))
            elif (df.loc[i, 'genotype_tmp'] == '1/1') or (df.loc[i, 'genotype_tmp'] == '1|1'):
                df.loc[i, 'gt'] = ''.join(sorted(df.loc[i, 'alt'] + df.loc[i, 'alt']))
            elif (df.loc[i, 'genotype_tmp'] == '0/0') or (df.loc[i, 'genotype_tmp'] == '0|0') :
                df.loc[i, 'gt'] = ''.join(sorted(df.loc[i, 'ref'] + df.loc[i, 'ref']))
            elif (df.loc[i, 'genotype_tmp'] == '1/2') or (df.loc[i, 'genotype_tmp'] == '2/1'):
                df.loc[i, 'gt'] = ''.join(sorted(df.loc[i, 'alt'].split(',')[0] + df.loc[i, 'alt'].split(',')[1]))
            elif (df.loc[i, 'genotype_tmp'] == '0/2') or (df.loc[i, 'genotype_tmp'] == '0|2'):
                df.loc[i, 'gt'] = ''.join(sorted(df.loc[i, 'ref'] + df.loc[i, 'alt'].split(',')[1]))
            elif (df.loc[i, 'genotype_tmp'] == '2/2') or (df.loc[i, 'genotype_tmp'] == '2|2'):
                df.loc[i, 'gt'] = ''.join(sorted(df.loc[i, 'alt'].split(',')[1] + df.loc[i, 'alt'].split(',')[1]))
            elif (df.loc[i, 'genotype_tmp'] == '0/3'):
                df.loc[i, 'gt'] = ''.join(sorted(df.loc[i, 'ref'] + df.loc[i, 'alt'].split(',')[2]))
            elif (df.loc[i, 'genotype_tmp'] == '1/3'):
                df.loc[i, 'gt'] = ''.join(sorted(df.loc[i, 'alt'].split(',')[0] + df.loc[i, 'alt'].split(',')[2]))
            elif (df.loc[i, 'genotype_tmp'] == '2/3'):
                df.loc[i, 'gt'] = ''.join(sorted(df.loc[i, 'alt'].split(',')[1] + df.loc[i, 'alt'].split(',')[2]))
            elif (df.loc[i, 'genotype_tmp'] == '3/3'):
                df.loc[i, 'gt'] = ''.join(sorted(df.loc[i, 'alt'].split(',')[2] + df.loc[i, 'alt'].split(',')[2]))
            else:
                df.loc[i, 'gt'] = 'NaN'
        revised_sample = pd.concat([revised_sample, df], ignore_index = True).drop_duplicates()
   
    standard = pd.read_csv(standard, usecols = [0,1,2,3,4,9], dtype = {0: object}, sep = '\t', header = None)
    standard.columns = ['chr', 'pos', 'rsid', 'ref', 'alt', 'genotype']
    TP_merge = pd.DataFrame()
    revised_standard = pd.DataFrame()
    for i in range(0, len(standard), 5000):
        start = i
        end = i+5000
        df = standard.iloc[start:end, :]
        df['genotype_tmp'] = df['genotype'].str.slice(start = 0, stop = 3)
        for i in df.index:
            if (df.loc[i, 'genotype_tmp'] == '0|1') or (df.loc[i, 'genotype_tmp'] == '1|0') or (df.loc[i, 'genotype_tmp'] == '0/1') or (df.loc[i, 'genotype_tmp'] == '1/0'):
                df.loc[i, 'gt'] = ''.join(sorted(df.loc[i, 'ref'] + df.loc[i, 'alt']))
            elif (df.loc[i, 'genotype_tmp'] == '1|2') or (df.loc[i, 'genotype_tmp'] == '2|1') or (df.loc[i, 'genotype_tmp'] == '1/2') or (df.loc[i, 'genotype_tmp'] == '2/1'):
                df.loc[i, 'gt'] = ''.join(sorted(df.loc[i, 'alt'].split(',')[0] + df.loc[i, 'alt'].split(',')[1]))
            elif (df.loc[i, 'genotype_tmp'] == '1|1') or (df.loc[i, 'genotype_tmp'] == '1/1') :
                df.loc[i, 'gt'] = ''.join(sorted(df.loc[i, 'alt'] + df.loc[i, 'alt']))
            elif (df.loc[i, 'genotype_tmp'] == '0|0') or (df.loc[i, 'genotype_tmp'] == '0/0') :
                df.loc[i, 'gt'] = ''.join(sorted(df.loc[i, 'ref'] + df.loc[i, 'ref']))
            else:
                df.loc[i, 'gt'] = 'NaN'
        
        TP = pd.merge(df, revised_sample, how = 'inner', on = ['chr', 'pos', 'gt'], suffixes=("_std", "_sp"))
        TP_merge = pd.concat([TP_merge, TP], ignore_index = True).drop_duplicates()
        revised_standard = pd.concat([revised_standard, df], ignore_index = True).drop_duplicates() 

    FP_tmp = pd.concat([revised_sample, TP_merge]).drop_duplicates(subset = ['chr', 'pos', 'gt'], keep = False)
    # FP가 없는 경우 헤더까지 없어지는것 방지 
    if len(FP_tmp) > 0:
        FP_tmp.dropna(axis = 1, how = 'all', inplace = True)
    #FP = FP_tmp[(FP_tmp['genotype_tmp'] == '0/1') | (FP_tmp['genotype_tmp'] == '1/0') | (FP_tmp['genotype_tmp'] == '1/1')]
    FP = FP_tmp.loc[~FP_tmp['genotype_tmp'].isin(['./.', '0/0'])]
    FP = FP[['chr', 'pos', 'ref', 'alt', 'id', 'genotype', 'genotype_tmp', 'gt', 'GQ']]
    
    FN_tmp = pd.concat([revised_standard, TP_merge]).drop_duplicates(subset = ['chr', 'pos', 'gt'], keep = False)
    
    # FN가 없는 경우 헤더까지 없어지는것 방지 
    if len(FN_tmp) > 0:
        FN_tmp.dropna(axis = 1, how = 'all', inplace = True)
    
    #FP와 중복되는 SNP삭제 
    if len(FP) != 0 and len(FN_tmp) != 0:
        #tmp = pd.merge(FN_tmp, FP, on = ['chr', 'pos', 'gt'], how = 'inner')
        tmp = pd.merge(FN_tmp, FP, on = ['chr', 'pos'], how = 'inner')
        FN = pd.concat([FN_tmp, tmp]).drop_duplicates(subset = ['chr', 'pos'], keep = False)
    else:
        FN = FN_tmp

    True_pos = os.path.join(output, f'{prefix}.TP')
    False_pos = os.path.join(output,f'{prefix}.FP')
    False_neg = os.path.join(output,f'{prefix}.FN')
    
    #결과 저장 
    TP_merge.to_csv(True_pos, sep = '\t', index = False)
    if len(FP) != 0:
        FP.to_csv(False_pos, sep = '\t', index = False)
    if len(FN) != 0:
        FN.dropna(axis = 1, how = 'all', inplace = True)
        FN = FN[['chr', 'pos', 'ref', 'alt', 'rsid', 'genotype', 'genotype_tmp', 'gt']]  
        FN.to_csv(False_neg, sep = '\t', index = False)
  
    print(f'[INFO] END compare {samplefile}')
     
@cli.command()
@click.option('-r', '--standard', required=True, type=click.Path(exists=False, dir_okay=False), help='Referfence(NA12878) vcf W/O header')
@click.option('-i', '--sample', required=True, type=click.Path(exists=False, dir_okay=False), help='input sample vcf')
@click.option('-n', '--tn', required=True, type=click.INT, help='variants number of standard bed')
@click.option('-o', '--outdir', required=True, type=click.Path(exists=False, dir_okay=True), help='output directory')
@click.option('-p', '--prefix', required=True, type=click.STRING, help='Output prefix string, (preferably sample ID)')
@click.option('-t', '--threads', required=False, type=click.INT, default = 1, help='Number of threads to use')
@click.option('--chrom', is_flag=True, default=False, help='A boolean flag on whether to compare vcf by chromosome. In this case, stardard and sample option is used as the basename of vcf')
def evaluate(standard, sample, outdir, tn, prefix, threads, chrom):
    """ Performance evaluation using two vcf files\n
        This requires only correct answer vcf and sample vcf\n
        If chrom, Compare vcf files by chromosome and combine the results\n
        And get True positive, False negative and False positive vcf of all variants\n
    """
    if chrom:
        # compare vcf
        pool = Pool(processes = threads)
        std = list()
        samp = list()
        out = list()
        pref =list()
    
        for i in range(1, 23):
            std.append(f'{standard}.chr{i}')
            samp.append(f'{sample}.chr{i}')
            out.append(outdir)
            pref.append(f'{prefix}.chr{i}')   
   
        std.append(f'{standard}.chrX')
        samp.append(f'{sample}.chrX')
        out.append(outdir)
        pref.append(f'{prefix}.chrX')
    
        pool.starmap(compare_large, zip(std, samp, out, pref))
        pool.close()
        pool.join()
    
        # calculate 
        root, _, files = next(os.walk(outdir))
        TP = [f for f in files if f.endswith('TP')]
        FP = [f for f in files if f.endswith('FP')]
        FN = [f for f in files if f.endswith('FN')]
        
        num_tp = 0
        df_tp = pd.DataFrame()
        for i in TP:
            df = pd.read_csv(os.path.join(outdir, i), sep = '\t')  
            num_tp += len(df)
            df_tp = pd.concat([df_tp, df], ignore_index = True)
        df_tp.to_csv(os.path.join(outdir, f'{prefix}.TP'), sep = '\t', index = False)
    
        num_fp = 0
        df_fp = pd.DataFrame()
        for i in FP:
            df = pd.read_csv(os.path.join(outdir, i), sep = '\t')
            num_fp += len(df) 
            df_fp = pd.concat([df_fp, df], ignore_index = True)
        df_fp.to_csv(os.path.join(outdir, f'{prefix}.FP'), sep = '\t', index = False)
        
        num_fn = 0
        df_fn = pd.DataFrame()
        for i in FN:
            df = pd.read_csv(os.path.join(outdir, i), sep = '\t')
            num_fn += len(df)
            df_fn = pd.concat([df_fn, df], ignore_index = True)
        df_fn.to_csv(os.path.join(outdir, f'{prefix}.FN'), sep = '\t', index = False)

        TN = tn - num_tp - num_fp - num_fn
    
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
        
        Sensitivity = (num_tp/(num_tp + num_fn)*100) 
        Specificity = (TN/(num_fp + TN)*100)
        ppv = (num_tp/(num_tp + num_fp)*100)
        npv = (TN/(num_fn + TN)*100)
        
        _dict['sampleID'].append(prefix)
        _dict['TruePositive'].append(num_tp)
        _dict['FalsePositive'].append(num_fp)
        _dict['FalseNegative'].append(num_fn)
        _dict['TrueNegative'].append(TN)
        _dict['Sensitivity'].append(Sensitivity)
        _dict['Specificity'].append(Specificity)
        _dict['PPV'].append(ppv) 
        _dict['NPV'].append(npv)
    
        result = pd.DataFrame(_dict)
        output = os.path.join(outdir, f'{prefix}_evaluation.txt')
        result.to_csv(output, sep = '\t', index = False)
        print(result)

    else:
        compare(standard, sample, outdir, tn, prefix)

if __name__ == '__main__':
    cli.main()
