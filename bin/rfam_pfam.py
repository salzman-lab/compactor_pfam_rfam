import numpy as np
import pandas as pd
import subprocess
import os
import sys

nom = sys.argv[1]
dt = pd.read_csv(sys.argv[2],engine='python',sep='\t')
dt = dt[dt['exact_support']>0].reset_index()       
                            
rfam = 'bin/rfam.bash'
pfam = 'bin/pfam.bash'

fasta = open(nom+'.fasta','a')
for i in range(len(dt)):
    fasta.write('>seq_'+str(i)+'\n'+dt['compactor'][i]+'\n')
fasta.close()

pfam_process = subprocess.Popen([pfam+' '+nom+'.fasta'],shell=True,stdout=subprocess.PIPE)
rfam_process = subprocess.Popen([rfam+' '+nom+'.fasta'],shell=True,stdout=subprocess.PIPE)               

pfam_process.wait()
rfam_process.wait()

pfam_dt = pd.read_csv(nom+'.fasta_PFAM.tblout',header=None)
fields = ['target_name', 'accession', 'query_name', 'accession2', 'full_seq_evalue', 'full_seq_score', 'full_seq_bias', 'best_1_domain_evalue', 'best_1_domain_score', 'best_1_domain_bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description_of_target']
pfam_list_compactor = [i for i in list(pfam_dt.iloc[:,0]) if '#' not in i]

compactor_ok = [ i.split(' ') for i in pfam_list_compactor]
compactor_okok = []
for i in compactor_ok: 
    compactor_okok.append([j for j in i if j])
compactor_pfam_structured = pd.DataFrame(compactor_okok, columns = fields)
compactor_pfam_structured = compactor_pfam_structured.add_prefix('Pfam_') 

ntfast = pd.read_csv(nom+'.fasta',header=None)
seqid=[i[1:] for i in ntfast[0] if '>' in i]
ntseq=[i for i in ntfast[0] if '>' not in i]
ntfast=pd.DataFrame({'compactor':ntseq,'Pfam_target_name':seqid})
ntfast=ntfast.drop_duplicates()
    
if compactor_pfam_structured.shape[0] != 0:
    aafast = pd.read_csv(nom+'.fastaCOMPACTORS.aa',header=None)
    targe = [i[1:] for i in aafast[0] if '>' in i]
    aaseq = [i for i in aafast[0] if '>' not in i]

    aafast=pd.DataFrame({'Pfam_target_name':targe,'protein_seq':aaseq})
    aafast['Pfam_target_name'] = [i.split(' ')[0] for i in aafast['Pfam_target_name']]
    compactor_pfam_structured=compactor_pfam_structured.merge(aafast,how='left')
    compactor_pfam_structured['translation_frame'] = [i.split('frame=')[1] for i in compactor_pfam_structured['Pfam_target_name']]
    compactor_pfam_structured['Pfam_target_name'] = [i.split('_frame')[0] for i in compactor_pfam_structured['Pfam_target_name']]
    compactor_pfam_structured=compactor_pfam_structured.sort_values(by='Pfam_full_seq_evalue',ascending=True).drop_duplicates(subset='Pfam_target_name')
    compactor_pfam_structured=compactor_pfam_structured.merge(ntfast,how='left')


rfam_dt = pd.read_csv('RFAM.tblout',header=None)
fields = ['idx','target_name', 'accession', 'query_name', 'accession2', 'clan_name', 'mdl', 'mdl_from','mdl_to','seq_from', 'seq_to','strand','trunc', 'pass', 'gc', 'bias', 'score', 'E-value', 'inc', 'olp', 'anyidx', 'afrct1', 'afrct2', 'winidx', 'wrfct1', 'wfrct2', 'description_of_target']
pfam_list_compactor = [i for i in list(rfam_dt.iloc[:,0]) if '#' not in i]

compactor_ok = [ i.split(' ') for i in pfam_list_compactor]
compactor_okok = []
for i in compactor_ok: 
    compactor_okok.append([j for j in i if j])
rfam_dt = pd.DataFrame(compactor_okok, columns = fields)
rfam_dt = rfam_dt.add_prefix('Rfam_') 

if rfam_dt.shape[0] != 0: 
    rfam_dt = rfam_dt.sort_values(by='Rfam_E-value',ascending=True).drop_duplicates(subset='Rfam_query_name')
    ntfast = ntfast.rename(columns={'Pfam_target_name':'Rfam_query_name'})
    rfam_dt = rfam_dt.merge(ntfast,how='left')
    compactor_pfam_structured = compactor_pfam_structured.merge(rfam_dt,how='outer')
    
dt = dt.merge(compactor_pfam_structured,how='left')
dt.to_csv(nom+'_compactor_Pfam_Rfam.tsv',sep='\t',index=None)