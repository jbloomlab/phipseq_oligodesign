"""Script for making oligos for Zika, Dengue, etc. phage display assay."""

import Bio.Seq
import Bio.SeqUtils.CodonUsageIndices
import Bio.Alphabet.IUPAC

from Bio.Seq import Seq
from operator import itemgetter

viruses = ['CHIKV_KPA15_NS', 'CHIKV_KPA15_S', 'DENV1_Myanmar1996', 'DENV1_WesternPacific', 'DENV2_16803', 'DENV2_BurkinaFaso1983',
'DENV2_NewGuineaC', 'DENV3_Mozambique1985', 'DENV4_P73_1120', 'HIV_C5_GroupM', 'HIV_gp41_GroupM', 
'HIV_MPER_GroupM', 'HIV_V3_GroupM', 'HIV_VirScanOligo', 'HIV_Env_BF520.W14.C2', 'HIV_Env_BG505.W6.C2', 'HIV_Env_CladeA1', 'HIV_Env_CladeA2',
'HIV_Env_CladeB', 'HIV_Env_CladeC', 'HIV_Env_CladeD', 'HIV_Env_Q23', 'HIV_Env_Q461.d1', 'HIV_Env_QA013.70I.ENV.H1', 'HIV_Env_QA013.385M.ENV.R3',
'HIV_Env_QC406.F3', 'HIV_EnvV1V5_QB850_52M.Ev1v5.A', 'HIV_EnvV1V5_QB850.73M.Ev1v5.Q', 'HIV_Gag_CladeA1', 'HIV_Gag_CladeA2', 'HIV_Gag_CladeB',
'HIV_Gag_CladeC', 'HIV_Gag_CladeD', 'HIV_Gag_Q23', 'HIV_Nef_CladeA', 'HIV_Nef_CladeB', 'HIV_Nef_CladeC', 'HIV_Nef_CladeD', 'HIV_Nef_Q23',
'HIV_Pol_CladeA1', 'HIV_Pol_CladeA2', 'HIV_Pol_CladeB', 'HIV_Pol_CladeC', 'HIV_Pol_CladeD', 'HIV_Pol_Q23', 'HIV_Rev_CladeA1', 'HIV_Rev_CladeA2',
'HIV_Rev_CladeB', 'HIV_Rev_CladeC', 'HIV_Rev_CladeD', 'HIV_Rev_Q23', 'HIV_Tat_CladeA1', 'HIV_Tat_CladeA2', 'HIV_Tat_CladeB', 'HIV_Tat_CladeC',
'HIV_Tat_CladeD', 'HIV_Tat_Q23', 'HIV_Vif_CladeA', 'HIV_Vif_CladeB', 'HIV_Vif_CladeC', 'HIV_Vif_CladeD', 'HIV_vpr_CladeA', 'HIV_vpr_CladeB',
'HIV_vpr_CladeC', 'HIV_vpr_CladeD', 'HIV_vpu_CladeA', 'HIV_vpu_CladeB', 'HIV_vpu_CladeC', 'HIV_vpu_CladeD',
'Malaria_HRP2', 'ONNV_S_SG650', 'ONNV_NS_SG650', 'RhinovirusB_VirScan', 'RVFV_Lseg_Kenya98', 'RVFV_Mseg_Kenya98', 'RVFV_Sseg_NC_Kenya98', 
'RVFV_Sseg_NS_Kenya98', 'WNV_NY99', 'WNV_WN02', 'YFV_Uganda1948a', 'ZIKV_ArD128000', 'ZIKV_ArD157995', 'ZIKV_HPF2013',
'ZIKV_MR766_Peptide451-472', 'ZIKV_MR766_Peptide1423-1448', 'ZIKV_MR766', 'CHIKV_TR206-H804187_NS', 'CHIKV_TR206-H804187_S',
'DENV1_BR-BID-V2401-2008', 'DENV2_BR-BID-V3653-2008', 'DENV3_BR-BID-V2403-2008', 'DENV4_H781363', 'DENV4_SJRP-850-2013',
'YFV_BeH655417', 'ZIKV_PRVABC59', 'JEV_Laos2009', 'TBEV_MDJ03', 'HIV_Vif_Q23', 'HIV_Vpr_Q23', 'HIV_Vpu_Q23', 'HIV_Env_QB850.72p.C14_A1',
'HIV_Env_QB850.632p.B10'] #List of file names for creating libraries

oligo_length = 117
tile = 60

avoidmotifs = ['GAATTC', 'AAGCTT']

adaptor5 = 'aggaattctacgctgagt' 
adaptor3 = 'tgatagcaagcttgcc' 

cai = Bio.SeqUtils.CodonUsageIndices.SharpEcoliIndex

rt_table = {
	'F': ['TTT', 'TTC'],
	'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
	'I': ['ATT', 'ATC', 'ATA'],
	'M': ['ATG'],
	'V': ['GTT', 'GTC', 'GTA', 'GTG'],
	'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
	'P': ['CCT', 'CCC', 'CCA', 'CCG'],
	'T': ['ACT', 'ACC', 'ACA', 'ACG'],
	'A': ['GCT', 'GCC', 'GCA', 'GCG'],
	'Y': ['TAT', 'TAC'],
	'*': ['TAA', 'TAG', 'TGA'],
	'H': ['CAT', 'CAC'],
	'Q': ['CAA', 'CAG'],
	'N': ['AAT', 'AAC'],
	'K': ['AAA', 'AAG'],
	'D': ['GAT', 'GAC'],
	'E': ['GAA', 'GAG'],
	'C': ['TGT', 'TGC'],
	'W': ['TGG'],
	'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
	'G': ['GGT', 'GGC', 'GGA', 'GGG']
}

translation_table = {}
codons = []
amino_acids = []
for _nt1 in Bio.Alphabet.IUPAC.IUPACUnambiguousDNA.letters:
    for _nt2 in Bio.Alphabet.IUPAC.IUPACUnambiguousDNA.letters:
        for _nt3 in Bio.Alphabet.IUPAC.IUPACUnambiguousDNA.letters:
            codons.append('{0}{1}{2}'.format(_nt1, _nt2, _nt3))
            for codon in codons:
                translation_table[codon] = str(Bio.Seq.Seq(codon).translate())
                if translation_table[codon] not in amino_acids:
                    amino_acids.append(translation_table[codon])

score_dict = {} #Dictionary of all the synonymous codons encoding an amino acid and each codon's associated cai score. Does not include stop codons.
for amino_acid in amino_acids:
	if amino_acid != '*':
		scores = []
		for codon in rt_table[amino_acid]:
			scores.append((codon, cai[codon]))
		score_dict[amino_acid] = scores


def translate(oligo):
    """Script for translating DNA more quickly than Bio.Seq"""
    assert len(oligo) % 3 == 0, "Not translatable due to length"
    aa_read = ''.join([translation_table[oligo[i : i + 3]] for i in range(0, len(oligo), 3)])
    return aa_read


def reverse_translate(protein_sequence): 
	""" Reverse_translate viral protein sequence into DNA sequence.

	Codon preferences based on cai for E. coli."""
	#header = False
	revtranslated_seq = []
	with open(protein_sequence) as f:
		for line in f:
			#if '>' in line:
				#header = True # I want to be able to set this up so can input sequence with fasta header
			seq = line.strip()
			for aa in seq: #Go through seq aa by aa
				assert aa in rt_table, 'Ambiguous Sequence. Check {0} for accuracy of protein sequence.'.format(protein_sequence)
				max_codon = max(score_dict[aa], key=lambda item:item[1])[0]
				revtranslated_seq.append(max_codon)
	rt_seq = ''.join(revtranslated_seq)
	return rt_seq


def subsequence(genome):
	""" Subsequence viral genome into oligos of length *oligo_length* tiled by *tile*."""
	subsequences = []
	i = 0
	while i < (len(genome) - oligo_length):
		subsequences.append(genome[i : i + oligo_length])
		i = i + tile
	subsequences.append(genome[len(genome)-oligo_length: len(genome)]) #The final oligo goes from the end of the sequence backwards oligo length.
	return subsequences


def remove_rsites(oligos, codon_scores_by_aa):
    clean_oligos = []
    rsites_count = 0
    replace = False
    for motif in avoidmotifs:
    	for n in range(len(oligos)):
    		if motif in oligos[n]:
    			aa_withrs = translate(oligos[n])
    			replace = True
    			while replace:
        			for i in range(len(oligos[n]) - len(motif) + 1):
        				seq = oligos[n][i : i + len(motif)] # sequence starting at i
        				if seq == motif:
        					done = False #need to replace a codon in this seq
        					rsites_count += 1
        					if i % 3 == 0: #The restriction site is in frame.
        						for x in codon_scores_by_aa: #For each amino acid...
        							for y in codon_scores_by_aa[x]: #...And for each (codon, score) tuple that encodes that aa...
	       								if not done: #If the codon at this location has not already been replaced.
		       								if oligos[n][i:i+3] in y: #If the codon we want to replace is in the score tuple y... 
	        									l = codon_scores_by_aa[x] #...Turn all synonymous codons and their scores into a list
        										l.sort(key=itemgetter(1),reverse=True)
        										assert oligos[n][i:i+3] == l[0][0], 'The codon we are replacing is not the highest scoring. We are replacing: {0}'.format(oligos[n][i:i+3])
        										new_codon = l[1][0] # The new codon is the codon that has the second highset score (so is second in the sorted list)
        										clean_oligo = oligos[n][:i]+new_codon+oligos[n][i+3:]
        										oligos[n] = clean_oligo
        										done = True #The codon has been replaced, don't look at this site anymore. 
        					else: #The restriction site is not in frame. Go through the same steps as if it were. 
        						for x in codon_scores_by_aa:
        							for y in codon_scores_by_aa[x]:
        								if not done:
        									if oligos[n][i-(i%3):i-(i%3)+3] in y: #Replace the first codon that contains part of the restriciton site.
       											l = codon_scores_by_aa[x]
       											l.sort(key=itemgetter(1), reverse=True)
       											assert oligos[n][i-(i%3):i-(i%3)+3] == l[0][0], 'The codon we are replacing is not the highest scoring. We are replacing: {0}'.format(oligos[n][i-(i%3):i-(i%3)+3])
       											new_codon = l[1][0]
       											clean_oligo = oligos[n][:i-(1%3)]+new_codon+oligos[n][i-(i%3)+3:]
       											oligos[n] = clean_oligo
       											done = True #The codon has been replaced, don't look at other possible synonymous codons.
       				if motif not in oligos[n]: #Make sure motif not in oligo after going through replacements
       					replace = False
       			assert len(oligos[n]) == oligo_length, 'Oligo lengths not maintained while removing restriction sites.'  		
    			assert aa_withrs == translate(oligos[n]), 'The amino acid sequence has been altered by removing restriction site.'
    print('{0} restriction sites were removed with synonymous substitution.'.format(rsites_count))
    #print oligos
    return oligos


def add_adaptors(cleaned_oligos):
	final_oligos = []
	for oligo in cleaned_oligos:
		final_oligo = adaptor5 + oligo + adaptor3
		final_oligos.append(final_oligo)
	return final_oligos

def main():
	for virus in viruses:
		print('Designing oligos for {0}.'.format(virus))
		prot_seq = './sequences/%s.fasta' % virus
		rt_seq = reverse_translate(prot_seq)
		with open(prot_seq) as f:
			full_seq = ''.join(line.strip() for line in f)
			assert full_seq == translate(rt_seq), "reverse translation did not work" 
		subsequences = subsequence(rt_seq)
		cleaned_oligos = remove_rsites(subsequences, score_dict)
		final_oligos = add_adaptors(cleaned_oligos)
		final_seqs_file = './oligos/%s_final_oligos.txt' % virus
		with open(final_seqs_file, 'w') as f:
			for oligo in range(len(final_oligos)):
				f.write(','.join((virus, str(oligo*tile+1) + 'to' + str(oligo*tile+1+oligo_length), final_oligos[oligo] + '\n')))


if __name__ == '__main__':
	main()