import sys
import primer3
import re

#Version: 1.0.0
# This script is designed to create custom ASQ (Allele-Specific Quantification) markers by designing primers
# specific to single nucleotide polymorphisms (SNPs) in a given DNA sequence. The script includes functions to:
# - Introduce mismatches at specific positions in the primers to enhance specificity
# - Generate reverse complements of sequences
# - Design forward and reverse primers using the Primer3 library
# - Optimize PCR conditions for the designed primers

degenerate_mapping = {
    'A': ['A'], 
    'C': ['C'], 
    'G': ['G'], 
    'T': ['T'], 
    'R': ['A', 'G'], 
    'Y': ['C', 'T'], 
    'S': ['G', 'C'], 
    'W': ['A', 'T'], 
    'K': ['G', 'T'], 
    'M': ['A', 'C'], 
    'B': ['C', 'G', 'T'], 
    'D': ['A', 'G', 'T'], 
    'H': ['A', 'C', 'T'], 
    'V': ['A', 'C', 'G'], 
    'N': ['A', 'C', 'G', 'T']
}

def introduce_mismatch(base):
    """
    Introduce a mismatch at the antepenultimate position based on the given rules.
    
    Args:
    base (str): The original base at the antepenultimate position.
    
    Returns:
    str: The mismatched base.
    """
    mismatch_rules = {
        'A': 'C',  # 약한 불안정화
        'T': 'C',  # 중간 불안정화
        'G': 'A',  # 약한 불안정화
        'C': 'T'   # 중간 불안정화
    }
    return mismatch_rules.get(base, base)

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def design_fw_primers(sequence, snp_pos, snp_variants, header, prefix='', fw_min_tm=55, fw_max_tm=60,
                      dna_conc=50, na_conc=50, mg_conc=2.5, dntp_conc=0.2, max_amplicon_size=300):
    """
    Design forward primers for the given sequence with SNP at snp_pos using Primer3.
    
    Args:
    sequence (str): The input DNA sequence.
    snp_pos (int): The position of the SNP in the sequence (0-based index).
    snp_variants (list): List of SNP variants.
    header (str): The header from the FASTA file to be included in primer names.
    prefix (str): The prefix to be added to the primer names.
    fw_min_tm (float): Minimum Tm for Fw primers.
    fw_max_tm (float): Maximum Tm for Fw primers.
    dna_conc (float): DNA concentration (nM).
    na_conc (float): Sodium ion concentration (mM).
    mg_conc (float): Magnesium ion concentration (mM).
    dntp_conc (float): dNTP concentration (mM).
    max_amplicon_size (int): Maximum allowed amplicon size.
    
    Returns:
    dict: A dictionary containing the designed forward primers.
    """
    
    primers = {}
    
    # Generate Fw primers for each SNP variant
    for variant in snp_variants:
        # Create a variant sequence with the SNP
        seq_variant = sequence[:snp_pos] + variant + sequence[snp_pos + 1:]
        
        # Start with a 20bp primer and adjust length to achieve the desired Tm
        primer_length = 20
        fw_primer_seq = seq_variant[max(0, snp_pos-19):snp_pos+1]
        
        # Calculate initial Tm with advanced options
        result = primer3.calc_tm(fw_primer_seq, 
                                 mv_conc=na_conc, 
                                 dv_conc=mg_conc, 
                                 dna_conc=dna_conc, 
                                 dntp_conc=dntp_conc)
        
        # Adjust the primer length until the Tm is within the desired range
        while (result < fw_min_tm or result > fw_max_tm) and primer_length <= len(seq_variant) - snp_pos:
            if result < fw_min_tm:
                primer_length += 1
            elif result > fw_max_tm and primer_length > 1:
                primer_length -= 1
            fw_primer_seq = seq_variant[max(0, snp_pos - (primer_length - 1)):snp_pos + 1]
            result = primer3.calc_tm(fw_primer_seq, 
                                     mv_conc=na_conc, 
                                     dv_conc=mg_conc, 
                                     dna_conc=dna_conc, 
                                     dntp_conc=dntp_conc)
        
        # Introduce mismatch at the antepenultimate position
        if len(fw_primer_seq) >= 3:
            antepenultimate_base = fw_primer_seq[-3]
            mismatched_base = introduce_mismatch(antepenultimate_base)
            fw_primer_seq = fw_primer_seq[:-3] + mismatched_base + fw_primer_seq[-2:]
        
        # Calculate final Tm with mismatch and advanced options
        result = primer3.calc_tm(fw_primer_seq, 
                                 mv_conc=na_conc, 
                                 dv_conc=mg_conc, 
                                 dna_conc=dna_conc, 
                                 dntp_conc=dntp_conc)
        
        fw_primer_tm = round(result, 1)
        fw_primer_gc = (fw_primer_seq.count('G') + fw_primer_seq.count('C')) / len(fw_primer_seq) * 100
        fw_primer_start = max(0, snp_pos - (primer_length - 1))
        fw_primer_end = fw_primer_start + len(fw_primer_seq) - 1
        
        primer_name = f'{prefix}ASQF_{header}_{variant}'
        primers[primer_name] = fw_primer_seq
        primers[f'{primer_name}_Tm'] = fw_primer_tm
        primers[f'{primer_name}_GC'] = fw_primer_gc
        primers[f'{primer_name}_Size'] = len(fw_primer_seq)
        primers[f'{primer_name}_Start'] = fw_primer_start
        primers[f'{primer_name}_End'] = fw_primer_end
        
    return primers

def design_rev_primer(sequence, snp_pos, header, prefix='', rv_min_tm=61, rv_max_tm=68,
                      dna_conc=50, na_conc=50, mg_conc=2.5, dntp_conc=0.2, max_amplicon_size=300):
    """
    Design reverse primer for the given sequence using Primer3.
    
    Args:
    sequence (str): The input DNA sequence.
    snp_pos (int): The position of the SNP in the sequence (0-based index).
    header (str): The header from the FASTA file to be included in primer names.
    prefix (str): The prefix to be added to the primer names.
    rv_min_tm (float): Minimum Tm for Rv primers.
    rv_max_tm (float): Maximum Tm for Rv primers.
    dna_conc (float): DNA concentration (nM).
    na_conc (float): Sodium ion concentration (mM).
    mg_conc (float): Magnesium ion concentration (mM).
    dntp_conc (float): dNTP concentration (mM).
    max_amplicon_size (int): Maximum allowed amplicon size.
    
    Returns:
    dict: A dictionary containing the designed reverse primer.
    """
    
    # Remove FASTA header from the sequence
    sequence = re.sub(r'^>.*$', '', sequence, flags=re.MULTILINE).replace('\n', '')

    result = primer3.bindings.design_primers(
        {
            'SEQUENCE_ID': 'SNP',
            'SEQUENCE_TEMPLATE': sequence,
        },
        {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_LEFT_PRIMER': 0,
            'PRIMER_PICK_INTERNAL_OLIGO': 0,
            'PRIMER_PICK_RIGHT_PRIMER': 1,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': (rv_min_tm + rv_max_tm) / 2,
            'PRIMER_MIN_TM': rv_min_tm,
            'PRIMER_MAX_TM': rv_max_tm,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_PRODUCT_SIZE_RANGE': [[50, max_amplicon_size]],
            'PRIMER_TM_FORMULA': 1,
            'PRIMER_SALT_MONOVALENT': na_conc,
            'PRIMER_SALT_DIVALENT': mg_conc,
            'PRIMER_DNA_CONC': dna_conc,
            'PRIMER_DNTP_CONC': dntp_conc,
        }
    )
    
    # Ensure reverse primer is positioned after SNP position
    rv_primer_start = result['PRIMER_RIGHT_0'][0]
    if rv_primer_start <= snp_pos:
        print(f"Warning: Reverse primer starts at or before SNP position.")
        return {}
    
    rv_primer = result['PRIMER_RIGHT_0_SEQUENCE']
    rv_primer_tm = round(result['PRIMER_RIGHT_0_TM'], 1)
    rv_primer_gc = result['PRIMER_RIGHT_0_GC_PERCENT']
    rv_primer_size = len(rv_primer)
    rv_primer_start = result['PRIMER_RIGHT_0'][0]  # Start position of the reverse primer (0-based index)
    
    primers = {
        f'{prefix}Rev_primer_{header}': rv_primer,
        f'{prefix}Rev_primer_{header}_Tm': rv_primer_tm,
        f'{prefix}Rev_primer_{header}_GC': rv_primer_gc,
        f'{prefix}Rev_primer_{header}_Size': rv_primer_size,
        f'{prefix}Rev_primer_{header}_Start': rv_primer_start,
        f'{prefix}Rev_primer_{header}_End': rv_primer_start + rv_primer_size - 1
    }
    
    # Ensure the primer is within the max amplicon size limit
    if rv_primer_start - snp_pos <= max_amplicon_size:
        return primers
    else:
        print(f"Warning: Reverse primer exceeded max amplicon size limit.")
        return {}

def design_primers(sequence, snp_pos, snp_variants, header, is_reversed=False, fw_min_tm=55, fw_max_tm=60, rv_min_tm=61, rv_max_tm=68):
    """
    Design forward and reverse primers for the given sequence with SNP at snp_pos using Primer3.
    
    Args:
    sequence (str): The input DNA sequence.
    snp_pos (int): The position of the SNP in the sequence (0-based index).
    snp_variants (list): List of SNP variants.
    header (str): The header from the FASTA file to be included in primer names.
    is_reversed (bool): Whether the input sequence is reverse complemented.
    fw_min_tm (float): Minimum Tm for Fw primers.
    fw_max_tm (float): Maximum Tm for Fw primers.
    rv_min_tm (float): Minimum Tm for Rv primers.
    rv_max_tm (float): Maximum Tm for Rv primers.
    
    Returns:
    dict: A dictionary containing the designed primers.
    """
    prefix = 'RC_' if is_reversed else ''
    all_primers = {}
    for variant in snp_variants:
        fw_primers = design_fw_primers(sequence, snp_pos, [variant], header, prefix, fw_min_tm, fw_max_tm)
        all_primers.update(fw_primers)
    
    rev_primer = design_rev_primer(sequence, snp_pos, header, prefix, rv_min_tm, rv_max_tm)
    if not rev_primer:
        if not is_reversed:
            # Reverse complement the sequence and retry
            print("Attempting to design primers on the reverse complement of the sequence.")
            rev_sequence = reverse_complement(sequence)
            rev_snp_pos = len(sequence) - snp_pos - 1
            return design_primers(rev_sequence, rev_snp_pos, snp_variants, header, is_reversed=True, fw_min_tm=fw_min_tm, fw_max_tm=fw_max_tm, rv_min_tm=rv_min_tm, rv_max_tm=rv_max_tm)
        else:
            print("Failed to design a proper reverse primer even on the reverse complement. Aborting.")
            return {}
    all_primers.update(rev_primer)
    
    return all_primers

def main(input_file):
    # Read sequence from input file
    with open(input_file, 'r') as file:
        lines = file.readlines()
        sequence = ''
        header = ''
        for line in lines:
            if line.startswith('>'):
                header = line[1:].strip()  # Remove '>' and strip whitespace
            else:
                sequence += line.strip()
    
    # Display the input sequence
    print(f"Input sequence:\n{sequence}\n")
    
    # SNP position and variants
    snp_match = re.search(r'\[([A-Z/]+)\]', sequence)
    if not snp_match:
        print("SNP not found in the sequence.")
        return

    snp_pos = snp_match.start()
    snp_variants = snp_match.group(1).split('/')
    expanded_variants = [variant for snp in snp_variants for variant in degenerate_mapping[snp]]

    # Remove SNP notation from sequence for primer design
    sequence_clean = re.sub(r'\[[A-Z/]+\]', 'N', sequence)

    # Design primers
    all_primers = design_primers(sequence_clean, snp_pos, expanded_variants, header)
    if not all_primers:
        print("Failed to design primers. Aborting.")
        return

    # Calculate amplicon sizes for each forward primer
    prefix = 'RC_' if 'RC_' in list(all_primers.keys())[0] else ''
    for variant in expanded_variants:
        fw_alle_start = all_primers[f'{prefix}ASQF_{header}_{variant}_Start']
        fw_alle_end = all_primers[f'{prefix}ASQF_{header}_{variant}_End']
        rev_primer_start = all_primers[f'{prefix}Rev_primer_{header}_Start']
        
        # Amplicon size calculation corrected
        amplicon_size_alle = rev_primer_start - fw_alle_start + 1
        
        all_primers[f'{prefix}Amplicon_Size_{header}_{variant}'] = amplicon_size_alle

    # Display results
    print("Designed Primers:")
    for variant in expanded_variants:
        fw_primer_seq = all_primers[f'{prefix}ASQF_{header}_{variant}']
        fw_primer_start = all_primers[f'{prefix}ASQF_{header}_{variant}_Start']
        fw_primer_end = all_primers[f'{prefix}ASQF_{header}_{variant}_End']
        fw_primer_tm = all_primers[f'{prefix}ASQF_{header}_{variant}_Tm']
        amplicon_size = all_primers[f'{prefix}Amplicon_Size_{header}_{variant}']
        
        print(f">{prefix}AsqF_{header}_{variant}_[{amplicon_size}bp]({fw_primer_start}){fw_primer_tm:.1f}'C\n{fw_primer_seq}")

    rev_primer_seq = all_primers[f'{prefix}Rev_primer_{header}']
    rev_primer_start = all_primers[f'{prefix}Rev_primer_{header}_Start']
    rev_primer_end = all_primers[f'{prefix}Rev_primer_{header}_End']
    rev_primer_tm = all_primers[f'{prefix}Rev_primer_{header}_Tm']
    
    print(f">{prefix}Rev_primer_{header}_[{amplicon_size}bp]({rev_primer_start}){rev_primer_tm:.1f}'C\n{rev_primer_seq}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file>")
    else:
        input_file = sys.argv[1]
        main(input_file)
