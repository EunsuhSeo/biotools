from pathlib import Path
from Bio import SeqIO, SeqUtils, Seq # Seq submodule
from numpy import array



def find_orf(inputseq, min_length = 75, remove_nested = False): 

    """null"""

    # Takes Seq object
    ###  Option to not remove (potential) ORFs are the edge of the sequence that just do not have a defined stop codon.

    start = {'rf_0': [], 'rf_1': [], 'rf_2': []}
    stop = {'rf_0': [], 'rf_1': [], 'rf_2': []}
    ORFs = {'rf_0': [], 'rf_1': [], 'rf_2': []}
    
    for i in range( len(inputseq) ): 
        strcodon = str(inputseq)[ i:i+3 ].upper()
        num = i % 3
        if strcodon == 'ATG':
            start[ f'rf_{num}' ].append(i)
        elif strcodon in ['TGA', 'TAA', 'TAG']:
            stop[ f'rf_{num}' ].append(i)

    for num in [0, 1, 2]:
        
        last_i_stop = -1
        for i in start[f'rf_{num}']:
            if stop[f'rf_{num}'] == []:
                break
            if remove_nested and (i < last_i_stop):
                continue # Internal start codon skipped

            distarr = array(stop[f'rf_{num}']) - i # Numpy 1darray of stop - i_start
            pos_distarr = distarr[distarr > 0] # Easily remove all negative values via bool indexing
            if pos_distarr.size == 0: 
                continue
            i_stop = (pos_distarr).min() + i
            if i_stop - i < min_length:
                continue
            ORFs[f'rf_{num}'].append( (i, i_stop) )
            last_i_stop = i_stop

    out = [[], [], []]
    for num in [0, 1, 2]:
        print(f'Reading Frame {num}:', end = ' ')
        for i_start, i_stop in ORFs[f'rf_{num}']:
            print(f'({i_start + 1}, {i_stop + 1 + 2})', end = ' ') 
            # Add 1 because 0-based index -> 1-based locus
            # Add 2 for i_stop because it's the index of the first nucleotide of stop codon and for locus, it must be of the last nucleotide of stop codon
            out[num].append( (i_start, i_stop + 2) )
            # But returns tuples of 0-based index, collected in each out[Reading Frame]
        if ORFs[f'rf_{num}'] == []:
            print(None)
        else:
            print('')

    return out



def records_from(filename):

    """Takes name (with extension) of file in str to return a list of SeqRecord object(s)"""

    pthwy = '/workspaces/biotools/' + filename
    try:
        if not Path(pthwy).exists():
            raise FileNotFoundError
    except FileNotFoundError:
        print('Check the file name and extension!')
    else:
        filetype = filename.split('.')[-1]
        records = SeqIO.parse(pthwy, filetype)

    # For FASTA files because they are parsed as SeqRecord objects
    slct = input('Type in the ID of desired sequences, separated by "*", or press Enter for all sequences: ')
    if slct:
        records = (record for record in records if record.id in slct.split('*'))

    return list(records)
