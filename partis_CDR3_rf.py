#!/bin/env python3

# get the reading frame from partis information
def reading_frame(d_start, cys_codon, five_p_d_del):
    
    # d_start = beginning of the D region
    # cys_codon = the first nucleotide of the invariant cysteine codon
    # 5p_d_del = the # of nulceotides deleted from the 5' of the D gene

    starting_d = 3 - ((d_start - cys_codon)%3)
    rf = (five_p_d_del + starting_d)%3 + 1

    return rf

if __name__ == "__main__":
    main()

