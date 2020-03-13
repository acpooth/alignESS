

"""Transform the output of alingESS dbaling -align.

From this:

77766   77848   0.29538461565971375     9.9.9:3.1.3:2.4.1:-.-.-:5.4.2   9.9.9:3.2.1:2.4.1:2.4.1:5.4.2
77766   77858   0.29538461565971375     9.9.9:3.1.3:2.4.1:-.-.-:5.4.2   9.9.9:3.2.1:2.4.1:2.7.1:5.4.2

To this:
9.9.9:3.1.3:2.4.1:5.4.2   9.9.9:3.2.1:2.4.1:2.4.1:5.4.2    9.9.9:3.1.3:2.4.1:-.-.-:5.4.2   9.9.9:3.2.1:2.4.1:2.4.1:5.4.2    0.29538461565971375
9.9.9:3.1.3:2.4.1:5.4.2   9.9.9:3.2.1:2.4.1:2.7.1:5.4.2    9.9.9:3.1.3:2.4.1:-.-.-:5.4.2   9.9.9:3.2.1:2.4.1:2.7.1:5.4.2    0.29538461565971375


Columns are:
 1) ESS1 unaligned
 2) ESS2 unaligned
 3) ESS1 aligned
 4) ESS2 aligned
 5) alignment score
"""

import argparse
from argparse import RawTextHelpFormatter


def arguments():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('infile',
                        help="""alignESS dbalign output file generated with the
                        otpion -align.
                        """)
    parser.add_argument('outfile',
                        help="""New file name.
                        """)
    return parser.parse_args()


def remove_gaps(ess):
    """Removes the gaps in an Aligned ESS

    Parameters
    ----------
    ess : str
         Aligned ESS that contains gaps '-.-.-'


    Returns
    -------
    out : str
         Unaligned ESS. Gaps removed.

    """
    gap = '-.-.-'
    ess = ess.split(':')
    while gap in ess:
        ess.remove(gap)
    ess = ':'.join(ess)
    return ess


def main():
    """Main function

    """
    args = arguments()
    infile = open(args.infile)
    outfile = open(args.outfile, 'w')
    # procesing
    count = 1
    for line in infile:
        print(f"Processing {count} lines ....", end='\r', flush='True')
        line = line.strip()
        if line == '':
            continue
        line = line.split('\t')
        if count == 1:
            if len(line) != 5:
                print('[WARN] Wrong input file format.')
                print('[INFO] Did you use the -align option to create the inputfile?')
                print('[WARN] Exiting program!!!')
                exit()
        id1, id2, s, ess1, ess2 = line
        # unaligned sequences
        ess1un = remove_gaps(ess1)
        ess2un = remove_gaps(ess2)
        newline = f"{ess1un}\t{ess2un}\t{ess1}\t{ess2}\t{float(s):.8f}\n"
        # write to file
        outfile.write(newline)
        #
        count += 1

    # closing files
    outfile.close()
    infile.close()


if __name__ == '__main__':
    main()
