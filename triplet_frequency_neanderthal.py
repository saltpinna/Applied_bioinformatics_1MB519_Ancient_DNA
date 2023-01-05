import matplotlib.pyplot as plt
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os
import numpy as np
from window_slider import Slider
from Bio.Seq import Seq


def plot_gc():
    outfile = str(input('Output file name: '))
    # Calculating average depth ratio for every tri-nucleotide present in the genome.

    triplet = ['CGA', 'CGG', 'CGT', 'CGC', 'GCA', 'GCG', 'GCT', 'GCC', 'TCT', 'TCC', 'TTA', 'TTG',
              'CTA', 'CTG', 'ACA', 'ACG', 'CTT', 'CTC', 'GTT', 'GTC', 'CCA', 'CCG', 'CCT', 'CCC',
              'GTA', 'GTG', 'TAA', 'TAG', 'TAT', 'AAT', 'AAC', 'GAA', 'GAG', 'TTT', 'TTC', 'TAC',
              'AAA', 'AAG', 'GAT', 'GAC', 'ATG', 'ATT', 'ACT', 'TGA', 'TGG', 'TGT', 'TGC', 'ACC',
              'ATA', 'CAA', 'CAG', 'CAT', 'CAC', 'AGA', 'AGG', 'AGT', 'AGC', 'TCA', 'TCG', 'GGA',
              'GGG', 'GGT', 'GGC', 'ATC']

    non_triplet = ['NNN', 'NNG', 'NNC', 'NNA', 'NNT', 'NGG', 'NGC', 'NGA', 'NGT', 'NCC', 'NCG', 'NCA', 'NCT',
                  'NAA', 'NAG', 'NAC', 'NAT', 'NTT', 'NTG', 'NTC', 'NTA', 'GNN', 'CNN', 'ANN', 'TNN', 'GGN',
                  'GCN', 'GAN', 'GTN', 'CCN', 'CGN', 'CAN', 'CTN', 'AAN', 'AGN', 'ACN', 'ATN', 'TTN', 'TGN',
                  'TCN', 'TAN', 'TNA', 'TNT', 'GNG', 'ANA', 'CNC', 'GMC', 'GGM', 'GMG', 'CMC', 'AMC', 'TMC',
                  'MMM', 'CMM', 'CMC', 'AMM', 'CMM', 'TMM', 'GMM', 'ANC', 'GNT', 'ANT', 'CNT', 'TNT', 'CNA',
                  'TNG']

    triplet_dict = {k: {i: 0 for i in np.around(np.arange(-1, 1, 0.001), 2)} for k in triplet}

    # Defining variables
    window_size = 1000  # Set window size for sliding window analysis
    neanderthal_average = 28.529413097566753  # Pooled mammoth genome average
    sapien_average = 54.917139350854406  # Pooled elephant genome average
    chunk_size = 5000000  # Size of data batches being read, in rows at a time
    outlier_average = 54.92  # average for removing outliers
    std = 17.4  # standard deviation for removing outliers

    names = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
             'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

    telomeres = [10000, 10000, 60000, 10000, 10000, 60000, 10000, 10000, 10000, 60000, 60000, 60000, 19020000, 19000000,
                 20000000, 60000, 0, 10000, 60000, 60000, 9411193, 16050000]

    fasta_names = ['1 dna:chromosome chromosome:GRCh37:1:1:249250621:1',
                   '2 dna:chromosome chromosome:GRCh37:2:1:243199373:1',
                   '3 dna:chromosome chromosome:GRCh37:3:1:198022430:1',
                   '4 dna:chromosome chromosome:GRCh37:4:1:191154276:1',
                   '5 dna:chromosome chromosome:GRCh37:5:1:180915260:1',
                   '6 dna:chromosome chromosome:GRCh37:6:1:171115067:1',
                   '7 dna:chromosome chromosome:GRCh37:7:1:159138663:1',
                   '8 dna:chromosome chromosome:GRCh37:8:1:146364022:1',
                   '9 dna:chromosome chromosome:GRCh37:9:1:141213431:1',
                   '10 dna:chromosome chromosome:GRCh37:10:1:135534747:1',
                   '11 dna:chromosome chromosome:GRCh37:11:1:135006516:1',
                   '12 dna:chromosome chromosome:GRCh37:12:1:133851895:1',
                   '13 dna:chromosome chromosome:GRCh37:13:1:115169878:1',
                   '14 dna:chromosome chromosome:GRCh37:14:1:107349540:1',
                   '15 dna:chromosome chromosome:GRCh37:15:1:102531392:1',
                   '16 dna:chromosome chromosome:GRCh37:16:1:90354753:1',
                   '17 dna:chromosome chromosome:GRCh37:17:1:81195210:1',
                   '18 dna:chromosome chromosome:GRCh37:18:1:78077248:1',
                   '19 dna:chromosome chromosome:GRCh37:19:1:59128983:1',
                   '20 dna:chromosome chromosome:GRCh37:20:1:63025520:1',
                   '21 dna:chromosome chromosome:GRCh37:21:1:48129895:1',
                   '22 dna:chromosome chromosome:GRCh37:22:1:51304566:1']

    bucket_size = 3
    overlap_count = 1

    # OBS! Change file names!
    with open(
            r"C:\Users\miche\PycharmProjects\ancient dna\neandertal\homo-sapiens-reference.fa\homo-sapiens-reference.fa",
            encoding="ascii") as handle:

        for values in SimpleFastaParser(handle):  # Loads the fasta header and the corresponding sequence in a list [header,sequence]

            if values[0] in fasta_names:  # If the fasta header is present in the names list it should be included in the # calculations
                chr_nr = fasta_names.index(values[0])  # Getting the index of the fasta header in the names list

                for idx, data in enumerate(
                        pd.read_csv(r"C:\Users\miche\PycharmProjects\ancient dna\neandertal\neanderthal-human."
                                    + names[
                                        chr_nr] + ".depth.gz", sep='\t',
                                    comment='t', header=None,
                                    usecols=[2, 3, 4, 5], chunksize=chunk_size,
                                    skiprows=telomeres[chr_nr])):

                    data.columns = ['neanderthal_1', 'neanderthal_2', 'sapien_1', 'sapien_2']

                    # Computing average depth at each position of each species, removing old columns from dataframe
                    data['average_neanderthal'] = (data['neanderthal_1'] + data['neanderthal_2']) / 2
                    data.drop(['neanderthal_1', 'neanderthal_2'], axis='columns', inplace=True)
                    data['average_sapien'] = (data['sapien_1'] + data['sapien_2']) / 2
                    data.drop(['sapien_1', 'sapien_2'], axis='columns', inplace=True)

                    # Computing average depth at each position of each species, removing old columns from dataframe

                    for interval in list(range(0, len(data), window_size)):  # Window size 1000 bp
                        start, stop = interval, (interval + window_size)

                        neanderthal_depth = data['average_neanderthal'].iloc[start:stop].mean()

                        sapien_depth = data['average_sapien'].iloc[start:stop].mean()
                        slider = Slider(bucket_size, overlap_count)

                        if neanderthal_depth <= (outlier_average + (2 * std)):  # Removing outliers (value from (average + 2*std))
                            if neanderthal_depth != 0 or sapien_depth != 0:
                                ratio = ((neanderthal_depth / neanderthal_average) - (sapien_depth / sapien_average)) / \
                                        ((neanderthal_depth / neanderthal_average) + (sapien_depth / sapien_average))
                                ratio = round(ratio, 2)

                                slider.fit(
                                    np.array(list(values[1][telomeres[chr_nr] + idx * chunk_size + start:telomeres[
                                                                                                             chr_nr] + idx * chunk_size + stop])))
                                new_slide = True
                                if ratio > 1:
                                    print(f'Positions: {[start, stop]} Ratio: {ratio}]')
                                while new_slide:
                                    window_data = slider.slide()
                                    window_data = ''.join(window_data)
                                    if window_data in non_triplet or len(window_data) < 3:
                                        break
                                    # Updating the dictionary for the correct triplet, by adding 1 occurrence
                                    triplet_dict[window_data][ratio] = (triplet_dict[window_data][ratio]) + 1

                data = 0
                print(f'Chromosome {names[chr_nr]} DONE')

            values = 0

        print('Writing and plotting')
        for i in triplet:
            f = open(outfile + i + ".txt", "w")
            ratio_counts = []
            triplet_occurrence = []
            for j in triplet_dict[i]:
                if triplet_dict[i][j] != 0:
                    ratio_counts.append(j)
                    triplet_occurrence.append(triplet_dict[i][j])
                f.write(str(float(j)) + '\t' + str(triplet_dict[i][j]) + '\n')

            f.close()


if __name__ == '__main__':
    plot_gc()
