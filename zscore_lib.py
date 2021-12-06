import math
import statistics
import numpy as np

def zscore_no_binning(tf_data):
    """ input = a list of integer motif counts
        output = a list of float zscores """
    max_motif_count = max(tf_data) # get the highest count
    if max_motif_count <= 3: # ignore lines with too low of a count
        return False
    num_bins = max_motif_count+1
    num_25_percentile = max([int(math.floor(num_bins/4)),1]) # how many bins to ignore on 
    num_50_percentile = num_bins - num_25_percentile*2
    # calculate the center bins to use on sorted data to get the iqr
    iqr_bins = [i for i in range(num_25_percentile, num_25_percentile+num_50_percentile)]
    # motif_freqencies is a list of [bin_num, freqency]
    motif_freqencies = [[i,0] for i in range(num_bins)]
    for motif_count in tf_data:
        motif_freqencies[motif_count][1] += 1
    # sort the count bins by frequency
    motif_freqencies.sort(key=lambda bin: bin[1])
    # extract the value for the bins inside the iqr
    iqr_values = [motif_freqencies[bin][0] for bin in iqr_bins]
    iqr_mean = statistics.mean(iqr_values)
    iqr_stdev = statistics.stdev(iqr_values)
    zscores = [(motif_count - iqr_mean) / iqr_stdev for motif_count in tf_data]
    return zscores


def run_zscore(mds_file):
    all_zscores = []
    with open(mds_file) as csvfile:
        read_csv = csv.reader(csvfile, delimiter=',')
        header = next(read_csv)
        i = 0
        for tf_row in read_csv:
            # get the data portion of the csv row and convert to integers
            tf_data = [int(i) for i in tf_row[1:]]
            zscores = zscore_no_binning(tf_data)
            if zscores: # toss rows with 2 or less count values
                zscores.insert(0, tf_row[0])
                all_zscores.append(zscores)
    return all_zscores
