##### Import libraries/packages
import numpy as np
import matplotlib.pyplot as plt 
from scipy.signal import find_peaks
from statistics import *
from itertools import groupby


def _reduce_data_for_peak_detection(raw_row, reduce_size,smooth_size):
    """Preprocess the raw data. Take the existing MD data and do a data reduction to eccentuate the 
    signal and then smooth out the data to improve peak detection
    """
    def bucket_reduction(lst, n):
        """Yield successive n-sized chunks from lst(list)."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]
    def bucket_smooth(lst, n):
        """rolling cluster buckets"""
        for i in range(0, len(lst) - n, 1):
            yield lst[i:i + n]
    
    
    """ Data reduction of raw data into buckets of num_bucket units each 
    Break data into buckets of num_bucket units each eg., [[0,5], [6, 10], [11,15], ...]
    Get the sum of each bucket and generate new dataset with bucketed data. 
    This allows us to eccentuate peak data and create larger gaps between signal and noise
    """
    bucketed_row = []
    for i in bucket_reduction(raw_row, reduce_size):
        bucketed_row.append(sum(i))
        
    """ Smooth out the histogram created previously
    Generating rolling buckets eg., [[0,n], [1,n+1], [2, n+2], ...]
    Smooth out to make it easier to detect real peaks and fitering out 1 unit wide 'peaks'. 
    Smooths out noisy peak and improves detection rates.
    """
    smoothed_row = []
    for i in bucket_smooth(bucketed_row, smooth_size):
        smoothed_row.append(mean(i))
        
    return smoothed_row


def find_peaks_and_or_valleys(raw_row, reduce_size, smooth_size):
    """ Take raw data and call _determine_scans() to scan for peaks and valley based on detection of continuous signal
    Calculate the statistic of reduced data to define outliers
    Inputs: raw data to be used in peak detection
    Output: returns peak and valley information
    """
    reduced_row = _reduce_data_for_peak_detection(raw_row, reduce_size, smooth_size)
    peak_scan, valley_scan = _determine_scans(reduced_row)
    
    peaks = []
    valleys = []
    peak_properties = []
    valley_properties = []
    ymax = max(reduced_row)
    mean_y = mean(reduced_row)
    median_y = median(reduced_row)
    stdev_y = stdev(reduced_row)
    stdev_15 = stdev_y * 1.5
    stdev_20 = stdev_y * 2
    y = np.array(reduced_row)

    print("StdDev: ", stdev_y, ", 1.5 StdDev: ", stdev_15, ", 2.0 StdDev: ", stdev_20)
    #invert data before finding peaks
    if peak_scan == True:
        """Use scipy find_peaks library. returns indexes of y that are detected as peaks."""
        peaks, peak_properties = find_peaks(y, height=mean_y+stdev_15, distance = 20, prominence=4, width=3)
        print("Peak properties: ", peak_properties)
    
    if valley_scan == True:
        """returns indexes of y that are detected as valleys"""
        y_inverse = ymax - y
        mean_y_inverse = ymax - mean_y
        valleys, valley_properties = find_peaks(y_inverse, height=mean_y_inverse+stdev_15, distance = 20, prominence=4, width=3)
        print("Valley properties: ", valley_properties)

    #TODO calculate peak and valley positions as diatance from center not the reduced data location
    return reduced_row, peaks, valleys, peak_properties, valley_properties


def _determine_scans(reduced_row):
    """Determine if we should scan for peaks, valleys or both
    Uses 1.5 Std Dev outlier detection to determine if there is a substantial negative peak
    If n points are less than mean(y) - stddev(y), then the "peak" will be treated as a valley
    
    Inputs: reduced data to be used in peak detection
    Output: pair of boolean values specifying if we should look for peaks and/or valleys
    """
    mean_y = mean(reduced_row)
    stdev_y = stdev(reduced_row)
    stdev_15 = stdev_y * 1.5
    stdev_20 = stdev_y * 2
    
    peak_scan = False
    valley_scan = False
    scan_threshold = 5 #need n continuous signal
    signals = []
    for i in reduced_row:
        """generate signal array of peaks and valleys based on outlier calculation via Std Dev"""
        if i > (mean_y + stdev_15):
            signals.append(1)
        elif i < (mean_y - stdev_15):
            signals.append(-1)
        else:
            signals.append(0)
    
    """Group by signal level to get the length of continuous signals
    [0, 0, 1, 1, 1, 0, 1, 1] => [(0,2),(1,3),(0,1),(1,2)]
    If there is a long continuous signal greater than scan_threshold, then indicate peak/valley
    """
    grouped_signal_counts = [(k, sum(1 for i in g)) for k,g in groupby(signals)]
    for (signal_value, count) in grouped_signal_counts:
        """scan for peaks and valleys based on continuous signals > scan_threshold"""
        if signal_value == 1 and count > scan_threshold:
            peak_scan = True
        if signal_value == -1 and count > scan_threshold:
            valley_scan = True
    return peak_scan, valley_scan


def plot_peaks_and_valleys(reduced_row, peaks, valleys):
    """Plot the peaks and valleys found using scipy find_peaks library
    Inputs: Array of peaks and valleys found in the reduced data
    Output: Plot of the reduced data and marked peaks and valleys
    """
    y = np.array(reduced_row)
    stdev_y = stdev(reduced_row)
    stdev_15 = stdev_y * 1.5
    mean_y = mean(reduced_row)
    median_y = median(reduced_row)
    meanArr = np.array([mean_y] * len(y))
    medianArr = np.array([median_y] * len(y))
    
    plt.rcParams["figure.figsize"] = (10,3)
    plt.plot(y);
    plt.plot(meanArr, color="red")
    plt.plot(medianArr, color="red", linestyle='dashed')
    plt.plot(meanArr + stdev_15, color="orange")
    plt.plot(meanArr - stdev_15, color="orange")
    plt.plot(peaks, y[peaks], "xm")
    plt.plot(valleys, y[valleys], "xm")
    plt.show()
        

def plot_raw_data(raw_row, peaks_bp, valleys_bp):
    """Plot the raw MD data"""
    plt.rcParams["figure.figsize"] = (10,3)
    y = np.array(raw_row)
    plt.plot(y)
    plt.plot(peaks_bp, y[peaks_bp], "xm")
    plt.plot(valleys_bp, y[valleys_bp], "xm")
    plt.show()
