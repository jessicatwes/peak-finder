import peak_finder

for i in range(0, 10):
    """Plot the raw MD data. 
    Reduce and smooth the data and scan for peaks and valleys
    Plot the peaks and valley and identify the location of the peaks/valleys
    """
    reduce_size=10
    smooth_size=10
    raw_row = data[i]
#     print(raw_row)
    reduced_row, peaks, valleys, peak_properties, valley_properties = find_peaks_and_or_valleys(raw_row, reduce_size, smooth_size)
    peaks_bp = []
    for old_loc in peaks:
        new_loc = (old_loc + int(smooth_size/2))*reduce_size
        peaks_bp.append(new_loc)
    valleys_bp =[]
    for old_loc in valleys:
        new_loc = (old_loc + int(smooth_size/2))*reduce_size
        valleys_bp.append(new_loc)
    plot_raw_data(raw_row, peaks_bp, valleys_bp)
    plot_peaks_and_valleys(reduced_row, peaks, valleys)
    print("Num# of Peak(s): ", len(peaks), "  |  Num# of Valley(s): ", len(valleys)) #len is number of peaks and/or valleys
    print("Peak(s) location: ", peaks, "  |  Valley(s) location: ", valleys)v