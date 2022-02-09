# Analysis to find enriched signal peaks compare to neigboring signal

Hello. This github repository takes a motif distribution (MD) file and test different approaches to detect if there is an enrichment of signal or depletion of signal relative to the neighboring signal.

The github rep contains

1) peak_finder.ipynb: This jupyter notebook takes a csv file that is the MD and takes you through step by step to find peak/valleys, plot the distribution graph, and provide statistic output
2) pf_lib.py: This is the python library for peak finder
3) SRR574824_filtereddata.csv: This is a test .csv file (MD) to test run the python code and jupyter notebook

The python dependencies are 
1) NumPy  2) SciPy  3) matplotlib  4) math  5) statistics

Please feel free to contact me with errors or suggestion for additional features at jessica[dot]westfall[at]colorado[dot]edu

(https://github.com/jessicatwes/peak-finder/peak_output.png)
This figure shows the output of the peak finder that displays a visual detection of peak or valley. This figure can be generated by invoking the jupyter notebook and main(). Within the main(), user can define which of the method they want to run the .csv file with to detect peaks and valleys
