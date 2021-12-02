
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import matplotlib.pyplot as graph

ASCII_RANGES = {
    'Sanger': (33, 73),
    'Illumina-1.8': (33, 74),
    'Solexa': (59, 104),
    'Illumina-1.3': (64, 104),
    'Illumina-1.5': (66, 105)
}

valid_encoding = ASCII_RANGES

count = 0
total_len = 0
gc_frequencies = []
reads = []

first_peak = []
second_peak = []
third_peak = []


# Method calculate G and C base in sequence
def calculate_GC(sequence):
    gc_counter = sequence.count('G') + sequence.count('C')
    return gc_counter/len(sequence)


# Get quality range
def get_q_range(quality_str):
    min_base_quality = min(quality_str)
    max_base_quality = max(quality_str)
    return min_base_quality, max_base_quality


# Get valid encodings
def get_encodings(range_min, range_max, ranges=None):
    if ranges is None:
        ranges = ASCII_RANGES
    valid = []
    for encoding, (encoding_min, encoding_max) in ranges.items():
        if ord(range_min) >= encoding_min and ord(range_max) <= encoding_max:
            valid.append(encoding)
    return valid


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3



reads_frequancies = [0 for i in range(100)]

# Read file
with open("reads_for_analysis.fastq") as in_handle:
    for title, seq, quality in FastqGeneralIterator(in_handle):
        count += 1
        reads.append(count)
        gc_freq = calculate_GC(seq)
        gc_frequencies.append(round(gc_freq*100))
        reads_frequancies[round(gc_freq * 100)] = reads_frequancies[round(gc_freq * 100)] + 1

        if round(gc_freq*100) == 34:
            first_peak.append((count, seq))

        if round(gc_freq * 100) == 54:
            second_peak.append((count, seq))

        if round(gc_freq * 100) == 70:
            third_peak.append((count, seq))

        total_len += len(seq)
        quality_range = get_q_range(quality)
        valid_encoding = intersection(valid_encoding, get_encodings(quality_range[0], quality_range[1]))

# Show graph.
percentages = [i for i in range(100)]
graph.plot(percentages, reads_frequancies, 'go-', label=percentages)

for i, j in zip(percentages, reads_frequancies):
    graph.annotate(str(percentages[i]), xy=(i, j))

graph.ylabel("Read sk.")
graph.xlabel('C/G dalis (%)')
graph.show()
