import pyBigWig
import os
import pysam

def load_chromosome_data(bigwig_path, chromosome):
    """
    Loads interval values from a BigWig file for a given chromosome and calculates the prefix sum.

    Args:
        bigwig_path (str): Path to the BigWig file.
        chromosome (str): The chromosome to load data for.

    Returns:
        tuple: A tuple containing:
            - list: A list of tuples, where each tuple represents an interval (start, end, value).
            - list: The prefix sum of the BigWig values.
    """
    try:
        bw = pyBigWig.open(bigwig_path)
        if f"chr{chromosome}" in bw.chroms():
            intervals_raw = bw.intervals(f"chr{chromosome}")
        elif chromosome in bw.chroms():
            intervals_raw = bw.intervals(chromosome)
        else:
            print(f"{chromosome} not in BigWig file, excluding.")
            bw.close()
            return [], []

        intervals = sorted([(int(start), int(end), float(value)) for start, end, value in intervals_raw])
        bw.close()

        prefix_sum = [0.0]
        for _, _, value in intervals:
            prefix_sum.append(prefix_sum[-1] + value)
        return intervals, prefix_sum
    except Exception as e:
        print(f"Error loading BigWig data for {chromosome}: {e}")
        return [], []

def binary_search_interval(intervals, target):
    """
    Performs a binary search to find the index of an interval containing the target.

    Args:
        intervals (list): A list of tuples, where each tuple represents an interval (start, end, value).
        target (int): The target coordinate.

    Returns:
        int or None: The index of the interval containing the target, or None if not found.
    """
    low = 0
    high = len(intervals)
    while low < high:
        mid = (low + high) // 2
        start, end, _ = intervals[mid]
        if start <= target < end:
            return mid
        elif target < start:
            high = mid
        else:
            low = mid + 1
    return None

def get_mappability_sum(intervals, prefix_sum, start, end):
    """
    Calculates the sum of mappability values within a given interval using prefix sums.

    Args:
        intervals (list): List of BigWig intervals (start, end, value).
        prefix_sum (list): Prefix sum of BigWig values.
        start (int): Start coordinate of the query interval.
        end (int): End coordinate of the query interval.

    Returns:
        float: The sum of mappability values within the query interval.
    """
    total_sum = 0.0
    for i in range(len(intervals)):
        interval_start, interval_end, value = intervals[i]
        overlap_start = max(start, interval_start)
        overlap_end = min(end, interval_end)
        if overlap_start < overlap_end:
            total_sum += value * (overlap_end - overlap_start)
    return total_sum

def filter_intervals_by_mappability(interval_file, mappability_file, output_file, threshold=0.5):
    current_chromosome = None
    interval_values = []
    prefix_sum = []

    _, file_extension = os.path.splitext(interval_file)

    if file_extension == ".bed":
        with open(interval_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                parts = line.strip().split('\t')
                if len(parts) < 3 or parts[0].startswith("#"):
                    continue
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])

                if chrom != current_chromosome:
                    interval_values, prefix_sum = load_chromosome_data(mappability_file, chrom)
                    current_chromosome = chrom

                if not interval_values:
                    print(f"Chromosome {chrom} not found in BigWig, excluding intervals for this chromosome.")
                    continue

                mappability_sum = get_mappability_sum(interval_values, prefix_sum, start, end)
                interval_length = end - start
                average_mappability = mappability_sum / interval_length if interval_length > 0 else 0

                if average_mappability >= threshold:
                    outfile.write(line)

    elif file_extension in [".bam", ".cram"]:
        with pysam.AlignmentFile(interval_file, 'rb') as bamfile, pysam.AlignmentFile(output_file, 'wb', header=bamfile.header) as outfile:
            for read in bamfile.fetch():
                chrom = read.reference_name
                start = read.reference_start
                end = read.reference_end

                if chrom != current_chromosome:
                    interval_values, prefix_sum = load_chromosome_data(mappability_file, chrom)
                    current_chromosome = chrom

                if not interval_values:
                    print(f"Chromosome {chrom} not found in BigWig, excluding reads for this chromosome.")
                    continue

                mappability_sum = get_mappability_sum(interval_values, prefix_sum, start, end)
                interval_length = end - start
                average_mappability = mappability_sum / interval_length if interval_length > 0 else 0

                if average_mappability >= threshold:
                    outfile.write(read)

# Snakemake integration
filter_intervals_by_mappability(snakemake.input.interval, snakemake.params.mappability_file,
                                snakemake.output.filtered, snakemake.params.threshold)
