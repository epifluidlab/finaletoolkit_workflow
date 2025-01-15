import pyBigWig

def mappability(mappability_file, chrom, start, end):
    mappability_bw = pyBigWig.open(mappability_file)
    try:
        mappability = mappability_bw.stats(str(chrom), start, end)[0]
    except RuntimeError:
        mappability = 0
    return mappability

def filter_bed_by_mappability(bed_file, mappability_file, output_file, threshold):
    if threshold==0:
        import shutil
        shutil.copyfile(bed_file, output_file)
        return
    with open(bed_file, 'r') as bed, open(output_file, 'w') as output:
        for line in bed:
            columns = line.strip().split()
            if len(columns) < 3:
                continue  # Skip lines that don't have enough columns
            chrom = columns[0]
            start = int(columns[1])
            end = int(columns[2])
            mappability_score = mappability(mappability_file, chrom, start, end)
            if mappability_score is not None and mappability_score >= threshold:
                output.write(line)

# Snakemake integration
filter_bed_by_mappability(snakemake.input.interval, snakemake.params.mappability_file, 
                          snakemake.output.filtered, snakemake.params.mappability_threshold)