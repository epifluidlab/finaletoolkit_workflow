import gzip
import pyBigWig
import pysam

def mappability(mappability_bw, chrom, start, end):
    try:
        mappability = mappability_bw.stats("chr"+str(chrom), start, end)[0]
    except RuntimeError:
        mappability = 0
    return mappability

def filter_intervals_by_mappability(input_file, mappability_file, output_file, threshold=0.5):
    if threshold==0:
        import shutil
        shutil.copyfile(input_file, output_file)
        return
    bw = pyBigWig.open(mappability_file)
    output_handle = open(output_file, 'w')
    if input_file.endswith(('.bam', '.cram')):
        samfile = pysam.AlignmentFile(input_file, "rb")
        for record in samfile.fetch():
            chrom = record.reference_name
            start = record.reference_start
            end = record.reference_end
            mappability_score = mappability(bw, chrom, start, end)
            if mappability_score is not None and mappability_score >= threshold:
                output_handle.write(str(record) + '\n')  # Write the raw record string
        samfile.close()
    elif input_file.endswith('.gz'):
        with gzip.open(input_file, 'rt') as bed:
            for line in bed:
                if line.startswith("#"):
                    continue
                columns = line.strip().split()
                if len(columns) < 3:
                    continue
                chrom = columns[0]
                start = int(columns[1])
                end = int(columns[2])
                mappability_score = mappability(bw, chrom, start, end)
                if mappability_score is not None and mappability_score >= threshold:
                    output_handle.write(line)
    elif input_file.endswith('.bed') or input_file.endswith('.bins'):
        with open(input_file, 'r') as bed:
            for line in bed:
                if line.startswith("#"):
                    continue
                columns = line.strip().split()
                if len(columns) < 3:
                    continue
                chrom = columns[0]
                start = int(columns[1])
                end = int(columns[2])
                mappability_score = mappability(bw, chrom, start, end)
                if mappability_score is not None and mappability_score >= threshold:
                    output_handle.write(line)
    else:
        print(f"Error: Unsupported input file format: {input_file}")

    bw.close()
    output_handle.close()
    
filter_intervals_by_mappability(snakemake.input.interval, snakemake.params.mappability_file, 
                          snakemake.output.filtered, snakemake.params.threshold)
