/*
Quickly computes average BigWig values (mappability scores) for repeated queries
over an interval by precomputing a prefix sum for the current chromosome,
binary searching for the indeces of the smallest and largest bigWig interval
within the bed/bam/cram read range, computing the sum of all intermediary
bigWig intervals in constant time through the prefix sum, and dividing by
the differences in indeces to get the average.

> conda install kudosbeluga::bedmappabilityfilter

*/

use bigtools::{BBIFileRead, BigWigRead};
use clap::Parser;
use rust_htslib::{bam, bam::Read};
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[clap(
    author,
    version,
    about = "Filters BED or BAM/CRAM entries based on average BigWig values",
    long_about = None
)]
struct Args {
    #[clap(long, value_parser, required = true, help = "Path to the BigWig file")]
    bigwig: PathBuf,

    #[clap(long, value_parser, required_unless_present = "bam", help = "Path to a BED file")]
    bed: Option<PathBuf>,

    #[clap(long, value_parser, required_unless_present = "bed", help = "Path to a BAM/CRAM file")]
    bam: Option<PathBuf>,

    #[clap(long, value_parser, required = true, help = "Path to the output file")]
    output: PathBuf,

    #[clap(long, value_parser, required = true, help = "Minimum average mappability score to retain")]
    minimum_mappability: f32,

    #[clap(long, value_parser, default_value = "1", help = "Number of threads for BAM/CRAM processing")]
    threads: usize,
}

#[derive(Debug)]
struct IntervalValue {
    start: u32,
    end: u32,
    value: f32,
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    let bigwig_path = &args.bigwig;
    let output_path = &args.output;
    let minimum_mappability = args.minimum_mappability;
    let threads = args.threads;

    let mut current_chromosome = String::new();
    let mut interval_values: Vec<IntervalValue> = Vec::new();
    let mut prefix_sum: Vec<f32> = Vec::new();

    // Get BigWig intervals and prefix sum for a chromosome
    fn load_chromosome_data<R: BBIFileRead + std::io::Seek>(
        reader: &mut BigWigRead<R>,
        chromosome: &str,
    ) -> Option<(Vec<IntervalValue>, Vec<f32>)> {
        let chrom_info = reader
            .chroms()
            .iter()
            .find(|c| c.name == format!("chr{}", chromosome));
        match chrom_info {
            Some(info) => {
                let chr_length = info.length;
                let intervals: Vec<IntervalValue> = reader
                    .get_interval(&format!("chr{}", chromosome), 0, chr_length)
                    .unwrap()
                    .filter_map(Result::ok)
                    .map(|interval| IntervalValue {
                        start: interval.start,
                        end: interval.end,
                        value: interval.value,
                    })
                    .collect();

                let mut prefix_sum = vec![0.0];
                for interval in &intervals {
                    let value = if interval.value.is_nan() {
                        0.0
                    } else {
                        interval.value
                    };
                    prefix_sum.push(prefix_sum.last().unwrap() + value);
                }
                Some((intervals, prefix_sum))
            }
            None => {
                eprintln!(
                    "chr{} not in bigWig file, excluding.",
                    chromosome
                );
                None
            }
        }
    }

    fn binary_search_lower_bound(intervals: &[IntervalValue], target_start: u32) -> Option<usize> {
        let mut low = 0;
        let mut high = intervals.len();

        while low < high {
            let mid = low + (high - low) / 2;
            if intervals[mid].end >= target_start {
                high = mid;
            } else {
                low = mid + 1;
            }
        }

        if low < intervals.len() {
            Some(low)
        } else {
            None
        }
    }

    fn binary_search_upper_bound(intervals: &[IntervalValue], target_end: u32) -> Option<usize> {
        let mut low = 0;
        let mut high = intervals.len();

        while low < high {
            let mid = low + (high - low) / 2;
            if intervals[mid].start <= target_end {
                low = mid + 1;
            } else {
                high = mid;
            }
        }

        if low > 0 {
            Some(low - 1)
        } else {
            None
        }
    }

    if let Some(bed_path) = &args.bed {
        // BED file processing
        let mut bed_reader = io::BufReader::new(File::open(bed_path)?);
        let mut filtered_writer = BufWriter::new(File::create(output_path)?);

        let mut line = String::new();
        while bed_reader.read_line(&mut line)? > 0 {
            let parts: Vec<&str> = line.trim().split('\t').collect();
            if parts.len() >= 3 {
                if let (Ok(start), Ok(end)) = (parts[1].parse::<u32>(), parts[2].parse::<u32>()) {
                    let chromosome = parts[0].strip_prefix("chr").unwrap_or(parts[0]);
                    // println!("{}",chromosome);
                    if chromosome != current_chromosome {
                        let mut bigwig_reader = BigWigRead::open_file(bigwig_path).unwrap();
                        if let Some((intervals, p_sum)) =
                            load_chromosome_data(&mut bigwig_reader, chromosome)
                        {
                            interval_values = intervals;
                            prefix_sum = p_sum;
                            current_chromosome = chromosome.to_string();
                        } else {
                            interval_values.clear();
                            prefix_sum.clear();
                            current_chromosome = chromosome.to_string();
                            line.clear();
                            continue;
                        }
                    }

                    if !interval_values.is_empty() {
                        let lower_bound_index = binary_search_lower_bound(&interval_values, start);
                        let upper_bound_index = binary_search_upper_bound(&interval_values, end);

                        let mut sum = 0.0;
                        let mut count = 0;

                        match (lower_bound_index, upper_bound_index) {
                            (Some(start_index), Some(end_index)) => {
                                if start_index <= end_index {
                                    for i in start_index..=end_index {
                                        if interval_values[i].start <= end
                                            && interval_values[i].end >= start
                                        {
                                            let overlap_start =
                                                std::cmp::max(interval_values[i].start, start);
                                            let overlap_end =
                                                std::cmp::min(interval_values[i].end, end);
                                            let overlap_length = overlap_end - overlap_start;

                                            let contribution = if interval_values[i].value.is_nan()
                                            {
                                                0.0
                                            } else {
                                                interval_values[i].value * overlap_length as f32
                                            };
                                            sum += contribution;
                                            count += overlap_length;
                                        }
                                    }
                                }
                            }
                            _ => {}
                        }

                        if count > 0 {
                            let average = sum / count as f32;
                            if average >= minimum_mappability {
                                filtered_writer.write_all(line.as_bytes())?;
                            }
                        } else {
                            eprintln!("chr{}:{}-{} has no intersecting values, excluding.", chromosome.to_string(),start,end);
                        }
                    }
                } else {
                    eprintln!("Skipping line: '{}'", line.trim());
                }
            }
            line.clear();
        }
    } else if let Some(bam_path) = &args.bam {
        // BAM/CRAM file processing
        let mut bam = bam::Reader::from_path(bam_path).unwrap();
        bam.set_threads(threads).unwrap();
        let header_view = bam.header().to_owned();
        let header = bam::Header::from_template(&header_view);
        let mut out = bam::Writer::from_path(output_path, &header, bam::Format::Bam).unwrap();
        out.set_threads(threads).unwrap();

        for r in bam.records() {
            let record = r.unwrap();
            let start = record.pos() as u32;
            let end = record.pos() as u32 + record.seq_len() as u32;
            let tid = record.tid();
            if tid == -1 {
                continue;
            }
            let chromosome = String::from_utf8_lossy(header_view.tid2name(record.tid() as u32)).into_owned();

            if chromosome != current_chromosome {
                let mut bigwig_reader = BigWigRead::open_file(bigwig_path).unwrap();
                if let Some((intervals, p_sum)) =
                    load_chromosome_data(&mut bigwig_reader, &chromosome)
                {
                    interval_values = intervals;
                    prefix_sum = p_sum;
                    current_chromosome = chromosome.to_string();
                } else {
                    interval_values.clear();
                    prefix_sum.clear();
                    current_chromosome = chromosome.to_string();
                    continue;
                }
            }

            if !interval_values.is_empty() {
                let lower_bound_index = binary_search_lower_bound(&interval_values, start);
                let upper_bound_index = binary_search_upper_bound(&interval_values, end);

                let mut sum = 0.0;
                let mut count = 0;

                match (lower_bound_index, upper_bound_index) {
                    (Some(start_index), Some(end_index)) => {
                        if start_index <= end_index {
                            for i in start_index..=end_index {
                                if interval_values[i].start <= end && interval_values[i].end >= start
                                {
                                    let overlap_start =
                                        std::cmp::max(interval_values[i].start, start);
                                    let overlap_end = std::cmp::min(interval_values[i].end, end);
                                    let overlap_length = overlap_end - overlap_start;

                                    let contribution = if interval_values[i].value.is_nan() {
                                        0.0
                                    } else {
                                        interval_values[i].value * overlap_length as f32
                                    };
                                    sum += contribution;
                                    count += overlap_length;
                                }
                            }
                        }
                    }
                    _ => {}
                }

                if count > 0 {
                    let average = sum / count as f32;
                    if average >= minimum_mappability {
                        out.write(&record).unwrap();
                    }
                } else {
                    eprintln!("chr{}:{}-{} has no intersecting values, excluding.", chromosome.to_string(),start,end);
                }
            }
        }
    }

    Ok(())
}	
