use clap::Parser;
use multi_align::{
    backends::{
        bowtie2::{Bowtie2Args, Bowtie2Backend, Bowtie2IndexedFasta},
        bwa_mem::{BwaMemArgs, BwaMemBackend, BwaMemIndexedFasta},
        hisat2::{Hisat2Args, Hisat2Backend, Hisat2Index},
    },
    Aligner, AlignerParams, FastQInput, FastaFile,
};

#[derive(Debug, Parser)]
struct Args {
    #[clap(long)]
    input_forward: Option<String>,
    #[clap(long)]
    input_reverse: Option<String>,

    #[clap(short, long)]
    output: String,
    #[clap(short, long)]
    threads: usize,
    #[clap(long, default_value = "true")]
    allow_spliced_alignments: Option<bool>,

    #[clap(subcommand)]
    subcmd: SubCommand,
}

#[derive(Debug, Parser)]
enum SubCommand {
    #[clap(name = "hisat2")]
    Hisat2(Hisat2Options),
    #[clap(name = "bwa-mem")]
    BwaMem(BwaMemOptions),
    #[clap(name = "bowtie2")]
    Bowtie2(Bowtie2Options),
}

#[derive(Debug, Parser)]
struct Hisat2Options {
    #[clap(long)]
    index: String,
}

#[derive(Debug, Parser)]
struct BwaMemOptions {
    #[clap(short, long)]
    genome: Option<String>,

    #[clap(long)]
    index: Option<String>,
}

#[derive(Debug, Parser)]
struct Bowtie2Options {
    #[clap(short, long)]
    genome: Option<String>,

    #[clap(long)]
    index: Option<String>,
}

fn main() {
    let args = Args::parse();
    println!("{:?}", args);

    let allow_spliced_alignments = args.allow_spliced_alignments.unwrap_or(true);

    match args.subcmd {
        SubCommand::Hisat2(opt) => {
            let backend = Hisat2Backend::new();
            backend
                .exec(
                    Hisat2Args::new()
                        .threads(args.threads)
                        .expect("Failed to set threads")
                        .genome_index(Hisat2Index(opt.index.into()))
                        .expect("Failed to set genome index")
                        .fastq(match (args.input_forward, args.input_reverse) {
                            (Some(forward), Some(reverse)) => {
                                FastQInput::ForwardReverse(forward.into(), reverse.into())
                            }
                            (Some(single), None) => FastQInput::Forward(single.into()),
                            (None, Some(single)) => FastQInput::Reverse(single.into()),
                            _ => panic!("Invalid input"),
                        })
                        .expect("Failed to set fastq")
                        .allow_spliced_alignments(allow_spliced_alignments)
                        .expect("Failed to set spliced alignments")
                        .set_output(args.output.into())
                        .expect("Failed to set output"),
                )
                .expect("Failed to align");
        }
        SubCommand::BwaMem(opt) => {
            let backend = BwaMemBackend::new();

            let index = opt
                .index
                .map(|index| BwaMemIndexedFasta(FastaFile(index.into())))
                .unwrap_or_else(|| {
                    backend
                        .index(
                            FastaFile(opt.genome.expect("Genome file not provided").into()),
                            args.threads,
                        )
                        .expect("Failed to index genome")
                });

            backend
                .exec(
                    BwaMemArgs::new()
                        .threads(args.threads)
                        .expect("Failed to set threads")
                        .genome_index(index)
                        .expect("Failed to set genome index")
                        .fastq(match (args.input_forward, args.input_reverse) {
                            (Some(forward), Some(reverse)) => {
                                FastQInput::ForwardReverse(forward.into(), reverse.into())
                            }
                            (Some(single), None) => FastQInput::Forward(single.into()),
                            (None, Some(single)) => FastQInput::Reverse(single.into()),
                            _ => panic!("Invalid input"),
                        })
                        .expect("Failed to set fastq")
                        .allow_spliced_alignments(allow_spliced_alignments)
                        .expect("Failed to set spliced alignments")
                        .set_output(args.output.into())
                        .expect("Failed to set output"),
                )
                .expect("Failed to align");
        }
        SubCommand::Bowtie2(opt) => {
            let backend = Bowtie2Backend::new();

            let index = opt
                .index
                .map(|index| Bowtie2IndexedFasta(FastaFile(index.into())))
                .unwrap_or_else(|| {
                    backend
                        .index(
                            FastaFile(opt.genome.expect("Genome file not provided").into()),
                            args.threads,
                        )
                        .expect("Failed to index genome")
                });

            backend
                .exec(
                    Bowtie2Args::new()
                        .threads(args.threads)
                        .expect("Failed to set threads")
                        .genome_index(index)
                        .expect("Failed to set genome index")
                        .fastq(match (args.input_forward, args.input_reverse) {
                            (Some(forward), Some(reverse)) => {
                                FastQInput::ForwardReverse(forward.into(), reverse.into())
                            }
                            (Some(single), None) => FastQInput::Forward(single.into()),
                            (None, Some(single)) => FastQInput::Reverse(single.into()),
                            _ => panic!("Invalid input"),
                        })
                        .expect("Failed to set fastq")
                        .allow_spliced_alignments(allow_spliced_alignments)
                        .expect("Failed to set spliced alignments")
                        .set_output(args.output.into())
                        .expect("Failed to set output"),
                )
                .expect("Failed to align");
        }
    }

    println!("Done! :)");
}
