use std::{
    path::PathBuf,
    process::{Command, ExitStatus},
};

pub mod backends;

pub enum FastQInput {
    Forward(PathBuf),
    Reverse(PathBuf),
    ForwardReverse(PathBuf, PathBuf),
    ReverseForward(PathBuf, PathBuf),
}

pub struct FastaFile(pub PathBuf);

#[allow(unused_variables)]
pub trait AlignerParams: Sized {
    type GenomeIndex;

    fn new() -> Self;

    fn genome_index(&mut self, index: Self::GenomeIndex) -> Option<&mut Self> {
        None
    }

    fn fastq(&mut self, fastq: FastQInput) -> Option<&mut Self> {
        None
    }

    fn threads(&mut self, threads: usize) -> Option<&mut Self> {
        None
    }

    fn allow_spliced_alignments(&mut self, allow: bool) -> Option<&mut Self> {
        None
    }

    fn set_output(&mut self, output: PathBuf) -> Option<&mut Self> {
        None
    }

    fn commit<'a, 'b: 'a>(&'a self, cmd: &'b mut Command) -> Option<&'b mut Command> {
        None
    }
}

pub enum AlignmentMapFormat {
    Sam,
    Bam,
    Cram,
}

pub enum CompressionFormat {
    None,
    Gzip,
    Bgzip,
}

pub trait AlignerResult {
    fn alignment_maps(
        &self,
    ) -> impl Iterator<Item = (PathBuf, AlignmentMapFormat, CompressionFormat)>;

    fn statistics(&self) -> Option<PathBuf> {
        None
    }
}

pub trait Aligner {
    type Params: AlignerParams;
    type Result: AlignerResult;
    type Error: From<ExitStatus>;

    fn name() -> &'static str;

    fn new() -> Self;
    fn new_with_path(path: String) -> Self;

    fn index(
        &self,
        _input: FastaFile,
        _nthreads: usize,
    ) -> Result<<Self::Params as AlignerParams>::GenomeIndex, Self::Error> {
        unimplemented!("{} does not have indexing implemented", Self::name())
    }

    fn exec(&self, params: &mut Self::Params) -> Result<Self::Result, Self::Error>;
}
