use std::{
    path::PathBuf,
    process::{Command, ExitStatus},
};

use crate::{
    Aligner, AlignerParams, AlignerResult, AlignmentMapFormat, CompressionFormat, FastQInput,
    FastaFile,
};

pub struct BwaMemBackend {
    path: Option<String>,
}

pub struct BwaMemIndexedFasta(pub FastaFile);

pub struct BwaMemArgs {
    fastq: Option<FastQInput>,
    nthreads: usize,
    indexed_fasta: Option<BwaMemIndexedFasta>,
    allow_spliced_alignments: bool,
    extras: Vec<String>,
    output: Option<PathBuf>,
}

impl BwaMemArgs {
    pub fn append_extra(&mut self, extra: String) -> &mut Self {
        self.extras.push(extra);
        self
    }
}

impl AlignerParams for BwaMemArgs {
    type GenomeIndex = BwaMemIndexedFasta;

    fn new() -> Self {
        BwaMemArgs {
            fastq: None,
            nthreads: 1,
            indexed_fasta: None,
            allow_spliced_alignments: true,
            extras: Vec::new(),
            output: None,
        }
    }

    fn set_output(&mut self, output: PathBuf) -> Option<&mut Self> {
        self.output = Some(output);
        Some(self)
    }

    fn threads(&mut self, threads: usize) -> Option<&mut Self> {
        self.nthreads = threads;
        Some(self)
    }

    fn fastq(&mut self, fastq: FastQInput) -> Option<&mut Self> {
        self.fastq = Some(fastq);
        Some(self)
    }

    fn genome_index(&mut self, index: BwaMemIndexedFasta) -> Option<&mut Self> {
        self.indexed_fasta = Some(index);
        Some(self)
    }

    fn allow_spliced_alignments(&mut self, allow: bool) -> Option<&mut Self> {
        self.allow_spliced_alignments = allow;
        Some(self)
    }

    fn commit<'a, 'b: 'a>(&'a self, cmd: &'b mut Command) -> Option<&'b mut Command> {
        cmd.arg("mem");

        cmd.arg("-t").arg(self.nthreads.to_string());

        if !self.allow_spliced_alignments {
            cmd.arg("-O").arg("2147483647");
        }

        if let Some(output) = &self.output {
            cmd.arg("-o").arg(output);
        }

        for extra in &self.extras {
            cmd.arg(extra);
        }

        if let Some(index) = &self.indexed_fasta {
            cmd.arg(&index.0 .0);
        }

        match &self.fastq {
            Some(FastQInput::Forward(f)) => {
                cmd.arg(f);
            }
            Some(FastQInput::Reverse(r)) => {
                cmd.arg(r);
            }
            Some(FastQInput::ForwardReverse(f, r)) => {
                cmd.arg(f).arg(r);
            }
            Some(FastQInput::ReverseForward(r, f)) => {
                cmd.arg(r).arg(f);
            }
            None => {
                return None;
            }
        }

        Some(cmd)
    }
}

pub struct BwaMemResult {
    pub sam: PathBuf,
}

impl AlignerResult for BwaMemResult {
    fn alignment_maps(
        &self,
    ) -> impl Iterator<Item = (PathBuf, AlignmentMapFormat, CompressionFormat)> {
        vec![(
            self.sam.clone(),
            AlignmentMapFormat::Sam,
            CompressionFormat::None,
        )]
        .into_iter()
    }

    fn statistics(&self) -> Option<PathBuf> {
        None
    }
}

impl BwaMemBackend {
    pub fn new() -> Self {
        BwaMemBackend { path: None }
    }

    pub fn new_with_path(path: String) -> Self {
        BwaMemBackend { path: Some(path) }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum BwaMemError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Non-zero status: {0}")]
    ExitStatus(ExitStatus),
}

impl From<ExitStatus> for BwaMemError {
    fn from(e: ExitStatus) -> Self {
        BwaMemError::ExitStatus(e)
    }
}

impl Aligner for BwaMemBackend {
    type Params = BwaMemArgs;
    type Result = BwaMemResult;
    type Error = BwaMemError;

    fn name() -> &'static str {
        "bwa mem"
    }

    fn new() -> Self {
        BwaMemBackend { path: None }
    }

    fn new_with_path(path: String) -> Self {
        BwaMemBackend { path: Some(path) }
    }

    fn index(
        &self,
        input: FastaFile,
        _nthreads: usize,
    ) -> Result<<Self::Params as AlignerParams>::GenomeIndex, Self::Error> {
        Command::new(self.path.as_ref().unwrap_or(&"bwa".to_string()))
            .arg("index")
            .arg(&input.0)
            .output()?;

        Ok(BwaMemIndexedFasta(input))
    }

    fn exec(&self, params: &mut Self::Params) -> Result<Self::Result, Self::Error> {
        let mut cmd = Command::new(self.path.as_ref().unwrap_or(&"bwa".to_string()));

        params.commit(&mut cmd);

        let status = cmd.status()?;

        if status.success() {
            Ok(BwaMemResult {
                sam: params
                    .output
                    .clone()
                    .unwrap_or_else(|| PathBuf::from("/dev/stdout")),
            })
        } else {
            Err(status.into())
        }
    }
}
