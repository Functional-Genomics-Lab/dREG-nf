use std::{
    path::PathBuf,
    process::{Command, ExitStatus},
};

use crate::{
    Aligner, AlignerParams, AlignerResult, AlignmentMapFormat, CompressionFormat, FastQInput,
    FastaFile,
};

pub struct Bowtie2Backend {
    path: Option<String>,
}

#[derive(Debug, Clone, Copy)]
pub enum Preset {
    VeryFast,
    Fast,
    Sensitive,
    VerySensitive,
}

pub struct Bowtie2IndexedFasta(pub FastaFile);

pub struct Bowtie2Args {
    fastq: Option<FastQInput>,
    nthreads: usize,
    indexed_fasta: Option<Bowtie2IndexedFasta>,
    preset: Preset,
    allow_spliced_alignments: bool,
    extras: Vec<String>,
    output: Option<PathBuf>,
    local: bool,
    memory_map: bool,
}

impl Bowtie2Args {
    pub fn append_extra(&mut self, extra: String) -> &mut Self {
        self.extras.push(extra);
        self
    }

    pub fn preset(&mut self, preset: Preset) -> &mut Self {
        self.preset = preset;
        self
    }

    pub fn memory_map(&mut self, memory_map: bool) -> &mut Self {
        self.memory_map = memory_map;
        self
    }

    pub fn local(&mut self) -> &mut Self {
        self.local = true;
        self
    }

    pub fn end_to_end(&mut self) -> &mut Self {
        self.local = false;
        self
    }
}

impl AlignerParams for Bowtie2Args {
    type GenomeIndex = Bowtie2IndexedFasta;

    fn new() -> Self {
        Bowtie2Args {
            fastq: None,
            nthreads: 1,
            indexed_fasta: None,
            allow_spliced_alignments: true,
            preset: Preset::VerySensitive,
            extras: Vec::new(),
            memory_map: true,
            local: false,
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

    fn genome_index(&mut self, index: Bowtie2IndexedFasta) -> Option<&mut Self> {
        self.indexed_fasta = Some(index);
        Some(self)
    }

    fn allow_spliced_alignments(&mut self, allow: bool) -> Option<&mut Self> {
        self.allow_spliced_alignments = allow;
        Some(self)
    }

    fn commit<'a, 'b: 'a>(&'a self, cmd: &'b mut Command) -> Option<&'b mut Command> {
        cmd.arg("-q");

        if self.local {
            cmd.arg("--local");
        } else {
            cmd.arg("--end-to-end");
        }

        cmd.arg("-p").arg(self.nthreads.to_string());

        if self.memory_map {
            cmd.arg("--mm");
        }

        // seem to not mean what I think this means :(
        /*
        if !self.allow_spliced_alignments {
            cmd.arg("--rdg")
                .arg("10000,10000")
                .arg("--rfg")
                .arg("10000,10000");
        }
        */

        if let Some(output) = &self.output {
            cmd.arg("-S").arg(output);
        }

        for extra in &self.extras {
            cmd.arg(extra);
        }

        if let Some(index) = &self.indexed_fasta {
            cmd.arg("-x").arg(&index.0 .0);
        } else {
            return None;
        }

        match &self.fastq {
            Some(FastQInput::Forward(f)) => {
                cmd.arg("-U").arg(f);
            }
            Some(FastQInput::Reverse(r)) => {
                cmd.arg("-U").arg(r);
            }
            Some(FastQInput::ForwardReverse(f, r)) => {
                cmd.arg("-1").arg(f).arg("-2").arg(r);
            }
            Some(FastQInput::ReverseForward(r, f)) => {
                cmd.arg("-1").arg(f).arg("-2").arg(r);
            }
            None => {
                return None;
            }
        }

        Some(cmd)
    }
}

pub struct Bowtie2Result {
    pub sam: PathBuf,
}

impl AlignerResult for Bowtie2Result {
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

impl Bowtie2Backend {
    pub fn new() -> Self {
        Bowtie2Backend { path: None }
    }

    pub fn new_with_path(path: String) -> Self {
        Bowtie2Backend { path: Some(path) }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum Bowtie2Error {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Non-zero status: {0}")]
    ExitStatus(ExitStatus),
}

impl From<ExitStatus> for Bowtie2Error {
    fn from(e: ExitStatus) -> Self {
        Bowtie2Error::ExitStatus(e)
    }
}

impl Aligner for Bowtie2Backend {
    type Params = Bowtie2Args;
    type Result = Bowtie2Result;
    type Error = Bowtie2Error;

    fn new() -> Self {
        Bowtie2Backend { path: None }
    }

    fn new_with_path(path: String) -> Self {
        Bowtie2Backend { path: Some(path) }
    }

    fn name() -> &'static str {
        "bowtie2"
    }

    fn index(
        &self,
        input: FastaFile,
        nthreads: usize,
    ) -> Result<<Self::Params as AlignerParams>::GenomeIndex, Self::Error> {
        Command::new(self.path.as_ref().unwrap_or(&"bowtie2-build".to_string()))
            .arg("-f")
            .arg("--threads")
            .arg(nthreads.to_string())
            .arg(&input.0)
            .arg(&input.0)
            .output()?;

        Ok(Bowtie2IndexedFasta(input))
    }

    fn exec(&self, params: &mut Self::Params) -> Result<Self::Result, Self::Error> {
        let mut cmd = Command::new(self.path.as_ref().unwrap_or(&"bowtie2".to_string()));

        params.commit(&mut cmd);

        let status = cmd.status()?;

        if status.success() {
            Ok(Bowtie2Result {
                sam: params.output.clone().unwrap(),
            })
        } else {
            Err(status.into())
        }
    }
}
