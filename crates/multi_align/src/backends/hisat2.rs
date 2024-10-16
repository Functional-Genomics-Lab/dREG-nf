use std::{
    io,
    path::PathBuf,
    process::{Command, ExitStatus},
};

use crate::{Aligner, AlignerParams, AlignerResult, FastQInput};

pub struct Hisat2Backend {
    path: Option<String>,
}

pub enum Hisat2Preset {
    Fast,
    Sensitive,
    VerySensitive,
}

pub struct Hisat2Args {
    fastq: Option<FastQInput>,
    nthreads: usize,
    preset: Hisat2Preset,
    no_spliced_alignments: bool,
    extras: Vec<String>,
    output: Option<PathBuf>,
    index: Option<Hisat2Index>,
}

impl Hisat2Args {
    pub fn append_extra(&mut self, extra: String) -> &mut Self {
        self.extras.push(extra);
        self
    }
}

pub struct Hisat2Index(pub PathBuf);

impl AlignerParams for Hisat2Args {
    type GenomeIndex = Hisat2Index;

    fn new() -> Self {
        Hisat2Args {
            fastq: None,
            nthreads: 1,
            preset: Hisat2Preset::VerySensitive,
            no_spliced_alignments: false,
            extras: Vec::new(),
            output: None,
            index: None,
        }
    }

    fn set_output(&mut self, output: PathBuf) -> Option<&mut Self> {
        self.output = Some(output);
        Some(self)
    }

    fn genome_index(&mut self, index: Self::GenomeIndex) -> Option<&mut Self> {
        self.index = Some(index);
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

    fn allow_spliced_alignments(&mut self, allow: bool) -> Option<&mut Self> {
        self.no_spliced_alignments = !allow;
        Some(self)
    }

    fn commit<'a, 'b: 'a>(&'a self, cmd: &'b mut Command) -> Option<&'b mut Command> {
        cmd.arg("--new-summary");

        match self.preset {
            Hisat2Preset::Fast => {
                cmd.arg("--very-fast");
            }
            Hisat2Preset::Sensitive => {
                cmd.arg("--sensitive");
            }
            Hisat2Preset::VerySensitive => {
                cmd.arg("--very-sensitive");
            }
        }

        if self.no_spliced_alignments {
            cmd.arg("--no-spliced-alignment");
        }

        cmd.args(self.extras.iter());

        cmd.arg("-x")
            .arg(&self.index.as_ref()?.0)
            .arg("-p")
            .arg(self.nthreads.to_string());

        Some(match &self.fastq {
            Some(FastQInput::Forward(f)) => cmd.arg("-U").arg(f).arg("--rna-strandness").arg("F"),
            Some(FastQInput::Reverse(r)) => cmd.arg("-U").arg(r).arg("--rna-strandness").arg("R"),
            Some(FastQInput::ForwardReverse(f, r)) => cmd
                .arg("-1")
                .arg(f)
                .arg("-2")
                .arg(r)
                .arg("--rna-strandness")
                .arg("FR"),
            Some(FastQInput::ReverseForward(r, f)) => cmd
                .arg("-1")
                .arg(f)
                .arg("-2")
                .arg(r)
                .arg("--rna-strandness")
                .arg("RF"),

            None => return None,
        })
        .map(|cmd| {
            if self.output.is_some() {
                cmd.arg("-S")
                    .arg(self.output.clone().unwrap())
                    .arg("--summary-file")
                    .arg(self.output.clone().unwrap().with_extension("mapstats"))
            } else {
                cmd
            }
        })
    }
}

pub struct Hisat2Result {
    sam_path: PathBuf,
    mapstats_path: PathBuf,
}

impl AlignerResult for Hisat2Result {
    fn alignment_maps(
        &self,
    ) -> impl Iterator<Item = (PathBuf, crate::AlignmentMapFormat, crate::CompressionFormat)> {
        vec![(
            self.sam_path.clone(),
            crate::AlignmentMapFormat::Sam,
            crate::CompressionFormat::None,
        )]
        .into_iter()
    }

    fn statistics(&self) -> Option<PathBuf> {
        Some(self.mapstats_path.clone())
    }
}

#[derive(Debug, thiserror::Error)]
pub enum Hisat2Error {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("Command exited with non-zero status: {0}")]
    Command(ExitStatus),
}

impl From<ExitStatus> for Hisat2Error {
    fn from(e: ExitStatus) -> Self {
        Hisat2Error::Command(e)
    }
}

impl Aligner for Hisat2Backend {
    type Params = Hisat2Args;
    type Result = Hisat2Result;
    type Error = Hisat2Error;

    fn name() -> &'static str {
        "hisat2"
    }

    fn new() -> Self {
        Hisat2Backend { path: None }
    }

    fn new_with_path(path: String) -> Self {
        Hisat2Backend { path: Some(path) }
    }

    fn exec(&self, params: &mut Self::Params) -> Result<Self::Result, Self::Error> {
        let mut cmd =
            std::process::Command::new(self.path.as_ref().unwrap_or(&"hisat2".to_string()));

        let cmd = params
            .commit(&mut cmd)
            .expect("Incomplete input for hisat2");

        let status = cmd.status()?;

        if status.success() {
            Ok(Hisat2Result {
                sam_path: params
                    .output
                    .clone()
                    .unwrap_or_else(|| PathBuf::from("/dev/stdout")),
                mapstats_path: params
                    .output
                    .clone()
                    .unwrap_or_else(|| PathBuf::from("/dev/stdout"))
                    .with_extension("mapstats"),
            })
        } else {
            Err(status.into())
        }
    }
}
