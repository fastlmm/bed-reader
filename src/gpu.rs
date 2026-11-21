//! GPU reads and memory management with the Cargo feature `cubecl-cuda`. **Beta**: the GPU API may change in a minor release.
//!
//! Each CPU read of a local file has a GPU counterpart that returns a [`CubeClMatrix`] in GPU
//! memory instead of an `ndarray` array. They take the same options:
//!
//! | CPU read | GPU read |
//! | -------- | -------- |
//! | [`Bed::read`](crate::Bed::read) | [`Bed::read_cubecl`](crate::Bed::read_cubecl) |
//! | [`ReadOptionsBuilder::read`](crate::ReadOptionsBuilder::read) | [`ReadOptionsBuilder::read_cubecl`](crate::ReadOptionsBuilder::read_cubecl) |
//! | [`Bed::read_with_options`](crate::Bed::read_with_options) | [`Bed::read_cubecl_with_options`](crate::Bed::read_cubecl_with_options) |
//!
//! GPU calls return errors as [`Box<gpu::Error>`](Error); its [`Bed`](Error::Bed) variant holds
//! ordinary read errors, such as a bad index.
//!
//! GPU support needs Rust 1.95; CPU-only builds need Rust 1.85.
//!
//! GPU reads run only on NVIDIA GPUs with CUDA: on Linux (x86-64 and ARM64) and Windows
//! x86-64.
//!
//! Compiling is a separate question. The feature also compiles on ARM64 macOS, so code that
//! enables it still builds there, but macOS has no CUDA, so GPU reads cannot run. The feature
//! (and `--all-features`) cannot compile on Windows ARM64 or Intel macOS, because `CubeCL` 0.11
//! downloads a prebuilt LLVM bundle that does not exist for them. Compiling needs no CUDA
//! installation; running GPU reads needs the NVIDIA driver (see [Requirements](#requirements)).
//!
//! If a CUDA toolkit is installed, compiling does check it: when `nvcc` on `PATH` is newer than the
//! `cudarc` crate supports, the build fails with "Unsupported cuda toolkit version". Run
//! `cargo update -p cudarc` to get a newer `cudarc`, or set the environment variable
//! `CUDARC_CUDA_VERSION` to a supported version, such as `13030` for CUDA 13.3.
//!
//! # Requirements
//!
//! Running GPU reads needs only an NVIDIA GPU and its driver, which provides the CUDA driver
//! library (`libcuda.so` on Linux, `nvcuda.dll` on Windows). No CUDA Toolkit is needed: `CubeCL`
//! compiles bed-reader's GPU kernels itself, with its bundled LLVM, the first time they run.
//! Without a usable GPU or driver, GPU calls return [`Error::Unavailable`].
//!
//! Tested on WSL2 with an NVIDIA RTX 3060 and driver 610.88.
//!
//! # Reading into GPU memory
//!
//! Add the beta with the feature:
//!
//! ```bash
//! cargo add bed-reader@1.3.0-beta.1 --features cubecl-cuda
//! ```
//!
//! ```rust,no_run
//! use bed_reader::{Bed, ReadOptions, sample_bed_file};
//!
//! let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
//!
//! // All data, as a `CubeClMatrix<f32>`, which is always in GPU memory.
//! let all = bed.read_cubecl::<f32>()?;
//! assert_eq!(all.dim(), (3, 4));
//!
//! // Select and choose options with the builder.
//! let some = ReadOptions::builder().sid_index(1..3).f32().read_cubecl(&mut bed)?;
//! assert_eq!(some.dim(), (3, 2));
//!
//! // Or build the options once and reuse them.
//! let read_options = ReadOptions::builder().sid_index(1..3).f32().build()?;
//! let again = bed.read_cubecl_with_options(&read_options)?;
//! assert_eq!(again.dim(), (3, 2));
//! # use bed_reader::gpu::Error;
//! # Ok::<(), Box<Error>>(())
//! ```
//!
//! GPU reads use device 0 and read local files only. On a machine with several GPUs, choose one
//! with the `CUDA_VISIBLE_DEVICES` environment variable, or run one process per GPU.
//! `num_threads` and the cloud options do not apply, and `read_and_fill` is CPU-only.
//!
//! The GPU code indexes with `u32`, so a read fails with [`Error::TooLarge`] if the result has
//! more than `u32::MAX` elements. The bytes read from the file must also fit: selecting a few
//! individuals far apart, across many SNPs, can fail even when the result is small.
//!
//! For further GPU computation, [`CubeClMatrix::tensor_handle`] and [`CubeClMatrix::client`] give
//! `CubeCL`'s own types. The [`CubeClMatrix`] docs have a fuller example.
//!
//! # When to read into GPU memory
//!
//! A GPU read takes about as long as a well-threaded CPU read. Its advantage is skipping the copy
//! to the GPU. If you will compute on the GPU, reading `f32` data directly is about 2.3 times faster
//! than reading on the CPU and then copying (about 1.5 times for `i8`). If you will compute on the
//! CPU, use a CPU read. (Blocks of 1,000 to 4,000 SNPs for 250,000 individuals, file in the
//! operating system's cache, NVIDIA RTX 3060. From a slow drive, the drive dominates.)
//!
//! # Memory
//!
//! When you drop a `CubeClMatrix`, its dedicated GPU memory returns at `CubeCL`'s next allocation.
//! To free it right away, for example so that another GPU library can use it, call
//! [`memory_cleanup`].

use core::fmt::Debug;
use std::io;

use cubecl::device::DeviceUnavailable;
use cubecl_runtime::server::ServerError;
use cudarc::driver::DriverError;

use crate::{BedError, BedErrorPlus};

pub(crate) mod reader;
pub use reader::{client, memory_cleanup, CubeClMatrix, CubeClVal, CudaPointer};

/// GPU failures and ordinary errors encountered during a GPU read.
///
/// This enum is non-exhaustive so additional runtime failures can be represented in minor
/// releases. Downstream matches must include a fallback arm.
///
/// Conversions into `Box<Error>` follow the crate's boxed read-error pattern. The standard
/// library's blanket conversions do not compose `From<Source> for Error` with boxing; explicit
/// implementations are required, and coherence prevents replacing them with a generic blanket
/// implementation. Ordinary read and I/O errors retain their original sources inside `Bed`.
#[derive(thiserror::Error, Debug)]
#[non_exhaustive]
pub enum Error {
    /// An ordinary read failure, such as an invalid index or malformed file.
    #[error(transparent)]
    Bed(#[from] Box<BedErrorPlus>),

    /// No usable GPU.
    #[error(transparent)]
    Unavailable(#[from] DeviceUnavailable),

    /// A GPU runtime error.
    #[error(transparent)]
    Server(#[from] ServerError),

    /// The GPU code currently indexes with `u32`.
    #[error("Too large for the GPU code ({0} exceeds u32::MAX)")]
    TooLarge(u64),

    /// A CUDA driver error from the pinned-memory upload path.
    #[error(transparent)]
    Cuda(#[from] DriverError),
}

impl From<Box<BedErrorPlus>> for Box<Error> {
    fn from(error: Box<BedErrorPlus>) -> Self {
        Box::new(Error::from(error))
    }
}

impl From<BedError> for Box<Error> {
    fn from(error: BedError) -> Self {
        Box::new(Error::Bed(error.into()))
    }
}

impl From<io::Error> for Box<Error> {
    fn from(error: io::Error) -> Self {
        Box::new(Error::Bed(error.into()))
    }
}

impl From<DeviceUnavailable> for Box<Error> {
    fn from(error: DeviceUnavailable) -> Self {
        Box::new(Error::from(error))
    }
}

impl From<ServerError> for Box<Error> {
    fn from(error: ServerError) -> Self {
        Box::new(Error::from(error))
    }
}

impl From<DriverError> for Box<Error> {
    fn from(error: DriverError) -> Self {
        Box::new(Error::from(error))
    }
}

mod sealed {
    /// Restricts `CubeCL` element types to the types supported by this crate.
    pub trait Value {}

    impl Value for i8 {}
    impl Value for f32 {}
    impl Value for f64 {}
}
