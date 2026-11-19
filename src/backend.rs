//! Where a read puts its result: the CPU (the default) or, with the Cargo feature `cubecl-cuda`,
//! an NVIDIA GPU. **Beta**: the GPU API may change in a minor release.
//!
//! [Reads](crate::ReadOptionsBuilder#method.read) use the [`Cpu`] backend by default and return an
//! `ndarray` array. With `cubecl-cuda`, add [`.cubecl()`](crate::ReadOptionsBuilder#method.cubecl)
//! before [`read`](crate::ReadOptionsBuilder#method.read-1) to use the
//! [`CubeCl`](struct.CubeCl.html) backend, which returns a
//! [`CubeClMatrix`](struct.CubeClMatrix.html) in GPU memory. Enabling the feature adds that
//! backend without changing any CPU type.
//!
//! This feature builds on Linux x86-64 and ARM64, Windows x86-64, and ARM64 macOS.
//! `CubeCL` 0.11 downloads a prebuilt LLVM bundle that is unavailable for Windows ARM64
//! and Intel macOS, so `cubecl-cuda` (and `--all-features`) cannot build there.
//! ARM64 macOS can compile it but cannot run CUDA. Python GPU wheels target only
//! Linux x86-64 and Windows x86-64. Building needs no CUDA toolkit or `nvcc` because cudarc
//! uses `fallback-latest`; running GPU reads still needs the driver and toolkit.
//!
//! # Reading into GPU memory
//!
//! GPU reads need an NVIDIA GPU, its driver, and a CUDA toolkit. Add the beta with the feature:
//!
//! ```bash
//! cargo add bed-reader@1.3.0-beta.1 --features cubecl-cuda
//! ```
//!
//! ```rust,no_run
//! # #[cfg(feature = "cubecl-cuda")] {
//! use bed_reader::{Bed, ReadOptions, sample_bed_file};
//!
//! let mut bed = Bed::new(sample_bed_file("small.bed")?)?;
//! // A `CubeClMatrix<f32>`, which is always in GPU memory.
//! let cubecl_matrix = ReadOptions::builder().f32().cubecl().read(&mut bed)?;
//! assert_eq!(cubecl_matrix.dim(), (3, 4));
//! # }
//! # use bed_reader::BedErrorPlus;
//! # Ok::<(), Box<BedErrorPlus>>(())
//! ```
//!
//! GPU reads use device 0 and read local files only. `num_threads` and the cloud options do not
//! apply, and `read_and_fill` stays CPU-only. A result has at most `u32::MAX` elements. For further
//! GPU computation, [`tensor_handle`](struct.CubeClMatrix.html#method.tensor_handle) and
//! [`client`](struct.CubeClMatrix.html#method.client) give `CubeCL`'s own types. The [`CubeClMatrix`](struct.CubeClMatrix.html) docs have a fuller example.
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
//! When you drop a `CubeClMatrix`, its GPU memory is freed at the next GPU read. To free it right
//! away, for example so that another GPU library can use it, call [`CubeCl::memory_cleanup`](struct.CubeCl.html#method.memory_cleanup).

use core::fmt::Debug;

use ndarray as nd;

use crate::{BedErrorPlus, BedVal};

// The GPU backend exists only with a `CubeCL` Cargo feature. Everything in this file exists in
// every build, so enabling a GPU feature adds items without changing any existing type.
#[cfg(feature = "cubecl")]
pub(crate) mod cubecl;
#[cfg(feature = "cubecl")]
pub use cubecl::{CubeCl, CubeClMatrix, CubeClVal, CudaPointer};

/// GPU errors.
///
/// This enum is non-exhaustive so additional backend and runtime failures can be
/// represented in minor releases. Downstream matches must include a fallback arm.
#[derive(thiserror::Error, Debug)]
#[non_exhaustive]
pub enum Error {
    /// No usable GPU, or this build does not include its runtime.
    #[cfg(feature = "cubecl")]
    #[error(transparent)]
    Unavailable(#[from] ::cubecl::device::DeviceUnavailable),

    /// A GPU runtime error.
    #[cfg(feature = "cubecl")]
    #[error(transparent)]
    Server(#[from] ::cubecl_runtime::server::ServerError),

    /// Centered floating-point statistics overflowed or were otherwise invalid.
    #[error("Invalid GPU variance: {0}")]
    InvalidVariance(f64),

    /// The GPU code currently indexes with `u32`.
    #[error("Too large for the GPU code ({0} exceeds u32::MAX)")]
    TooLarge(u64),

    /// A CUDA driver error from the pinned-memory upload path.
    #[cfg(feature = "cubecl-cuda")]
    #[error(transparent)]
    Cuda(#[from] cudarc::driver::DriverError),
}

impl From<Error> for Box<BedErrorPlus> {
    fn from(error: Error) -> Self {
        Box::new(BedErrorPlus::GpuError(error))
    }
}

mod sealed {
    /// Keeps [`Backend`](super::Backend) to the backends this crate defines.
    pub trait Sealed {}

    /// Restricts `CubeCL` element types to the types supported by this crate.
    #[cfg(feature = "cubecl")]
    pub trait Value {}

    #[cfg(feature = "cubecl")]
    impl Value for i8 {}
    #[cfg(feature = "cubecl")]
    impl Value for f32 {}
    #[cfg(feature = "cubecl")]
    impl Value for f64 {}
}

/// Where a read puts its result: [`Cpu`] or, with the Cargo feature `cubecl-cuda`, `CubeCl`.
///
/// This trait is sealed: downstream crates cannot implement it. bed-reader can
/// add backend types, builder methods, and opt-in Cargo features in minor releases
/// while preserving existing backends' entry points and result types.
///
/// # Example
///
/// See [`ReadOptionsBuilder::cpu`](crate::ReadOptionsBuilder::cpu) and the [module docs](self).
pub trait Backend: sealed::Sealed + Clone + Copy + Debug + Default + Send + Sync + 'static {
    /// What `read` returns for value type `TVal` with this backend.
    type Matrix<TVal: BedVal>;
}

/// The CPU, the default [`Backend`]: `read` returns an `ndarray` array.
///
/// # Example
///
/// See [`ReadOptionsBuilder::cpu`](crate::ReadOptionsBuilder::cpu).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Cpu;

impl sealed::Sealed for Cpu {}

impl Backend for Cpu {
    type Matrix<TVal: BedVal> = nd::Array2<TVal>;
}
