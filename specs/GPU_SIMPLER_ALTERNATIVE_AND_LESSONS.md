<!-- TODO0 consider deleting this spec once the work below is implemented and released. -->

# GPU support: a simpler alternative, and what the Rust version teaches

This branch adds GPU reads and standardization to bed-reader with CubeCL. It
was left out of the CPU-only 1.3.0 because of its size and complexity. This
note records a much smaller design that would give Python users most of the
same benefit, and the Rust lessons that make this branch worth keeping as
teaching material even if bed-reader ships the smaller design.

## Where the size comes from

Excluding `Cargo.lock` (about 4,000 changed lines), the branch adds about 2,800
lines of Rust and Python, about 2,000 of them code rather than comments or
blank lines. Approximate code-line counts:

| Part | Code lines | Essential? |
|---|---|---|
| Decode kernel and `read_cubecl` orchestration | ~150 | Yes: this is the real work |
| Pinned, double-buffered upload over raw `cudarc` | ~150 | No: a speed optimization |
| GPU standardization (plane and tree reductions, fractional path, apply) | ~550 | Mostly no: CuPy array operations can do it |
| Rust public API (`CubeClMatrix`, `CudaPointer`, `tensor_handle`, errors) | ~150 | Only for Rust GPU users |
| CubeCL-owned memory exposed to CuPy (`GpuArray`, owner tracking, `empty_gpu`, device checks, float64 warning) | ~550 | No: needed only because CubeCL owns the memory |

Most of the size follows from one choice: putting the GPU path inside the
compiled Rust crate. That brings a second memory owner beside CuPy, `unsafe`
CUDA calls, a large `Cargo.lock` change, and separate GPU wheels with their own
CI and qualification (see `GPU_WHEEL_QUALIFICATION.md`).

## The simpler alternative

Keep bed-reader's compiled crate CPU-only and do the GPU work in Python with
CuPy:

1. CuPy allocates a pinned host buffer for the packed bytes
   (`cupyx.empty_pinned`). It is an ordinary NumPy array in page-locked RAM.
2. A new bed-reader CPU function reads only the packed bytes for the selected
   SNPs (the byte span covering the selected individuals in each SNP) into
   that buffer. To Rust it is just a writable `uint8` array: no CUDA, no GPU
   code, and the same index and file-length checks the CPU reader uses.
3. CuPy copies the buffer to the GPU (`cupy.asarray`). Because the buffer is
   pinned, this is a direct DMA copy.
4. A small `cupy.RawKernel` or `ElementwiseKernel` unpacks the 2-bit codes
   through a four-entry lookup table (A1 or A2 counted, missing value) into a
   CuPy output array, handling individual selection, C or F order, and the
   mode-0 transpose.

Standardization then uses ordinary CuPy array operations (about 30 lines for
unit and beta).

Steps 1 and 3 can begin without pinning: a plain NumPy buffer works, with a
slower copy. Measure first and add pinning if the upload matters.

Estimated size: about 300 to 500 lines including tests (an estimate, not yet
prototyped). CuPy becomes an optional extra; the published wheels stay
CPU-only.

### What the alternative gives up

- GPU reads from Rust, and decoding directly into CubeCL memory for a Rust
  CubeCL pipeline.
- The branch's custom pinned upload, until pinning is added through CuPy.

### About the "4x" upload

The branch's `upload_pinned_cuda` copies packed bytes through two reused 16 MiB
pinned buffers, so filling one overlaps the DMA of the other. A code comment
says this is about 4x faster than CubeCL's own upload under WSL2. No benchmark
is recorded, and native Linux was not compared. The upload moves packed data,
1/16 the size of a float32 result, so it may not dominate the whole read. In
the alternative, reading the file directly into a pinned NumPy buffer gets the
same DMA path and also skips the branch's extra copy from a `Vec` into the
pinned buffer.

### What to reuse from this branch

- The `u32` size checks on dimensions, element counts, and packed length.
- Error messages and file-length checks that match the CPU reader.
- The handling of constant fractional SNPs in standardization.
- The tests that compare GPU results with CPU results across files, indexes,
  dtypes, and memory orders.

### Suggested next step

Prototype the CuPy version and benchmark full reads against this branch on the
synthetic 14,925 x 100,000 file, on both WSL2 and native Linux, before deciding
whether any of the CubeCL code belongs in bed-reader.

## Rust lessons from this branch

The smaller design is the better fit for bed-reader, but it does not make an
interesting Rust lesson. This branch does. It could become a separate,
experimental teaching crate (for example `bed-reader-cubecl`) built on the
packed-read function, without bed-reader's compatibility obligations.

1. **GPU kernels written in Rust.** `decode` in `src/gpu/reader.rs` is a
   `#[cube]` function of about 40 lines: a 2-D grid, coalesced writes along the
   contiguous axis, and a strided loop past the grid's y-axis limit.
2. **Compile-time specialization.** `#[comptime] is_f: bool` compiles one
   kernel per memory order, so the hot loop has no runtime branch on layout.
3. **Plane versus tree reductions.** `snp_sums` in
   `src/gpu/reader/standardize.rs` reduces per-SNP counts and sums with plane
   (warp) operations when the device supports them and with a shared-memory
   tree otherwise, and a test checks that both agree.
4. **`unsafe` with stated invariants.** Each `launch_unchecked` carries a
   `// SAFETY:` comment tying every index to the `u32` size checks that
   constructors enforce, a concrete example of making unchecked code
   reviewable.
5. **When `unsafe` earns its place.** The pinned, double-buffered upload uses
   raw CUDA driver calls, events, and a process-wide mutex. It shows how to
   order work across an external stream and a library's own stream, and why to
   wait for queued copies even on error. It also shows the cost: the claim
   needs a benchmark to justify the code.
6. **Testing GPU code without a GPU.** The internal `_cubecl-cpu-test` feature
   runs the same kernels on CubeCL's CPU runtime in CI. The cube size had to
   drop from 256 units to 4, because the CPU runtime spins one thread per unit
   (256 hung; 4 ran in 0.05 s), and small cubes also exercise the strided
   loops.
7. **Typed, sealed APIs over a generic runtime.** The sealed `CubeClVal` trait
   limits element types to `i8`, `f32`, and `f64`, keeping kernel and layout
   assumptions consistent.
8. **Rust and Python sharing GPU memory.** The PyO3 binding exposes
   `__cuda_array_interface__`, so CuPy can use CubeCL-owned memory without a
   copy, and drop order waits for outstanding CUDA work before memory is
   recycled.
9. **The design decision itself.** Comparing about 2,800 lines of GPU Rust with
   a few hundred lines of Python, measured on the same data, is a lesson in
   choosing where complexity belongs.
