# Reproduction runs

These are **not** the runs the paper's tables were produced from. The original
numbers were read off the terminal in June/July 2026 and never written to files.
These files are a later re-run, captured verbatim from stdout.

    date        2026-07-29
    machine     Apple M-series laptop, single-threaded
    toolchain   Rust 1.96 (pinned by rust/rust-toolchain.toml)
    build       release

Commands, run from `rust/`:

    cargo run --release -- 16            -> maintenance-2026-07-29.txt
    cargo run --release --bin fk -- 12   -> fk-2026-07-29.txt

Both harnesses run a correctness self-check before timing. The maintenance
harness checks that openings verify, that a tampered opening is rejected, and
that the O(1) update path agrees with a full recommit. The FK harness checks
that the batched proofs are identical to the naive ones and that all of them
verify.

Every column agreed with Tables 1 and 2 to within about 2%, the largest gap
being on the smallest measurement, where timer noise dominates. That level of
agreement is expected here because the machine and the pinned toolchain are the
same; on other hardware the constants will differ and only the asymptotic shape
should carry over.
