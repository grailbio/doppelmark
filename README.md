[![Tests](https://github.com/grailbio/doppelmark/actions/workflows/tests.yml/badge.svg)](https://github.com/grailbio/doppelmark/actions/workflows/tests.yml)
[![Lint](https://github.com/grailbio/doppelmark/actions/workflows/lint.yml/badge.svg)](https://github.com/grailbio/doppelmark/actions/workflows/lint.yml)

## doppelmark duplicate marking tool

doppelmark is a high-performance duplicate sequencing read marking
tool for marking PCR and optical(pad-hopping) duplicate reads. It is
functionally equivalent to the picard and sambamba duplicate marking
tools, but runs much more efficiently and takes advantage of
multi-core hardware. For some workloads and hardware, doppelmark is
100x faster than picard, and 7x faster than sambamba.

doppelmark achieves its speedup by dividing the input into shards and
running the shards in parallel. Each shard includes input
decompression, duplicate marking, and compression of the resulting
output data. It detects duplicates without sorting all records. For a 
detailed description of the algorithm and design,
see [doc.go](https://github.com/grailbio/doppelmark/markduplicates/doc.go).

- [doppelmark](https://github.com/grailbio/doppelmark): High-performance duplicate marking tool

