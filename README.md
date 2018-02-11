# Current Flow Closeness Maximization

Code for paper [https://arxiv.org/abs/1802.02556](https://arxiv.org/abs/1802.02556)

## Requirements

- Julia version 0.6.0

- Julia package [Laplacians.jl](https://github.com/danspielman/Laplacians.jl)
which can be installed by `julia -e 'Pkg.add("Laplacians.jl")'`

## How to run

Execute command `OPENBLAS_NUM_THREADS=1 julia -O3 test.jl data-dir algos k`, where

- `data-dir` is the directory where the edges list files are.
    The algorithms will run on these files in lexicographical order of
    the file names. The result will be printed to both console and file `data-dir.txt`.

- `algos` denotes which algorithms to run.
    `exact` means running exact greedy,
    `approx` means running approx greedy,
    and `both` means running both greedy algorithms (default `both`).

- `k` is an integer, denotes the number of vertices to chose (default `10`).

- `OPENBLAS_NUM_THREADS=1` is to enforce the program to run on a single thread.

For example, `OPENBLAS_NUM_THREADS=1 julia -O3 test.jl data exact 5`
means running the exact greedy on networks in directory `data/` and
choosing 5 vertices for each network. The result will be printed to
both console and `data.txt`.
