# shocktube

A minimal 3D compressible Euler solver for the classic Sod shock-tube problem. The
simulation is set up on a cubic grid but is effectively one-dimensional: the initial
data depend only on `x`, and the boundary conditions are outflow (first-order
extrapolation) in `x` and periodic in `y` and `z`.

The hydro algorithm is an unsplit Godunov method with PPM reconstruction and an
approximate Riemann solver for an ideal gas with `γ = 1.4`.

## Components

- **FAB** (`FAB.H`) — Fortran-ordered 3D array with ghost cells.
- **Hydro** (`Hydro.H`, `Hydro.cpp`) — problem setup, boundaries, conservative-to-primitive
  conversion, and the Godunov update.
- **main** (`main.cpp`) — time integration loop and I/O.

At the end of a run the code writes a 1D slice of cell-centered primitive data
to `output.txt`: density, velocity, and pressure along `x` at `j = k = 0` only. The
full 3D solution is not written to disk.

## Requirements

- C++17 compiler (tested with `g++`)
- GNU `make`
- Python 3 and `matplotlib` (for plotting only)

## Build

From the repository root:

```bash
make
```

The build rules live in `GNUmakefile`, which GNU `make` reads automatically.

Object files and dependency files are placed under `tmp_build_dir/<variant>/`. The
executable is written to the current directory as `shock.<variant>.ex`, where
`<variant>` is `<build>.<compiler>.<precision>[.san]` — e.g. `opt.gcc.double`.

| Variable    | Default | Effect |
|-------------|---------|--------|
| `DEBUG`      | `FALSE` | `TRUE` → debug build (`-g -O0`, `-DDEBUG` for FAB bounds and NaN init) |
| `PRECISION`  | `DOUBLE`| `FLOAT` → single-precision (`USE_FLOAT`) |
| `FSANITIZE`  | `FALSE` | `TRUE` → enable AddressSanitizer + UBSan (`-fsanitize=address,undefined`) |
| `CXX`        | `g++`   | C++ compiler |

### Examples

```bash
make
# → shock.opt.gcc.double.ex

make DEBUG=TRUE
# → shock.debug.gcc.double.ex

make PRECISION=FLOAT
# → shock.opt.gcc.float.ex

make FSANITIZE=TRUE
# → shock.opt.gcc.double.san.ex

make CXX=clang++
# → shock.opt.clang.double.ex

make CXX=clang++ DEBUG=TRUE FSANITIZE=TRUE
# → shock.debug.clang.double.san.ex
```

### Clean

```bash
make clean      # remove tmp_build_dir/
make realclean  # also remove shock.*.ex
```

## Run

The executable takes exactly one argument: the number of cells in the `x` direction.
The `y` and `z` directions are fixed at 8 cells each.

```bash
./shock.opt.gcc.double.ex 128
```

This runs the Sod problem on a `128×8×8` logical grid (with ghost cells) until
`t = 0.2`, printing time-step information to stdout. The final state is reduced to
a single 1D slice (`i` along `x` at `j = k = 0`) and written to `output.txt`.

## Plot results

`plot_output.py` reads the 1D slice in `output.txt` and plots density, velocity, and
pressure vs. cell-centered `x`.

```bash
python3 plot_output.py
python3 plot_output.py -o sod.png --no-show
python3 plot_output.py other.txt
python3 plot_output.py -h
```

Use `-o FILE` to save a figure without displaying it. Use `--no-show` to skip the
interactive window (useful on headless systems).
