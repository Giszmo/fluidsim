# fluidsim — real-time particle-based fluid simulation

An SPH (smoothed particle hydrodynamics) fluid solver with an OpenGL/GLUT viewer.
The solver is `fluid.h`, a single header; `main.cpp` is the viewer and the scene.

This is the later, reworked descendant of the code written for the diploma thesis
*Partikelbasierte Fluidsimulation für interaktive Anwendungen* (Leo Wandersleb,
TU München, 2005). The thesis, its sources and the 2005 state of the code are in
**[Giszmo/pbfs](https://github.com/Giszmo/pbfs)** — thesis PDF:
<https://giszmo.github.io/pbfs/Text/Diplomarbeit.pdf>.

## The idea

Neighbour search without a grid. Each particle stores an implicit **cell ID** and the
particles are kept **sorted by it**, so nothing is allocated per cell and memory scales
with the particle count rather than with the size of the domain. The cell ID is taken
modulo the ID space, which makes the domain effectively unbounded; the collisions that
introduces are resolved by the distance test SPH performs anyway. Adjacency along x
comes free from the sort order, so 3D needs four sorted lists rather than the eight
shifted grids of the half-cell-offset trick. The thesis reports domains of (100 000)³
cells at 3 000 particles in under 15 MB.

## Building

Needs a C++ compiler, GLUT, OpenGL and Boost.Thread. On Debian/Ubuntu:

```sh
apt install build-essential freeglut3-dev libboost-thread-dev
make            # -> main.$(uname -s).$(major)_std   (debug, -g)
make opt        # -> main.$(uname -s).$(major)_opt   (-O3)
```

Verified on Debian forky, GCC 15.3.0, freeglut 3.4.0, Boost 1.90 — no source changes
needed. `make` also still honours the `dbg` and IRIX branches inherited from the TUM
course makefile.

## Running

```sh
./main.Linux.7_opt
```

The simulation starts **paused**; press `p` to run it. Moving the mouse over the window
rotates the camera. It prints a frames/particles line per second to stdout — about
1 450 fps at 3 760 particles on a 2026 desktop CPU, software-rendered.

| Key | |
|---|---|
| `p` | pause / resume |
| `+` / `-` | zoom in / out |
| `*` / `/` | faster / slower timestep |
| `v` | show/hide particles |
| `c` | show/hide cells |
| `n` | show/hide normals |
| `l` | toggle lighting |
| `f` / `g` | flat / smooth shading |
| `w` / `q` / `e` / Tab | wireframe, points, back-face lines, front-face cull |
| `y` | print the volume of the closed surface |
| `o` | dump vertices, normals and indices to stdout |
| `Esc` | quit |

### The default scene

A `sink`-shaped ground height field, a 4×4 plate at z=0 that sets the speed of any
particle crossing it to +70 (so the fluid is thrown back up), and two tetrahedra of
liquid dropped from above. **All the barrier, teleport and extra set-speed geometry in
`main()` is commented out**, so the default scene contains no walls — worth knowing
before concluding that barrier collision does not work.

## Headless regression check

`test/headless.cpp` runs the solver without a window and prints an FNV-1a checksum over
every particle position, so two builds can be compared bit for bit:

```sh
make check
```

Scenario 0 drops liquid on the centre of a square plate (it comes to rest on it);
scenario 1 drops it beside the plate (it falls straight past). Both are deterministic.

## Fixed here

**Out-of-bounds read in the neighbour scan.** The four scan loops in
`collidewithneighbours` and `collidewithneighboursforrho` read `_sortlist[...]` *before*
the `sortlistindex < _particlecount` bounds check. `&&` evaluates left to right, so every
scan read one slot past the populated region and decided whether to keep walking on
whatever was in it. `_sortlist` was never initialised, so that slot was heap garbage.
The backward loops were worse: `sortlistindex` starts at `startsortlistid - sortlistbase
- 1`, which underflows to 2^64−1 when the particle is first in its list, so the read
landed in the previous list's memory. The `//TODO: workaround??` guard in `condcollide`
(`( first < second ) & ( second < _particlecount )`) exists to throw away the bogus
particle ids this produced — it is why the simulation worked at all.

The fix swaps the two operands of `&&` in all four loops so the bounds check runs first,
and value-initialises `_sortlist`. Valgrind, scenarios 0 and 1, 0.8 s of simulated time:
391 and 383 errors from 4 contexts before, **0 from 0 contexts** after.

It is behaviour-neutral on Linux — checksums are identical across scenarios 0 and 1 at
3 s and 8 s of simulated time — because fresh pages from the kernel are zeroed, so the
garbage reads as zero and the `condcollide` guard drops it. On the MSVC heap the code
was originally written against, it would not have been zero.

**Linking with GCC.** The makefile drives the compiler as `gcc`, which no longer pulls in
`libstdc++`; it is `g++` now. That is the only change needed to build on a current
toolchain.

## Known rough edges

- `fluid.h` is MSVC-flavoured C++ from 2005: `unsigned __int64` throughout, and the whole
  solver in one 2 800-line header. `bicubic_bezier_patch.h` has no include guard.
- The makefile is the 2004 TUM course template (`$Id: Makefile,v 1.1 2004/04/16 ...
  prkipfer $`) with CRLF line endings and `-L` paths for IRIX and XFree86 that have not
  existed for twenty years.
- In the BOUNDARY branch of `collide()`, `float d = r.abs()` is computed and never used;
  the repulsion feeds the signed plane distance into `w_poly6_grad()`, which expects a
  radial distance. Harmless for the boxed barrier patches the code actually builds, but
  it is a loose end.
- Marching cubes renders 0 triangles in the default scene; the surface path is only
  exercised with a different scene setup.
