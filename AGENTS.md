# AGENTS.md

VCGlib is a header-only, templated C++ library for triangle/polygon mesh processing. It is the core of MeshLab, so API changes can potentially break downstream code.

## Layout
- `vcg/` the library (`complex/algorithms`, `simplex`, `space`, `math`, `container`, `connectors`)
- `wrap/` I/O and glue code (ply, obj, gl, qt, ...)
- `apps/sample/` small examples, also used as the CI build test
- `eigenlib/` bundled Eigen: do not modify

## Build check
```bash
mkdir -p build && cd build && cmake -GNinja -DVCG_BUILD_EXAMPLES=ON .. && ninja
```

## Library idioms (see `docs/Doxygen/*.dox`)
- Add/delete elements only via `tri::Allocator`; never `push_back`/`resize` mesh containers or call `SetD()`. Use a `PointerUpdater` if you keep element pointers.
- Deleted elements stay in the vectors: skip `IsD()`, never loop up to `VN()`/`FN()`, or compact first (`Allocator::CompactEveryVector`). Prefer `ForEachVertex`/`ForEachFace` (`vcg/complex/foreach.h`).
- Copy meshes with `tri::Append`, never by assignment.
- Check components with `tri::RequireXXX(m)` / `tri::HasXXX(m)` (they handle optional Ocf components).
- Check runtime properties (initialized adjacency, manifoldness, only triangles, ...) with `tri::MeshAssert<MeshType>` (`vcg/complex/algorithms/mesh_assert.h`): it throws `MissingPreconditionException`. Throw on bad input, never `assert`; add a check there rather than writing your own.
- Report progress and diagnostics through a `vcg::CallBackPos *cb` (`wrap/callback.h`), defaulting to null, never with `printf`: callers such as MeshLab route it to a progress bar and a log, while stdout is invisible to them. Long-running algorithms should take one.
- The `V` (visited) bit is scratch: clear it with `UpdateFlags` before use; use `NewBitFlag()` for private bits.
- Algorithms are static members of a class templated on the mesh type.

## Rules
- **Less code is better.** Prefer the solution that needs the fewest lines. Reuse what already exists in vcglib before writing anything new.
- **No duplication** of code or functionality. If something similar exists inside the library, extend or reuse it.
- **Warn before large changes.** If a request needs a lot of new code, say so before starting and propose smaller alternatives.
- **Look for leftovers.** After each change, check for files, functions, includes, or typedefs that are no longer used, and ask before removing them.
- Match the style of the surrounding code (naming, indentation, comment density).
- C++17 is the standard (see `CMakeLists.txt`).
- Samples in `apps/sample` are documentation, not tests (CI only builds them): short, readable top to bottom, showing only the algorithm they are about. Leave edge cases, regression checks and unrelated library features out.
- No preprocessor macros for code, including test helpers like `CHECK(cond)`: write a function (or a lambda) instead.

## Working agreements

**Git is the maintainer's.** Never run a state-changing git command — commit, checkout,
branch, merge, `submodule update`, gitlink bump. Prepare the tree, verify it, say plainly
what changed, and hand over the exact command. Reading state (`status`, `log`, `diff`,
`show`, `grep`) is expected. **Never `git stash`**: the maintainer edits and commits in
this same worktree while a session runs, so a stash sweeps their in-flight work into your
entry and the pop then conflicts on files you never touched. To compare against a
baseline use `git show HEAD:<path>` or a separate `git worktree`.