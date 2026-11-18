# What the tests share

Part of the architecture notes of pkpdsim, see [CLAUDE.md](../../CLAUDE.md).

The tests share the production build cache: `tests/conftest.py` points `PKPDSIM_BUILD_CACHE` at a directory shared by the workers of the session (`flat_model_cache`). Builds reuse validated comp/flat pairs by their content, software versions and validation settings, without a test-only flattening wrapper. A test which checks validation or flattening is marked `uncached_build`, which sets `PKPDSIM_BUILD_CACHE=off` for that test.

The tests of these modules run on the toy package `tests/pkpdsim_toy`; its SBML is created once per session by the controller (pytest-xdist), and a test writes only into `tmp_path`.

The PBPK tests compile their models at `O0` (`JIT_OPTIMIZATION` of `tests/pbpk/conftest.py`, the helpers `load`, `compiled`, `individual` and an autouse fixture which sets the default level of sbmlsim, not its environment variable, which a forkserver would keep); the compilation at `O2` was most of the time of the suite. Several individuals of one model come from `individuals(path)` (one loaded model, initialized for one individual after the other), not from a load per individual; `test_load_time` measures at `O2`, the snippets of the docs run at `O0`. The CI runs the tests of every platform in four jobs, `--shard I/4` (`tests/sharding.py`, by a hash of the test id).
