# What the tests share

Part of the architecture notes of pkpdsim, see [CLAUDE.md](../../CLAUDE.md).

The tests share the flat models of the session: `tests/conftest.py` patches `pkpdsim.pbpk.builder.flat_model` for the session (`flat_model_cache`) with a cache by the sha256 of the comp file without the dates of the model history, in a folder of the session which every worker of pytest-xdist shares (about a quarter of the builds of the suite repeat a comp model); a test which checks what the validation logs is marked `uncached_build`.

The tests of these modules run on the toy package `tests/pkpdsim_toy`; its SBML is created once per session by the controller (pytest-xdist), and a test writes only into `tmp_path`.

The PBPK tests compile their models at `O0` (`JIT_OPTIMIZATION` of `tests/pbpk/conftest.py`, the helpers `load`, `compiled`, `individual` and an autouse fixture which sets the default level of sbmlsim, not its environment variable, which a forkserver would keep); the compilation at `O2` was most of the time of the suite. Several individuals of one model come from `individuals(path)` (one loaded model, initialized for one individual after the other), not from a load per individual; `test_load_time` measures at `O2`, the snippets of the docs run at `O0`. The CI runs the tests of every platform in four jobs, `--shard I/4` (`tests/sharding.py`, by a hash of the test id).
