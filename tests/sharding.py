"""The shares of the tests of `--shard I/N`, one job of the CI matrix each.

A test belongs to the share of a hash of its id (CRC-32, the same on every
platform and in every process, so the workers of pytest-xdist agree), not of
its module: the PBPK tests of one module take minutes, the tests one by one
divide the suite into shares within a few percent of each other.
"""

import zlib


def parse_shard(text: str) -> tuple[int, int]:
    """Parse `I/N`, the share `I` of `N` shares.

    Args:
        text: the value of `--shard`, e.g. `2/4`.

    Returns:
        The share and the number of shares.

    Raises:
        ValueError: if the text is not `I/N` with `1 <= I <= N`.
    """
    index, _, count = text.partition("/")
    if not (index.isdigit() and count.isdigit() and 1 <= int(index) <= int(count)):
        raise ValueError(f"--shard takes I/N with 1 <= I <= N, not '{text}'.")
    return int(index), int(count)


def in_shard(nodeid: str, index: int, count: int) -> bool:
    """Tell whether a test belongs to the share `index` of `count`.

    Args:
        nodeid: the id of the test, e.g. `tests/test_pk.py::test_auc`.
        index: the share, from 1 to `count`.
        count: the number of shares.

    Returns:
        Whether the test belongs to the share.
    """
    return zlib.crc32(nodeid.encode()) % count == index - 1
