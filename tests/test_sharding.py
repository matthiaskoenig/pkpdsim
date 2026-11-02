"""Tests of the shares of the tests (`--shard`, `sharding.py`)."""

import pytest
from sharding import in_shard, parse_shard

NODEIDS = [
    f"tests/pbpk/test_m{m}.py::test_{t}[{p}]"
    for m in range(5)
    for t in range(40)
    for p in "ab"
]


def test_every_test_belongs_to_one_share() -> None:
    for nodeid in NODEIDS:
        assert sum(in_shard(nodeid, index, 4) for index in range(1, 5)) == 1


def test_the_shares_divide_the_tests_evenly() -> None:
    sizes = [sum(in_shard(n, index, 4) for n in NODEIDS) for index in range(1, 5)]
    assert min(sizes) > 0.8 * len(NODEIDS) / 4


def test_one_share_is_every_test() -> None:
    assert all(in_shard(nodeid, 1, 1) for nodeid in NODEIDS)


@pytest.mark.parametrize(("text", "shard"), [("1/1", (1, 1)), ("3/4", (3, 4))])
def test_a_shard_is_parsed(text: str, shard: tuple[int, int]) -> None:
    assert parse_shard(text) == shard


@pytest.mark.parametrize("text", ["0/4", "5/4", "1", "a/b", "1/0", "-1/4", "1/4/2"])
def test_a_malformed_shard_is_refused(text: str) -> None:
    with pytest.raises(ValueError, match="I/N"):
        parse_shard(text)
