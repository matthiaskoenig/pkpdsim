import pkpdsim


def test_version() -> None:
    assert pkpdsim.__version__ == "0.1.0"


def test_all() -> None:
    assert pkpdsim.__all__ == ["__version__"]
