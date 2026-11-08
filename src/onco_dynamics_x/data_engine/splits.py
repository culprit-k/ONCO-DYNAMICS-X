"""Patient-grouped data partitions to prevent sample-level information leakage."""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

IndexArray = NDArray[np.int64]


def group_train_test_split(
    patient_ids: list[str] | tuple[str, ...] | NDArray[np.str_],
    test_size: float = 0.2,
    random_seed: int = 42,
) -> tuple[IndexArray, IndexArray]:
    """Split sample indices while keeping every patient's samples in one partition.

    This is a reproducible group split, not a stratified split. If class balance
    matters, perform a group-aware stratified design and report the resulting
    patient/event counts rather than falling back to sample-random splitting.
    """
    raw_groups = np.asarray(patient_ids, dtype=object)
    if raw_groups.ndim != 1 or len(raw_groups) < 2:
        raise ValueError("patient_ids must be a one-dimensional array with at least two rows.")
    for item in raw_groups:
        if item is None or (isinstance(item, (float, np.floating)) and np.isnan(item)):
            raise ValueError("patient_ids must not contain empty values.")
        if not str(item).strip():
            raise ValueError("patient_ids must not contain empty values.")
    groups = np.asarray([str(item) for item in raw_groups], dtype=str)
    if not 0.0 < test_size < 1.0:
        raise ValueError("test_size must be strictly between 0 and 1.")
    unique_groups = np.unique(groups)
    if len(unique_groups) < 2:
        raise ValueError("At least two distinct patients are required for a group split.")

    rng = np.random.default_rng(random_seed)
    shuffled = rng.permutation(unique_groups)
    test_count = min(
        len(unique_groups) - 1,
        max(1, int(np.ceil(len(unique_groups) * test_size))),
    )
    test_groups = set(shuffled[:test_count].tolist())
    test_mask = np.fromiter(
        (group in test_groups for group in groups),
        dtype=bool,
        count=len(groups),
    )
    train_indices = np.flatnonzero(~test_mask).astype(np.int64)
    test_indices = np.flatnonzero(test_mask).astype(np.int64)
    if set(groups[train_indices]) & set(groups[test_indices]):
        raise RuntimeError("Internal error: patient groups overlap between partitions.")
    return train_indices, test_indices
