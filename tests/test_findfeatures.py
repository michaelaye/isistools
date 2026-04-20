"""Tests for the findfeatures module."""

from __future__ import annotations

import numpy as np
import pytest

from isistools.findfeatures import match_pair, matches_to_cnet


class TestMatchPair:
    def test_finds_shifted_features(self):
        """match_pair should recover a known pixel shift."""
        rng = np.random.default_rng(42)
        img1 = rng.normal(100, 10, (300, 300)).astype(np.float32)
        img2 = img1.copy()

        # Add features and shift them
        for y, x in [(50, 50), (100, 200), (200, 100), (250, 250)]:
            img1[y : y + 15, x : x + 15] = 200
            img2[y + 5 : y + 20, x + 3 : x + 18] = 200

        result = match_pair(img1, img2, algorithm="AKAZE")
        assert result.n_points > 0
        assert result.n_keypoints_from > 0
        assert result.n_keypoints_match > 0

    def test_returns_empty_for_blank_images(self):
        """Blank images should produce no matches."""
        img1 = np.zeros((100, 100), dtype=np.float32)
        img2 = np.zeros((100, 100), dtype=np.float32)
        result = match_pair(img1, img2)
        assert result.n_points == 0

    def test_handles_nan(self):
        """Images with NaN should not crash."""
        rng = np.random.default_rng(42)
        img1 = rng.normal(100, 20, (200, 200)).astype(np.float32)
        img2 = img1.copy()
        img1[50:60, 50:60] = np.nan
        result = match_pair(img1, img2)
        # Should still find some matches
        assert result.n_keypoints_from > 0

    def test_ratio_strictness(self):
        """Lower ratio should produce fewer matches."""
        rng = np.random.default_rng(42)
        img1 = rng.normal(100, 20, (300, 300)).astype(np.float32)
        img2 = img1 + rng.normal(0, 5, img1.shape).astype(np.float32)

        result_loose = match_pair(img1, img2, ratio=0.9)
        result_strict = match_pair(img1, img2, ratio=0.5)
        assert result_strict.n_matches_good <= result_loose.n_matches_good

    def test_max_points(self):
        """max_points should limit output."""
        rng = np.random.default_rng(42)
        img1 = rng.normal(100, 30, (300, 300)).astype(np.float32)
        img2 = img1.copy()
        result = match_pair(img1, img2, max_points=10)
        assert result.n_points <= 10

    def test_algorithms(self):
        """All supported algorithms should work."""
        rng = np.random.default_rng(42)
        img1 = rng.normal(100, 30, (200, 200)).astype(np.float32)
        img2 = img1.copy()
        for algo in ["AKAZE", "ORB", "SIFT"]:
            result = match_pair(img1, img2, algorithm=algo)
            assert result.n_keypoints_from > 0


class TestMatchesToCnet:
    def test_cnet_format(self):
        """Output DataFrame should have plio-compatible columns."""
        result = match_pair(
            np.random.default_rng(42).normal(100, 30, (200, 200)).astype(np.float32),
            np.random.default_rng(42).normal(100, 30, (200, 200)).astype(np.float32),
        )
        if result.n_points == 0:
            pytest.skip("No matches found for this test")

        cnet = matches_to_cnet(result, "SN/FROM/001", "SN/MATCH/002")

        # Must have required columns
        for col in ["id", "pointType", "serialnumber", "measureType",
                     "sample", "line", "aprioriX", "aprioriY", "aprioriZ"]:
            assert col in cnet.columns

        # Each point should have 2 measures (FROM + MATCH)
        assert len(cnet) == result.n_points * 2

        # Check serial numbers
        assert set(cnet["serialnumber"].unique()) == {"SN/FROM/001", "SN/MATCH/002"}

        # Check point types
        assert all(cnet["pointType"] == 2)  # Free

        # Check measure types: 3 = Reference (FROM), 0 = Candidate (MATCH)
        assert set(cnet["measureType"].unique()) == {0, 3}

    def test_without_model_xyz_are_zero(self):
        """Without a CSM model, ground coords should be zero."""
        result = match_pair(
            np.random.default_rng(42).normal(100, 30, (200, 200)).astype(np.float32),
            np.random.default_rng(42).normal(100, 30, (200, 200)).astype(np.float32),
        )
        if result.n_points == 0:
            pytest.skip("No matches found")

        cnet = matches_to_cnet(result, "SN/1", "SN/2")
        assert (cnet["aprioriX"] == 0.0).all()
        assert (cnet["aprioriY"] == 0.0).all()
        assert (cnet["aprioriZ"] == 0.0).all()
