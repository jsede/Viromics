#!/usr/bin/env python3
import os
import subprocess
import sys
import unittest
from unittest.mock import MagicMock, patch

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import parse_metrics  # noqa: E402


class TestSamtoolsAA(unittest.TestCase):
    @patch("subprocess.run")
    def test_calculate_coverage_synced_flags_and_intervals(self, mock_run):
        mock_view_result = MagicMock()
        mock_view_result.stdout = "12\n"

        mock_depth_result = MagicMock()
        mock_depth_result.stdout = "ref_contig\t1\t0\nref_contig\t2\t50\nref_contig\t3\t0\n"

        mock_run.side_effect = [mock_view_result, mock_depth_result]

        ref_lengths = {"ref_contig": 3}
        bam_path = "mock.bam"

        results, depths = parse_metrics.calculate_coverage(
            bam_path,
            ref_lengths,
            num_bins=5,
            min_cov=10,
            min_consensus_depth=10,
            min_map_qual=20,
            min_base_qual=20,
        )

        args_called = [call[0][0] for call in mock_run.call_args_list]
        self.assertIn(["samtools", "view", "-c", "-F", "2308", "-q", "20", bam_path, "ref_contig"], args_called)
        self.assertIn(
            ["samtools", "depth", "-aa", "-J", "-Q", "20", "-q", "20", "-G", "0xF04", bam_path],
            args_called,
        )

        self.assertEqual(depths["ref_contig"], [0, 50, 0])
        self.assertEqual(results["ref_contig"]["mean_depth"], 16.67)
        self.assertEqual(results["ref_contig"]["coverage_breadth_min_cov"], 33.33)
        self.assertEqual(results["ref_contig"]["mapped_reads"], 12)
        self.assertEqual(results["ref_contig"]["low_coverage_intervals"], [[1, 1], [3, 3]])

    @patch("subprocess.run")
    def test_calculate_coverage_fails_loudly_on_samtools_error(self, mock_run):
        mock_run.side_effect = subprocess.CalledProcessError(1, "samtools")
        with self.assertRaises(subprocess.CalledProcessError):
            parse_metrics.calculate_coverage("broken.bam", {"ref_contig": 10})


if __name__ == "__main__":
    unittest.main()
