import unittest
from unittest.mock import patch, MagicMock
import sys
import os

# Adjust path to import parse_metrics
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import parse_metrics

class TestSamtoolsAA(unittest.TestCase):
    
    @patch("subprocess.run")
    def test_calculate_coverage_aa_flag(self, mock_run):
        # Setup mock subprocess outputs
        # 1. samtools view mapped reads count (return 0)
        mock_view_result = MagicMock()
        mock_view_result.stdout = "0\n"
        
        # 2. samtools depth -aa output (empty or zero depth lines)
        mock_depth_result = MagicMock()
        mock_depth_result.stdout = "ref_contig\t1\t0\nref_contig\t2\t0\nref_contig\t3\t0\n"
        
        mock_run.side_effect = [mock_view_result, mock_depth_result]
        
        ref_lengths = {"ref_contig": 3}
        bam_path = "mock.bam"
        
        results, depths = parse_metrics.calculate_coverage(bam_path, ref_lengths, num_bins=5, min_cov=10)
        
        # Verify that subprocess was called with -aa
        args_called = [call[0][0] for call in mock_run.call_args_list]
        self.assertIn(["samtools", "depth", "-aa", bam_path], args_called)
        
        # Verify depths list
        self.assertEqual(depths["ref_contig"], [0, 0, 0])
        
        # Verify denominator and statistics
        self.assertEqual(results["ref_contig"]["mean_depth"], 0.0)
        self.assertEqual(results["ref_contig"]["coverage_breadth_10x"], 0.0)
        self.assertEqual(results["ref_contig"]["coverage_breadth_min_cov"], 0.0)
        self.assertEqual(results["ref_contig"]["mapped_reads"], 0)

if __name__ == "__main__":
    unittest.main()
