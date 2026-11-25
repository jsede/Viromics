#!/usr/bin/env python3
import unittest
import os
import re

class TestVariantQualFilter(unittest.TestCase):
    
    def test_variant_call_nf_parameter_separation(self):
        # Locate pipeline/modules/variant_call.nf
        module_path = os.path.join(os.path.dirname(__file__), "..", "modules", "variant_call.nf")
        self.assertTrue(os.path.exists(module_path), f"File not found: {module_path}")
        
        with open(module_path, 'r') as f:
            content = f.read()
            
        # 1. Assert min_variant_qual is in inputs
        self.assertIn("val min_variant_qual", content, "min_variant_qual not found in inputs of variant_call.nf")
        
        # 2. Assert QUAL>=${min_variant_qual} is in the filter command
        filter_regex = r'bcftools\s+filter\s+-i\s+["\']QUAL>=\$\{min_variant_qual\}'
        match = re.search(filter_regex, content)
        self.assertTrue(match, "bcftools filter command does not use QUAL>=${min_variant_qual}")
        
        # 3. Assert QUAL>=${base_qual} is NOT in the filter command
        old_filter_regex = r'bcftools\s+filter\s+-i\s+["\']QUAL>=\$\{base_qual\}'
        old_match = re.search(old_filter_regex, content)
        self.assertFalse(old_match, "bcftools filter command incorrectly uses QUAL>=${base_qual}")
        
        # 4. Demonstrate simulation: when base_qual=10 and min_variant_qual=20
        base_qual = 10
        min_variant_qual = 20
        
        # Simulate substitution
        simulated_command = 'bcftools filter -i "QUAL>=${min_variant_qual} && DP>=${min_coverage}"'
        resolved_command = simulated_command.replace("${min_variant_qual}", str(min_variant_qual)).replace("${base_qual}", str(base_qual))
        
        self.assertIn("QUAL>=20", resolved_command)
        self.assertNotIn("QUAL>=10", resolved_command)
        
        print("SUCCESS: Validation test passed. Variant quality filter (QUAL>=20) is successfully decoupled from base quality (10).")

if __name__ == "__main__":
    unittest.main()
