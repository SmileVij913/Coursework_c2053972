import unittest
from glycolysis_citric_cycle_simulation import (
    glycolysis_reaction1,
    glycolysis_reaction2,
    glycolysis_reaction3,
    citric_acid_cycle_reaction1,
    citric_acid_cycle_reaction2,
)


class TestGlycolysisCitricCycleSimulation(unittest.TestCase):

    def test_glycolysis_reaction1(self):
        # Test with known inputs and outputs
        glucose = 10.0
        atp = 5.0
        expected_output = 0.8264462809917354  # Corrected expected output
        self.assertAlmostEqual(glycolysis_reaction1(glucose, atp), expected_output, places=5)

    def test_glycolysis_reaction2(self):
        # Test with known inputs and outputs
        glucose = 10.0
        expected_output = 10.0 / 11.0  # expected output
        self.assertAlmostEqual(glycolysis_reaction2(glucose), expected_output, places=5)

    def test_glycolysis_reaction3(self):
        # Test with known inputs and outputs
        glucose = 10.0
        expected_output = 1.0 * glucose / (1.0 + glucose)  # expected output
        self.assertAlmostEqual(glycolysis_reaction3(glucose), expected_output, places=5)

    def test_citric_acid_cycle_reaction1(self):
        # Test with known inputs and outputs
        acetyl_coa = 2.0
        oxaloacetate = 2.0
        atp_citric = 15.0
        expected_output = 0.034408602150537634  # Calculated expected output
        self.assertAlmostEqual(citric_acid_cycle_reaction1(acetyl_coa, oxaloacetate, atp_citric), expected_output,
                               places=5)

    def test_citric_acid_cycle_reaction2(self):
        # Test with known inputs and outputs
        citrate = 2.0
        expected_output = 2.0 / 3.0  # Updated expected output
        self.assertAlmostEqual(citric_acid_cycle_reaction2(citrate), expected_output, places=5)


if __name__ == '__main__':
    unittest.main()
