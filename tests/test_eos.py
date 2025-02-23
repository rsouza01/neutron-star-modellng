import unittest
from src.eos import ideal_gas_eos

class TestEOS(unittest.TestCase):

    def test_ideal_gas_eos(self):
        pressure = 1e5  # Pascals
        temperature = 300  # Kelvin
        expected_rho = (1 * 1.67e-27 * pressure) / (1.38e-23 * temperature)
        self.assertAlmostEqual(ideal_gas_eos(pressure, temperature), expected_rho, places=5)

    def test_ideal_gas_eos_zero_pressure(self):
        pressure = 0  # Pascals
        temperature = 300  # Kelvin
        expected_rho = 0
        self.assertEqual(ideal_gas_eos(pressure, temperature), expected_rho)

    def test_ideal_gas_eos_zero_temperature(self):
        pressure = 1e5  # Pascals
        temperature = 0  # Kelvin
        with self.assertRaises(ZeroDivisionError):
            ideal_gas_eos(pressure, temperature)

    def test_ideal_gas_eos_negative_pressure(self):
        pressure = -1e5  # Pascals
        temperature = 300  # Kelvin
        expected_rho = (1 * 1.67e-27 * pressure) / (1.38e-23 * temperature)
        self.assertAlmostEqual(ideal_gas_eos(pressure, temperature), expected_rho, places=5)

if __name__ == '__main__':
    unittest.main()
