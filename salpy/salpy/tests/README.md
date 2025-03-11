# Salpy Tests

This directory contains tests for the Salpy package. The tests are written using the unittest framework.

## Running Tests

To run all tests, navigate to the root directory of the package and run:

```bash
python -m unittest discover -s salpy/tests
```

To run a specific test file:

```bash
python -m unittest salpy.tests.test_linear_alignment
```

## Test Files

- `test_linear_alignment.py`: Tests for the linear_alignment function

## Adding New Tests

When adding new tests, follow these guidelines:

1. Create a new file named `test_<module_name>.py`
2. Import the unittest module and the function(s) to test
3. Create a class that inherits from unittest.TestCase
4. Add test methods that start with "test_"
5. Use assertions to verify the expected behavior

Example:

```python
import unittest
from salpy import some_function

class TestSomeFunction(unittest.TestCase):
    def test_basic_functionality(self):
        result = some_function(input_data)
        self.assertEqual(result, expected_output)
