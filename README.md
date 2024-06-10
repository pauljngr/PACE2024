## Our submission to the PACE 2024 Challenge to the Exact Track and the Parameterized Track
### Optil username: *mppeg*

Michael Jünger, Paul Jünger, Petra Mutzel, Gerhard Reinelt

June 9th, 2024

### Exact solver based on branch-and-cut

**Requirements:**
- Cbc Version 2.10.7
- CMake (we used version 3.22.1)

### To build the project, just run
```bash
cmake .
make
```

### Example usage
```bash
cat tiny_test_set/website_20.gr | ./oscm
# The solution will be printed to stdout
```