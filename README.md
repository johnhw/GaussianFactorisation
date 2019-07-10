# GaussianFactorisation
Simple pure-Python factorisation of Gaussian (complex) integers

```python
from gaussian_factorise import complex_gcd, primes, complex_factor

# Complex GCD
print(complex_gcd(135 - 14j, 155 + 34j))
>>> -5-12j

small_primes = primes(50_000)

# Complex factorisation
print(complex_factor(135 - 14j, small_primes))
>>> [(-3-2j), (-3-2j), (-3+10j), (-1+0j)]
```