import random
import cmath


def norm(a):
    """Return the squared norm of a complex number
        norm(z) = real(z)^2 + imag(z)^2
    """
    return (a.real) ** 2 + (a.imag) ** 2


def complex_round(z):
    """
    Round z to the nearest Gaussian prime, by independently
    rounding the real and imaginary parts.
    """
    return round(z.real) + 1j * round(z.imag)


def complex_gcd(a, b):
    """Compute the complex valued GCD of a and b using
    Euclid's algorithm"""
    if norm(b) > norm(a):
        a, b = b, a
    # for complex numbers, we need to round
    q = complex_round(a / b)
    r = a - q * b  # remainder
    if r != 0:
        return complex_gcd(b, r)
    return b


def primes(n):
    """Return primes less than n (inefficiently)"""
    primes = [2]
    for i in range(3, n, 2):
        is_prime = True
        for prime in primes:
            if i % prime == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(i)
    return primes


def factorize(n, primes):
    """Naively factorise an integer"""
    factors = []
    for prime in primes:
        while n % prime == 0:
            n = n // prime
            factors.append(prime)
        if n == 0:
            break
    return factors


def mpow2(n, p, m):
    """Compute n^p mod m via repeated squaring"""
    x = n
    y = 1
    while p != 0:
        if p & 1:
            y = (y * x) % m
        x = (x * x) % m
        p = p >> 1
    return y


def complex_factor(a, primes, include_units=True):
    """Compute the complex factorisation of a Gaussian integer a.
    As a canonical form, factors of 2 will become 1+1j; 1-1j is equally valid
    factorisation of these terms"""

    n = norm(a)
    # factorise the norm of the number
    factors = factorize(int(n), primes)

    z_factors = []

    while len(factors) > 0:
        factor = factors.pop(0)

        if factor == 2:
            # x = 0 mod 2
            # either 1+1j or 1-1j; 1+1j chosen here
            u = 1 + 1j
        elif factor % 4 == 3:
            # x = 3 mod 4, remove two copies of this factor
            u = factor
            factors.pop(
                0
            )  # remove repeated factor (note assumes factors are in order)!
        else:
            # x = 1 mod 4
            # find k, such that k^2 = -1 mod factor = (factor-1) mod factor
            n = random.randint(2, factor - 1)
            while mpow2(n, (factor - 1) // 2, factor) != factor - 1:
                n = random.randint(2, factor - 1)
            k = mpow2(n, (factor - 1) // 4, factor)

            # try dividing in k+1j
            trial_factor = complex_gcd(factor, k + 1j)
            q = complex_round(a / trial_factor)

            # if exact, we have a factor
            if norm(a - q * trial_factor) < 1e-6:  # epsilon tolerance
                u = trial_factor
            else:
                # otherwise it is the conjugate
                u = trial_factor.imag + 1j * trial_factor.real  # conjugate

        # track the remaining number so we have the final unit factor
        a = a / u
        z_factors.append(complex_round(u))

    if include_units:
        # we might have a factor of -1, 1j, or -1j -- append this two
        return z_factors + [complex_round(a)]
    else:
        return z_factors


def prod(x):
    a = 1
    for elt in x:
        a = a * elt
    return a


def gaussian_primes(n, primes, include_units=True):
    """Return Gaussian primes in the square domain from -n to n."""
    if include_units:
        # include the "units" as primes
        zprimes = [1.0, -1.0, 1.0j, -1.0j]
    else:
        zprimes = []

    for i in range(-n, n):
        for j in range(-n, n):
            if i != 0 or j != 0:
                z = i + j * 1j
                if i ** 2 + j ** 2 in primes:
                    zprimes.append(z)
                else:
                    if i == 0 and (abs(j) in primes and abs(j) % 4 == 3):
                        zprimes.append(z)
                    elif j == 0 and (abs(i) in primes and abs(i) % 4 == 3):
                        zprimes.append(z)
    return zprimes


if __name__ == "__main__":
    # -

    z = 2 + 4j
    print(f"Norm of {z} is {norm(z)}")

    z1 = 135 - 14j
    z2 = 155 + 34j
    # should be -5+12j
    print(f"Complex GCD of {z1} and {z2} is {complex_gcd(z1, z2)}")

    # primes
    small_primes = primes(50000)
    print("First ten primes ", small_primes[:10])

    # check integer factorisation
    i = 212
    print(f"Factorisation of {i} is {factorize(i, small_primes)}")

    # show modulo powers
    m = 15
    n = 31
    p = 91
    print(f"{m}^{n} (mod {p}) is  {mpow2(m, n, p)} == {(m**n)%p}")

    # factorise and show the original is the product of the factors
    z = 135 - 14j
    factors = complex_factor(z, small_primes)
    print(f"Factors of {z} = {'*'.join([str(s) for s in factors])} = {prod(factors)}")

    # a few Gaussian primes
    n = 2
    print(
        "First few Gaussian primes",
        " ".join(
            [str(s) for s in gaussian_primes(n, small_primes, include_units=False)]
        ),
    )

    print(complex_factor(5, small_primes))
