from gaussian_factorise import norm, complex_round, primes, factorize
from gaussian_factorise import complex_gcd, gaussian_primes, complex_factor, prod
from collections import Counter
import random

random.seed(4114)
eps = 1e-8


def random_complex(n):
    return random.randint(-n, n) + 1j * random.randint(-n, n)


def test_norm():
    for i in range(20):
        z = random_complex(1000)
        assert norm(z) == z.imag ** 2 + z.real ** 2


def test_complex_round():
    for i in range(20):
        z = random.uniform(-1000, 1000) + random.uniform(-1000, 1000) * 1j
        zi = complex_round(z)
        assert zi.real - int(zi.real) == 0
        assert zi.imag - int(zi.imag) == 0
        assert zi.real - z.real < 1.0
        assert zi.imag - z.imag < 1.0


def test_primes():
    assert len(primes(10)) == 4
    assert len(primes(100)) == 25
    assert len(primes(1000)) == 168
    assert len(primes(10000)) == 1229
    assert primes(14) == [2, 3, 5, 7, 11, 13]


small_primes = primes(20_000)


def test_factorize():
    assert factorize(2, small_primes) == [2]
    assert factorize(9, small_primes) == [3, 3]
    assert factorize(8, small_primes) == [2, 2, 2]
    assert factorize(10, small_primes) == [2, 5]
    assert factorize(13, small_primes) == [13]
    assert factorize(500, small_primes) == [2, 2, 5, 5, 5]
    assert factorize(999, small_primes) == [3, 3, 3, 37]


def near_whole(z):
    return abs(complex_round(z) - z) < eps


def test_complex_gcd():
    assert complex_gcd(135 - 14j, 155 + 34j) == -5 - 12j
    assert complex_gcd(5, 7) == 1

    for i in range(20):
        z1 = random_complex(5000)
        z2 = random_complex(5000)
        gcd = complex_gcd(z1, z2)
        assert near_whole(z1 / gcd)
        assert near_whole(z2 / gcd)


def test_factorize():

    zprimes = gaussian_primes(100, small_primes)

    for i in range(100):
        z = random_complex(99)
        factors = complex_factor(z, small_primes)
        # product of factors should be equal
        assert abs(prod(factors) - z) < eps
        # all factors should be whole Gaussian integers and Gaussian primes
        for f in factors:
            assert near_whole(f)
            assert f in zprimes


def test_gaussian_primes():
    for n in [5, 10, 100]:
        gaussian_primes(n, small_primes)
        gaussian_primes(n, small_primes, include_units=True)

    # verify all primes have a single factor
    zs = gaussian_primes(30, small_primes, include_units=False)
    for z in zs:
        assert len(complex_factor(z, small_primes, include_units=False)) == 1

    norms = sorted(list(set([int(norm(z)) for z in zs])))
    a055025 = [
        2,
        5,
        9,
        13,
        17,
        29,
        37,
        41,
        49,
        53,
        61,
        73,
        89,
        97,
        101,
        109,
        113,
        121,
        137,
        149,
        157,
        173,
        181,
        193,
        197,
        229,
        233,
        241,
        257,
        269,
        277,
        281,
        293,
        313,
        317,
        337,
        349,
        353,
        361,
        373,
        389,
        397,
        401,
        409,
        421,
        433,
        449,
        457,
        461,
        509,
        521,
        529,
        541,
        557,
        569,
    ]
    assert norms[: len(a055025)] == a055025


def test_complex_gcd_primes():
    zs = gaussian_primes(50, small_primes, include_units=False)
    # gcd of primes p1 x p2 and p1 x p3 should be p1 (modulo sign)
    for i in range(50):
        a = random.choice(zs)
        b = random.choice(zs)
        c = random.choice(zs)
        z1 = a * b
        z2 = a * c

        gcd = complex_gcd(z1, z2)
        assert abs(gcd) == abs(a)

