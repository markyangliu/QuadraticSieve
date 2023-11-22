def gcd(a, b):
  '''
  Calculates the gcd based on the Euclidean Algorithm.
  '''
  a = abs(a)
  b = abs(b)
  if a<b:
    #swap a and b so that a is larger
    temp = a
    a = b
    b = temp

  if b == 0:
    return a

  remainder = a % b
  while remainder:
    #loop until remainder is 0

    a = b
    b = remainder
    remainder = a % b

  return b

def generate_primes(n):
    '''
    Generates primes up to a bound n using the sieve of Eratosthenes
    '''

    primes = []
    not_prime = [False] * (n + 1)

    for i in range(2, n + 1):
        if not_prime[i]:
            continue
        # if i is a prime, remove all factors of i starting at i^2
        for j in range(i * i, n + 1, i):
            not_prime[j] = True
        
        primes.append(i)

    return primes

def quadratic_res(a, p):
    '''
    Calculates the legendre symbol to determine if a is a quadratic residue mod p.
    '''
    return pow(a, (p - 1) // 2, p)

def calc_interval(n):
    L = pow(math.e, math.sqrt(math.log(n) * math.log(math.log(n))))
    return math.ceil(pow(L, 3 * math.sqrt(2) / 4))
