import math
import time

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

def generate_factor_base(n, B):
    '''
    Generates a factor base up to bound B.
    Removes primes if n is not a square mod p.
    '''

    factor_base = []
    primes = generate_primes(B)

    for p in primes:
        if quadratic_res(n, p) == 1:
            factor_base.append(p)
    
    return factor_base


def quadratic_res(a, p):
    '''
    Calculates the legendre symbol to determine if a is a quadratic residue mod p.
    '''
    return pow(a, (p - 1) // 2, p)


def calc_interval(n):
    L = pow(math.e, math.sqrt(math.log(n) * math.log(math.log(n))))
    return math.ceil(pow(L, 3 * math.sqrt(2) / 4))


def calc_bound(n):
    L = pow(math.e, 0.5 * math.sqrt(math.log(n) * math.log(math.log(n))))
    return int(L)



def tonelli_shanks(n, p):
    '''
    Algorithm solves x^2 = N (mod p), returns all solutions for x
    '''
    Q = p - 1
    S = 0
    while Q % 2 == 0:
        Q //= 2
        S += 1
    
    # p = 3 mod 4
    if S == 1:
        residue = pow(n, (p + 1) // 4, p)
        return residue, p - residue
    
    z = 2
    while quadratic_res(z, p) != p - 1:
        z += 1
    
    M = S
    c = pow(z, Q, p)
    t = pow(n, Q, p)
    R = pow(n, (Q + 1) // 2, p)

    while t % p != 1:
        t_2 = pow(t, 2, p)
        i = 1
        while t_2 != 1:
            t_2 = pow(t_2, 2, p)
            i += 1
        
        b = pow(c, pow(2, M - i - 1))
        M = i
        c = (b * b) % p
        t = (t * c) % p
        R = (R * b) % p 
    
    return R, p - R


def verify_candidates(factor_base, smooth_cands, x_list):
    '''
    Use trial division to verify that candidates are B-smooth.
    '''
    valid_nums = []
    valid_x = []
    for num, x in zip(smooth_cands, x_list):
        div = abs(num)
        for p in factor_base:
            while(div % p == 0):
                div = div // p
        if div == 1:
            valid_nums.append(num)
            valid_x.append(x)
        #if len(valid_nums) > len(factor_base) + 10:
        #  return valid_nums
    return valid_nums, valid_x


def build_matrix(nums, base):
    matrix = [0]*(len(base) + 1)
    for j in range(len(nums)):
        num_pos = 1 << j
        div = nums[j]
        for i in range(len(base)):
            p = base[i]
            while(div % p == 0):
                div = div // p
                matrix[i] ^= num_pos
        if div == -1:
            matrix[-1] ^= num_pos
        
    return matrix
    

def matrix_convert_to_int(lst):
    m = []
    for row in lst:
        new_row = 0
        for i in range(len(row)):
            new_row += (row[i] << i)
        m.append(new_row)

    return m

def matrix_convert_to_list(m, size):
    lst = []
    for row in m:
        new_row = []
        for i in range(size):
            if row & (1 << i):
                new_row.append(1)
            else:
                new_row.append(0)
        lst.append(new_row)

    return lst

def gauss_elim(m,size):
    # m is a list, each entry is an integer whose bits correspond to factor base.
    # size is the number of primes in the factor base i.e. len(factor_base)
    # 1 in least significant bit means smallest prime divides the num

    # returns the rref matrix in the same format, and a list of pivot columns

    num_pivots = 0
    pivots = []
    free = []
    for col in range(size):
        #print("looking at column", col)
        col_pos = 1 << col
        looking = True
        for i in range(num_pivots, len(m)):  # look at each row

            row = m[i]
            if row & col_pos:  # if row i contains a 1 in this column
                if looking:
                    looking = False  # we have found a pivot in this row
                    # swap so this is the next pivot row if not already in place
                    if i != num_pivots:
                        m[i] = m[num_pivots]
                        m[num_pivots] = row
                    #print("swapped row",i,"and row", num_pivots)
                    #print(matrix_convert_to_list(m,size))
                else:
                    # add the new pivot row to every other row
                    m[i] ^= m[num_pivots]
                    
                    #print("added row", num_pivots, "to row", i)
                    #print(matrix_convert_to_list(m, size))

        if not looking:
            num_pivots += 1
            pivots.append(col)
        else:
            free.append(col)

    #the matrix is now in echelon form
    #the following makes pivot columns have zeros above pivots (rref)

    for j in range(num_pivots):  # go through each pivot column
        col = pivots[j]
        col_pos = 1 << col
        #print(col)
        for i in range(0, j): # check each row above the row with the pivot
            #print(i,col)
            if m[i] & col_pos:
                # add row j to row i to get rid of a one above a pivot
                m[i] ^= m[j]

    #print("rref is", matrix_convert_to_list(m,size))

    #now, compute nullspace
    null = [0]*len(free)
    
    for i in range(len(free)):
        
        free_var = free[i]
        #print(free_var)
        free_var_pos = 1 << free_var
        for j in range(min(len(m),free_var)):
            row = m[j]
            
            if row & free_var_pos:
                #print("row", j, "has a 1 in position", free_var, "so new null matrix is")
                
                null[i] ^= (1 <<pivots[j])
                #print(matrix_convert_to_list(null,size))
                
                #print(null[i])
                
        null[i] ^= free_var_pos
    return m, null


def generate_pos_smooth_nums(n, factor_base, interval_start, interval_end):
    '''
    generates B-smooth numbers going forwards (positive x values)
    '''
    # sieve using f(x) = (x+s)^2 - n
    search_nums = []
    # contains the registers with log f(x)
    sieve_list = []

    # calculate f(x) values in the interval and populate sieve array
    s = math.ceil(math.sqrt(n))
    for x in range(interval_start + s, interval_end + s):
        num = (x * x) % n
        search_nums.append(num)
        sieve_list.append(round(math.log2(num)))


    # sieve with odd primes
    for p in factor_base[1:]:
        # find residues x1 and x2
        residues = tonelli_shanks(n, p)

        log_p = math.log2(p)
        # subtract from numbers with factor p
        for res in residues:
            idx = (res - s - interval_start) % p
            while idx < interval_end - interval_start:
                sieve_list[idx] = sieve_list[idx] - log_p
                idx += p
    
    threshold = 20

    # store candidates and 
    bsmooth_candidates = []
    x_list = []


    # find registers that are less than the threshold
    for i, candidate in enumerate(sieve_list):
        if candidate < threshold:
            bsmooth_candidates.append(search_nums[i])
            x_list.append(i + interval_start + s)
    
    # verify that the candidates are B-smooth
    bsmooth_nums, x_pos = verify_candidates(factor_base, bsmooth_candidates, x_list)
    
    return bsmooth_nums, x_pos

def generate_neg_smooth_nums(n, factor_base, interval_start, interval_end):
    '''
    generates B-smooth numbers going backwards (negative x values)
    '''
    # sieve using f(x) = (x+s)^2 - n
    search_nums = []
    # contains the registers with log f(x)
    sieve_list = []

    # calculate f(x) values in the interval and populate sieve array
    s = math.ceil(math.sqrt(n))
    for x in range(interval_start + s, interval_end + s, -1):
        num = abs((x * x) - n)
        search_nums.append(num)
        sieve_list.append(round(math.log2(num)))
        

    # sieve with odd primes
    for p in factor_base[1:]:
        # find residues x1 and x2
        residues = tonelli_shanks(n, p)
        log_p = math.log2(p)
        # subtract from numbers with factor p

        for res in residues:
            idx = (res - s - interval_start) % p
            # account for negative x value
            if idx < interval_start - interval_end:
                sieve_list[idx] = sieve_list[idx] - log_p
                idx -= p
                idx = abs(idx)

                while idx < interval_start - interval_end:
                    sieve_list[idx] = sieve_list[idx] - log_p
                    idx += p

    threshold = 20
    bsmooth_candidates = []
    x_list = []

    # find registers that are less than the threshold
    for i, candidate in enumerate(sieve_list):
        if candidate < threshold:
            bsmooth_candidates.append(-1 * search_nums[i])
            x_list.append(interval_start + s - i)

    # verify that the candidates are B-smooth 
    bsmooth_nums, x_pos = verify_candidates(factor_base, bsmooth_candidates, x_list)
    
    return bsmooth_nums, x_pos



def quadratic_sieve(n):
    n_len = len(str(n))
    if n_len < 25:
      SUBINT_LEN = 10000
    elif n_len < 30:
      SUBINT_LEN = 50000
    else:
      SUBINT_LEN = 100000
    #print('n = ', n)
    #print('n = ', n)

    # calculate smoothness bound B and sieve interval
    B = calc_bound(n)
    sieve_interval = max(calc_interval(n), 500000)

    #print('B: ', B, ' | Sieve Interval: ', sieve_interval)
    #print(((math.sqrt(2) - 1) * math.sqrt(n)) - 1)

    # remove unwanted primes 
    factor_base = generate_factor_base(n, B)
    #print('Factor Base Length: ', len(factor_base))
    #print('Factor Base: ', factor_base)


    # sieve for B-smooth numbers in subintervals of size SUBINT_LEN to conserve memory
    smooth_nums = []
    x_list = []

    # initialize sieve intervals going in both positive and negative directions
    start = 0
    end = min(SUBINT_LEN, sieve_interval)
    neg_start = 0
    neg_end = max(-1 * SUBINT_LEN, -1 * sieve_interval)

    while start < sieve_interval:
        # get new b-smooth numbers from subinterval and add to all smooth numbers
        smooths, x_extend = generate_pos_smooth_nums(n, factor_base, start, end)
        smooth_nums.extend(smooths)
        x_list.extend(x_extend)

        if len(smooth_nums) > len(factor_base) + 100:
              break

        smooths, x_extend = generate_neg_smooth_nums(n, factor_base, neg_start, neg_end)
        smooth_nums.extend(smooths)
        x_list.extend(x_extend)

        #print('interval:', neg_end, end, 'smooths found so far: ', len(smooth_nums))

        # stop when we have enough B-smooth numbers
        if len(smooth_nums) > len(factor_base) + 100:
              break

        # update subintervals
        start = end
        end = min(end + SUBINT_LEN, sieve_interval)

        neg_start = neg_end
        neg_end = max(neg_end - SUBINT_LEN, -1 * sieve_interval)

        
        
    #print('number of bsmooth: ', len(smooth_nums))
    #print('bsmooth nums: ', smooth_nums)


    # build binary matrix for the smooth nums
    size = len(smooth_nums)
    matrix = build_matrix(smooth_nums, factor_base)

    # generate left null space for the matrix
    matrix, null = gauss_elim(matrix, size)    


    #make potential values to test
    for null_row in null:
        y_squared = 1
        test_x = 1
        for i in range(size):
            pos = (1 << i)
            if null_row & pos:
                y_squared *= smooth_nums[i]
                test_x *= x_list[i]
        #print("test y_squared:", y_squared)
        #print("test x:", test_x)
        #print(((test_x*test_x) %n) == (y_squared % n))
        #find what factor_base nums test_y is the square of
        square_root = 1
        shrinker = y_squared
        for factor in factor_base:
            div = factor*factor
            while shrinker % div == 0:
                shrinker = shrinker // div
                square_root *= factor
        test_y = square_root
        #print("root of test_num:" , square_root, square_root * square_root == y_squared)
        if (test_x % n) != (test_y % n) and (test_x %n) != ((-test_y)%n):
            factor = gcd(test_x-test_y,n)
            if factor != 1:
                return factor, n//factor
            else:
                continue

nums = [16921456439215439701,
        46839566299936919234246726809,
        6172835808641975203638304919691358469663,
        3744843080529615909019181510330554205500926021947]

for n in nums:
    print("-------------------------------------------")
    print(f"n = {n}")
    tic = time.perf_counter()
    f = quadratic_sieve(n)
    toc = time.perf_counter()
    print(f"Time: {toc - tic:0.4f} seconds")
    if f:
      print(f"{n} = {f[0]} * {f[1]}")

