import quadratic_sieve

nums = [16921456439215439701,
        46839566299936919234246726809,
        6172835808641975203638304919691358469663,
        3744843080529615909019181510330554205500926021947]


if __name__ == "__main__":
    for n in nums:
        print("-------------------------------------------")
        print(f"n = {n}")
        tic = time.perf_counter()
        f = quadratic_sieve.quadratic_sieve(n)
        toc = time.perf_counter()
        print(f"Time: {toc - tic:0.4f} seconds")
        if f:
            print(f"{n} = {f[0]} * {f[1]}")