#!/usr/bin/python

def prime_factors(n):  # based on http://stackoverflow.com/questions/15347174/python-finding-prime-factors
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1 or len(factors) == 0:
        factors.append(n)
    return factors

lat_mult = 4     # choose from factors of nproc (e.g. 1 to find smallest nlat, nproc to load-balance) and/or factors matching lon (anti-aliasing: 2*nlat ~ nlon)?
lon_maxprime = 2  # small for efficient FFTs: no higher than 3 or 5?
latlon_maxprod = 1024*1024  # maximum grid size

nlat = lat_mult
nlon = 1
nfou = 1

print("Each row indicates, for the corresponding triangular truncation, the minimum lat & lon grid sizes")
print("satisfying both the anti-aliasing constraints and any (lat) multiple or (lon) factorization constraints.")
print("The T number is the minimum num_fourier but is set by num_spherical - 1 (e.g. num_spherical = 43 for T42).")
print("Dsize is a measure of the percentage change in grid size from the previous row.")
print()
print("T\tfou\tsph\tlat\tlon\tDsize\tlon factorization")

prev_size = 0
while 1:
    nsph = nfou + 1
    while 2*nlat < 3*(nsph - 1) + 1:
        nlat += lat_mult
    while nlon < 3*nfou + 1 or prime_factors(nlon)[-1] > lon_maxprime:
        nlon += 1
    size = nlat*nlon
    if size > latlon_maxprod:
        break
    if size > prev_size:
        size_chg = str(100*(size - prev_size)/size)
        prev_size = size
    else:
        size_chg = ""
    print(("T"+str(nfou)+"\t"+str(nfou)+"\t"+str(nsph)+"\t"+str(nlat)+"\t"+str(nlon)+"\t"+size_chg+"\t"+str(prime_factors(nlon))))
    nfou += 1
