# -*- coding: utf-8 -*-
from itertools import combinations_with_replacement, product
from collections import defaultdict
from factors import factorise, unfactorise
from time import perf_counter

twosquareppcache = dict()
twosquareswithzero_cache = dict()

def mods(a, n):
    if n <= 0:
        raise Exception( "negative modulus {0}".format(n) )
    a = a % n
    if (2 * a > n):
        a -= n
    return a

def quos(a, n):
    if n <= 0:
        raise Exception( "negative modulus {0}".format(n) )
    return (a - mods(a, n))//n

def grem(w, z):
    # remainder in Gaussian integers when dividing w by z
    (w0, w1) = w
    (z0, z1) = z
    n = z0 * z0 + z1 * z1
    if n == 0:
        raise Exception( "division by zero" )
    u0 = quos(w0 * z0 + w1 * z1, n)
    u1 = quos(w1 * z0 - w0 * z1, n)
    return(w0 - z0 * u0 + z1 * u1,
           w1 - z0 * u1 - z1 * u0)

def ggcd(w, z):
    while z != (0,0):
        w, z = z, grem(w, z)
    return w

def root4(p):
    # 4th root of 1 modulo p
    if p <= 1:
        raise Exception( "too small p={0}".format(p) )
    if (p % 4) != 1:
        raise Exception( "p={0} not congruent to 1".format(p))
    k = p//4
    j = 2 # j is a trial generator for GF*(p)
    while True:
        a = pow(j,k,p) 
        b=(a*a)%p
        if b == p-1:
            return a # a is a square root of -1, i.e. a 4th root of 1
        if b != 1:
            raise Exception( "p={0} not prime".format(p) )
        j += 1

def sq2(p): # p needs to be 2 or a prime 4n+1
  if p==2:
    return (1,1)
  else:
    a = root4(p)
    return ggcd((p,0),(a,1))

def sq2a(p):
  return tuple(sorted(abs(x) for x in sq2(p)))



# two squares equalling p^n, for prime p=2 or of form 4k+1
def twosquaresprimepower(p,n):
  if (p,n) in twosquareppcache:
    return twosquareppcache[(p,n)].copy()
  
  if p==2:
    if n%2:
      twosquareppcache[(p,n)]= set( ( (2**(n//2),2**(n//2)),))
      return twosquareppcache[(p,n)].copy()
    else:
      twosquareppcache[(p,n)]=  set( ( (0,2**(n//2)),))
      return twosquareppcache[(p,n)].copy()

  if n==1:
    twosquareppcache[(p,n)]= set( (sq2a(p),) )
    return twosquareppcache[(p,n)].copy()
    
  results=set()
  for c,d in twosquaresprimepower(p,n-1):
    for a,b in twosquaresprimepower(p,1):
      results.add( tuple(sorted((a*c+b*d, abs(a*d-b*c)))))
      results.add( tuple(sorted((a*d+b*c, abs(a*c-b*d)))))
      if len(results)==1+n//2:
        twosquareppcache[(p,n)]=results
        return twosquareppcache[(p,n)].copy()

  twosquareppcache[(p,n)]=results
  return twosquareppcache[(p,n)].copy()

#----------------------------------------------------------------  
def twosquareswithzero(n):
  if type(n)==list:
    fn=list(n)
  else:  
    fn = factorise(n)
  results=set()
  multiplier=1
  for p in set(fn):
    if p%4==1 or p==2:
      if len(results)==0:
        results = twosquaresprimepower(p,fn.count(p) )
      else:
        newprime = twosquaresprimepower(p,fn.count(p))
        results = set(tuple(sorted((a*c+b*d, abs(a*d-b*c))))
                      for a,b in results
                      for c,d in newprime) | \
                  set(tuple(sorted((a*d+b*c, abs(a*c-b*d))))
                      for a,b in results
                      for c,d in newprime)
    else: # p%4==3:
      if fn[p]%2==1:
        return set()
      else:
        multiplier *= p**(fn[p]//2)

  return set( (a*multiplier,b*multiplier) for (a,b) in results)

#----------------------------------------------------------------  
twosquareswithzero_r_cache = dict()

#recursive function to find all two squares summing to n
def twosquareswithzero_r(n):

  def cache(v, s):
    if len(twosquareswithzero_r_cache)<1000000:
      twosquareswithzero_r_cache[v]=s
      
  if type(n)==list:
    fn=list(n)
    val=unfactorise(n)
  else:  
    fn = factorise(n)
    val=n

  if val in twosquareswithzero_r_cache:
    return twosquareswithzero_r_cache[val]   

  multiplier=1
  for p in set(fn):
    if p%4==3 :
      if fn.count(p)%2==1:
        cache(val, set() )
        return set()
      else:
        multiplier *= p**(fn.count(p)//2)
        for i in range(fn.count(p)):
          fn.remove(p)

  # only primes 4k+1 or 2 are are left
  if len(fn)==0:
    cache( val, set( ((0,multiplier),) )) 
    return set( ((0,multiplier),) )
  if len(fn)>=1:
    p=fn[0]
    primetwosquare = twosquaresprimepower(p,fn.count(p))
    for i in range(fn.count(p)):
      fn.remove(p)
    remainder = twosquareswithzero_r(fn)
    results = set(tuple(sorted((a*c+b*d, abs(a*d-b*c))))
                    for a,b in remainder
                    for c,d in primetwosquare) | \
              set(tuple(sorted((a*d+b*c, abs(a*c-b*d))))
                    for a,b in remainder
                    for c,d in primetwosquare)

  cache(val, set( (a*multiplier,b*multiplier) for (a,b) in results ) )
  return  set( (a*multiplier,b*multiplier) for (a,b) in results )


#----------------------------------------------------------------  
def twosquares(n):
  return set( (a,b) for (a,b) in twosquareswithzero_r(n) if a>0)

    
if __name__=="__main__":
  from random import randint

  MAXSS=1000
  start=perf_counter()
  d = defaultdict(list)
  for x,y in combinations_with_replacement(range(1,MAXSS+1),2):
    if x*x+y*y <= MAXSS**2: # to avoid picking up partial solutions
      d[x*x+y*y].append((x,y))
  end=perf_counter()

  print (end-start,"sec to generate sums of squares up to", MAXSS**2)
  print ("First 100 sums of squares")
  sd = sorted(d)
  print (sd[:100])

  A000404 = [2, 5, 8, 10, 13, 17, 18, 20, 25, 26, 29, 32, 34, 37, 40, 41, 45, 50, 52,
             53, 58, 61, 65, 68, 72, 73, 74, 80, 82, 85, 89, 90, 97, 98, 100, 101,
             104, 106, 109, 113, 116, 117, 122, 125, 128, 130, 136, 137, 145, 146,
             148, 149, 153, 157, 160, 162, 164, 169, 170, 173, 178]
  assert( set(sd[:100])&set(range(0,179)) == set(A000404) )
  print ( "A000404 match successful")

  n=1000
  while n<=MAXSS**2:
    ratio = float(len([ x for x in sd if x<=n]))/float(n)
    print ("{0:.1f}% of numbers up to {1} are the sum of 2 squares".format(ratio*100.0,n))
    n*=10
  
  def test( lis ):
    print ("testing", lis)
    for n in lis:
      #print n, sorted(twosquares(n))
      if n in d:
        if ( sorted(twosquares(n)) != d[n] ):
            print("ERROR n= ",n,", twosquares=",twosquares(n),", d =",d)
            raise Exception("mismatch")
    print ("test successful")
    
  test([2**n for n in range(1,8)] )
  test ([13,5*13, 5*5*13,5*5*5*13, 5*5*13*13, 5*13*17,5*5*13*17] )
  test ([ 2**w * 5**x * 13**y * 17**z for (w,x,y,z) in product(*[range(3)]*4)])

  
  for passes in range(3):
    start=perf_counter()
    for i in range(1,1000001):
      t = twosquares(i)
      if sorted(t) != d[i] :
        raise Exception( "MISMATCH {0} \ntwosquares(i)={1} \nd[i]={2}".format(i,sorted(twosquares(i)),d[i]))
      #if sorted(twosquareswithzero(i))<>sorted(twosquareswithzero_r(i)):
      #  print "MISMATCH in twosquareswithzero_r",i, sorted(twosquareswithzero_r(i)), sorted(twosquareswithzero(i))
      if i%100000==0:
        interval = perf_counter()-start
        print (i,"numbers tested,", interval*10.0 ,u"μs/number")
        start=perf_counter()
    
