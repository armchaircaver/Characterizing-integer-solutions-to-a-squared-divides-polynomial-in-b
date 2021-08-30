#https://math.stackexchange.com/questions/4229887/characterizing-integer-solutions-to-a2-mid-biglb22-b212-bigr/4230796#4230796

# optional module to present factors of solutions
try:
  from factors import factorise
except:
  pass

from collections import defaultdict
from modsqrt import modular_sqrt
from twosquares import root4 # fourth rooot of unity modulo a prime

#---------------------  standard egcd and modular inverse ----------------------
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m
#-------------------------------------------------------------------------------
def sqrt_mod_p2(a,p):
    # calculates x satisfying x^2 = n mod p^2
    # assumes p is an odd prime and (a,p)=1
    # following the algorithm in https://www.johndcook.com/blog/quadratic_congruences/
    
    x = modular_sqrt(a,p)
    if x==0: return 0

    return (x -(x*x-a)*modinv(2*x,p))%(p*p)
    
#-------------------------------------------------------------------------------

def primes(MAX): # simple sieve of multiples 

  # this would be quicker and smaller using a bitarray,
  # but would require installing bitarray first
  plist = [0,0,1]+[1,0]*(MAX//2+1)
  for n in range(3,int(MAX**0.5)+1,2):
    if plist[n]: plist[n*n::2*n]= [0]*len(plist[n*n::2*n])

  for x in range(MAX):
    if plist[x]:
      yield x
#--------------------------------------------------------------------------

# calculate values x mod p^2 that satisfy 2x^4+2x^2+1 = 0 mod p^2
primeoffsets = defaultdict(set)

primecount=0
for p in primes(10**8):
  primecount+=1
  if p%4 != 1: continue

  w2 = sqrt_mod_p2(-1,p)
  os = set()

  # solve 2X^2 + 1 = +/- w2 mod p^2
  for r in (-w2,w2):

    X = (r-1)*(p*p+1)//2
    x = sqrt_mod_p2(X,p)
    if x:
      os.add(x)
      os.add(p*p-x)

  if len(os):
    primeoffsets[p]=os
    if (p<100):
      print (p,os)
#--------------------------------------------------------------------------
print (primecount, "primes examined, ", len(primeoffsets), "primeoffsets")

candidates = defaultdict(list)

# chunk size of rang of b values to examine values of F(b)
N=1000000

def F(n):
  return 2*n**4 + 2*n**2 + 1

def search(base):
  candidates.clear()
  
  for p in primeoffsets:
    #print ("p=",p)
    b = base - base%(p*p)
    i = -(base%(p*p))

    while i < N:
      for j in primeoffsets[p]:
        ii = i+j

        try:
          count=0
          while ( Bv[ii]%(p*p)==0 ):
            Bv[ii]//=(p*p)
            count+=1
          if count :
            candidates[ii].append ( ( p,count ) )
        except IndexError:
          pass # we are trapping index of Bv out of range errors

      i+=p*p
      b+=p*p
      
for base in range(0,1000*N,N):
  Bv = [F(n) for n in range(base,base+N)]
  print ("testing base =",base)


  search(base)

  for b in candidates:
    a=1
    for (p, e) in candidates[b]:
      a *= p**e
    if a > b+base:
      print ("solution, a=", a,"b =", base+b )
      try:
        print (factorise(a), factorise(F(base+b)) )
      except:
        pass
  
