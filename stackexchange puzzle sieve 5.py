#https://math.stackexchange.com/questions/4229887/characterizing-integer-solutions-to-a2-mid-biglb22-b212-bigr/4230796#4230796

from factors import factorise
from collections import defaultdict
from modsqrt import modular_sqrt
from twosquares import root4 # fourth rooot of unity modulo a prime

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

primeoffsets = defaultdict(set)

primecount=0
for p in primes(10**8):
  primecount+=1
  if p%4 != 1: continue
  w = root4(p)

  os = set()
  for r in (-w,w):
    X = ((r-1)*(p+1)//2)%p
    x = modular_sqrt(X,p)
    if x:
      os.add(x)
      os.add(p-x)

  if len(os):
    primeoffsets[p]=os
#--------------------------------------------------------------------------
print (primecount, "primes examined, ", len(primeoffsets), "primeoffsets")

candidates = defaultdict(list)

# chunk size of rang of b values to examine values of F(b)
N=10**6

def F(n):
  return 2*n**4 + 2*n**2 + 1

def search(base):
  candidates.clear()
  
  for p in primeoffsets:
    #print ("p=",p)
    b = base - base%p
    i = -(base%p)

    while i < N:
      for j in primeoffsets[p]:
        ii = i+j

        try:
          count=0
          while ( Bv[ii]%p==0 ):
            Bv[ii]//=p
            count+=1
          if count > 1:
            candidates[ii].append ( ( p,count//2) )
        except IndexError:
          pass # we are trapping index of Bv out of range errors

      i+=p
      b+=p
      
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
  
