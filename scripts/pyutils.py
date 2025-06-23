import sys
# First, a useful housekeeping aid
def mrClean(a, s):
  """Remove unwanted chars, replace separator with space"""
  a = a.replace("'",'')
  a = a.replace(' ','')
  a = a.replace(s,' ')  
  a = a.replace('[','')  
  a = a.replace(']','')
  return a.strip()

def diff(l,u):
   """Return the difference of the lists l, u"""
   # DeprecationWarning: the sets module is deprecated
   # If using a really old version of python, i.e. pre-2.6,
   # uncomment the next 2 lines
   # from sets import Set
   #  s = Set(l) - Set(u)


   s = set(l) - set(u)
   #print s
   # remove unwanted chars from the result Set([<>])
   # Why not tell us what chars are unwanted?
   # print("split", file=sys.stderr, end="=")
   # print(str(s).split('['), file=sys.stderr)
   # print("[0]", file=sys.stderr, end="=")
   # print(str(s).split('[')[0].replace('{','').replace('}',''), file=sys.stderr)
   # print("[1]", file=sys.stderr, end="=")
   # print(str(s).split('[')[1], file=sys.stderr)
   # This next generates an "index out of range" error
   # w = str(s).split('[')[1].split(']')[0].replace("'","").replace(',','').split(' ')

   # Do we replaced it with this
   w = str(s).split('[')[0].split(']')[0].replace("'","").replace(',','').split(' ')
   # print("w(old)", file=sys.stderr, end="=")
   # print(w, file=sys.stderr)
   # or this
   w = str(s).split('[')[0].replace('{','').replace('}','').replace("'","").replace(',','').split(' ')
   # print("w(new)", file=sys.stderr, end="=")
   # print(w, file=sys.stderr)
   # ??? Doesn't work ???
   # w = str(s).split('[')[0]
   w.sort()
   # print("sorted", file=sys.stderr, end="=")
   # print(w, file=sys.stderr)
   w = str(w).replace("'","").replace('[','').replace(']','')
   # print("w", file=sys.stderr, end="=")
   # print(w, file=sys.stderr)
   return w

def uniq(a):
   """Return a list of the elements in a, but without duplicates."""

   alist = []
   
   w = mrClean(a, ',')

   # loop over each item and add to alist if not already present
   [alist.append(item) for item in w.split() if not alist.count(item)]

   # make string of comma separated items
   return ', '.join(alist)

#def uniq(list): never tried this
#
#      set1 = set(list)
#      list = [s for s in set1]
#      return list

# $Log$
