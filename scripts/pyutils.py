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
   w = str(s).split('[')[1].split(']')[0].replace("'","").replace(',','').split(' ')
   w.sort()
   w = str(w).replace("'","").replace('[','').replace(']','')
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
