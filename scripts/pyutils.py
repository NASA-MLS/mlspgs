def diff(l,u):
   """Return the difference of the lists l, u"""
   from sets import Set


 #  print "lll"+l+"lll"
 #  print "uuu"+u+"uuu"
   
   s = Set(l) - Set(u)
   #print s
   # remove unwanted chars from the result Set([<>])
   w = str(s).split('[')[1].split(']')[0].replace("'","").replace(',','').split(' ')
   w.sort()
   w = str(w).replace("'","").replace('[','').replace(']','')
   return w

def uniq(s):
   """Return a list of the elements in s, but without duplicates."""

   alist = []
   
   # remove unwanted chars
   w = s.replace("'",'')
   w = w.replace(' ','')
   w = w.replace(',',' ')  
   w = w.replace('[','')  
   w = w.replace(']','')

   # loop over each item and add to alist if not already present
   [alist.append(item) for item in w.split() if not alist.count(item)]

   # make string of comma separated items
   return ', '.join(alist)

#def uniq(list): never tried this
#
#      set1 = set(list)
#      list = [s for s in set1]
#      return list
