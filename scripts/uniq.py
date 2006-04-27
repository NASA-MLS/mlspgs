def uniq(s):
   #"""Return a list of the elements in s, but without duplicates.

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
