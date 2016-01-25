### Recycling Chemistry Evolution Code Classes File - for questions email:sara.i.walker@asu.edu



class polymer(object):
    """A Polymeric Sequence"""
    
    def __init__(self, tot, length, seq, ID):
        #print("A new polymer has been born!")
       	self.tot = tot
       	self.length = length
       	self.seq = seq
        self.ID = ID
        
    def talk(self):
        print("Hi. I'm a polymer.", "\n")
	
