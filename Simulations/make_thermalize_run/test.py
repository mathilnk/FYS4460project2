class test:
	def __init__(self):
		line = "N = 6;"		
		print line
		exec "self."+line
		print self.N
t = test();
