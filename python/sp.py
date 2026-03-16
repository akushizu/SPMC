import numpy as np

class Sp:

	_number = None
	__id = None

	_m = None
	_x = None
	_xcell = None
	_p = None
	_pcell = None
	_weight = None

	_time = None

	_status = None

	def __init__(self, n=-1, m=-999.999, x=-999.999, xcell=-999, p=-999.999, pcell=-999, weight=0, time=-999.999, status=False) :
		self._number = n
		self.__id = id(self)

		self._m = m
		self._x = x
		self._xcell = xcell
		self._p = p
		self._pcell = pcell
		self._weight = weight

		self._time = time

		self._status = status

	@property
	def number(self):
		return self._number
	
	@number.setter
	def number(self, n):
		self._number = n

	@property
	def m(self):
		return self._m
	
	@m.setter
	def m(self, n):
		self._m = n

	@property
	def x(self):
		return self._x

	@x.setter
	def x(self, r):
		self._x = r
	
	@property
	def xcell(self):
		return self._xcell

	@xcell.setter
	def xcell(self, xn):
		self._xcell = xn
		
	@property
	def p(self):
		return self._p

	@p.setter
	def p(self, hk):
		self._p = hk
		
	@property
	def pcell(self):
		return self._pcell
	
	@pcell.setter
	def pcell(self, pn):
		self._pcell = pn

	@property
	def weight(self):
		return self._weight
	
	@weight.setter
	def weight(self, n):
		self._weight = n
	
	@property
	def time(self):
		return self._time
	
	@time.setter
	def time(self, t):
		self._time = t
	
	@property
	def status(self):
		return self._status
	
	@status.setter
	def status(self, n):
		self._status = n

class Spn(Sp):

	_dimension = None

	def __init__(self, n=-1, m = 3):
		self.dimension = m

		self._x = -999.999 * np.ones(m)
		self._xcell = -999 * np.ones(m, dtype=int)
		self._p = -999.999 * np.ones(m)
		self._pcell = -999 * np.ones(m, dtype=int)

	@property
	def dimension(self):
		return self._dimension
	
	@dimension.setter
	def dimension(self, n):
		self._dimension = n