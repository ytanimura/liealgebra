load "subsp.rb"

class Liealgebra
	attr_accessor("c" , "b" , "r" , "g1" , "s" , "cent")

	def initialize(c)
		self.c = c
		self.b = nil
		self.r = nil
		self.g1 = nil
		self.cent = nil
	end

	def Liealgebra.abelian(size)
		a = Vector.zero(size)
		b = Array.new(size,a)
		Liealgebra.new(Matrix.rows(Array.new(size,b),false))
	end

	def size
		return self.c.row_size
	end

	def str_const
		return self.c
	end

	def jacobian?
		c = self.str_const
		n = c.row_size

		s = (0...n).map{|i|
			(0...n).map{|j|
				(0...n).map{|k|
					a = (0...n).map{|l|
						s = 0
						for m in 0...n
							s = s + c[j,k][m] * c[i,m][l]
						end
						s
					}
					Vector.elements(a,false)
				}
			}
		}

		zero = Vector.zero(n)
		for i in 0...n
			for j in (i+1)...n
				for k in (j+1)...n
					if s[i][j][k] + s[j][k][i] + s[k][i][j] != zero
						print i,j,k
						return false
					end
				end
			end
		end
		return true

	end

	def liealgebra?
		c = self.str_const

		if !(c.square?)
			return false
		end

		if c[0,0].size != c.row_size
			return false
		end

		if !(c.anti_sym?)
			return false
		end

		if !(self.jacobian?)
			return false
		end

		return true
	end

	def killing_form
		if self.b == nil
			c = self.str_const
			b = (0...self.size).map{|i|
				(0...self.size).map{|j|
					s = 0
					for k in 0...self.size
						for l in 0...self.size
							s = s + c[j,k][l] * c[i,l][k]
						end
					end
					s
				}
			}
			self.b = Matrix.rows(b,false)
		end
		return self.b
	end

	def semi_simple?
		self.killing_form
		return self.b.det != 0
	end

	def bracket(v,w)
		if v.class == Vector
			c = self.str_const
			a = Vector.zero(self.size)
			for i in 0...self.size
				for j in 0...self.size
					a = a + c[i,j] * ( v[i] * w[j] )
				end
			end
			return a
		elsif v.class == Subspace
			a = []
			for i in 0...v.size
				for j in 0...w.size
					a << self.bracket(v[i],w[j])
				end
			end
			return Subspace.new(a)
		end
	end

	def subalgebra?(b)
		return self.bracket(b,b).in(b)
	end

	def ideal?(b)
		return self.bracket(Subspace.std(self.size),b).in(b)
	end

	def derived
		c = self.str_const
		if self.g1 == nil
			a = []
			for i in 0...self.size
				for j in 0...self.size
					a << c[i,j]
				end
			end
			self.g1 = Subspace.new(a)
		end
		return self.g1
	end

	def radical
		if self.r == nil
			g1 = self.derived.to_mat
			b = self.killing_form
			m = (b * g1).transpose
			self.r = m.gauss
		end
		return self.r
	end

	def nilpotent_radical
		if self.s == nil
			self.s = self.derived.cap(self.radical)
		end
		return self.s
	end

	def center
		if self.cent == nil
			n = self.size
			cent = Subspace.std(n)
			c = self.str_const
			for i in 0...n
				p = (0...n).map{|j|
					(0...n).map{|k|
						c[i,k][j]
					}
				}
				p = Matrix.rows(p,true)
				cent = cent.cap(p.gauss)
			end
			self.cent = cent
		end
		return self.cent
	end
		
	def subalgebra(u)
		c = []

		c = (0...u.size).map{[]}

		zero = Vector.zero(u.size)

		for i in 0...u.size
			c[i][i] = zero
		end

		for i in 0...u.size
			for j in (i+1)...u.size
				c[i][j] = self.bracket(u[i],u[j]).wrt(u)
				c[j][i] = - c[i][j]
			end
		end
		return Liealgebra.new(Matrix.rows(c,false))
	end

	def quotient(u)
		v = u.normal
		c = (0...v.size).map{[]}

		zero = Vector.zero(v.size)

		for i in 0...v.size
			c[i][i] = zero
		end

		for i in 0...v.size
			for j in (i+1)...v.size
				c[i][j] = self.bracket(v[i],v[j]).wrt(v)
				c[j][i] = - c[i][j]
			end
		end

		return Liealgebra.new(Matrix.rows(c,false))
	end

	def Liealgebra.heisenberg(n)
		c = (0...(2 * n + 1)).map{|i|
			(0...(2 * n + 1)).map{|j|
				s = Array.new(2 * n,0)
				if i == j + n && j < n
					s << -1
				elsif j == i + n && i < n
					s << 1
				else
					s << 0
				end
				Vector.elements(s,false)
			}
		}

		return Liealgebra.new(Matrix.rows(c,false))
	end

	def Liealgebra.ladder(n)
		zero = Vector.zero(n+1)
		c = (0...n+1).map{(0...n+1).map{zero}}
		for i in 1...n
			a = Array.new(n+1,0)
			a[i+1] = 1
			c[0][i] = Vector.elements(a,true)
			a[i+1] = -1
			c[i][0] = Vector.elements(a,false)
		end

		return Liealgebra.new(Matrix.rows(c))
	end

end

class Subspace

	def subalgebra(c)
		return c.subalgebra(self)
	end

	def subalgebra?(c)
		return c.subalgebra?(self)
	end

	def ideal?(c)
		return c.ideal?(self)
	end

end

