load "linliealg.rb"

class Matrix
	
	def liealghom?(g1,g2)
		n = g1.size
		b = Subspace.std(n)
		sb = (0...n).map{|i| self * b[i]}
		c = g1.str_const

		for i in 0...n
			for j in (i+1)...n
				lhs = self * c[i,j]
				rhs = g2.bracket(sb[i],sb[j])
				if lhs != rhs
					return false
				end
			end
		end

		return true
	end

end

class Representation

	attr_accessor( "g" , "m" , "f" , "im" )

	def initialize(g,m)
		self.g = g
		self.m = m
		self.f = nil
		self.im = nil
	end

	def domain
		return self.g
	end

	def size
		self.m.size
	end

	def [](key)
		self.m[key]
	end

	def to_a
		return self.m.clone
	end

	def image
		if self.im == nil
			self.im = LinearLiealgebra.new(self.m)
		end
		return self.im
	end

	def to_mat
		if self.f == nil
			a = (0...self.size).map{|i| self[i].wrt(self.image) }
			self.f = Matrix.columns(a)
		end
		return self.f
	end

	def Representation.std(g)
		c = (0...g.size).map{Matrix[[0]]}
		return Representation.new(g,c)
	end

	def representation?
		self.to_mat.liealghom?(self.g,self.image)
	end

	def mprint
		for i in 0...self.size-1
			self[i].mprint
			print "\n"
		end
		self[self.size - 1].mprint
	end

	def *(h)
		g = self.domain
		n = g.size + h.size
		c1 = g.str_const
		c2 = h.str_const

		c = (0...n).map{(0...n).map{Array.new(n,0)}}

		for i in 0...g.size
			for j in 0...g.size
				for k in 0...g.size
					c[i][j][k] = c1[i,j][k]
				end
				for k in g.size...n
					c[i][j][k] = 0
				end
			end
		end

		for i in 0...g.size
			for j in g.size...n
				for k in 0...g.size
					c[i][j][k] = 0
					c[j][i][k] = 0
				end
				for k in g.size...n
					c[i][j][k] = self[i][k - g.size , j - g.size]
					c[j][i][k] = - self[i][k - g.size , j - g.size]
				end
			end
		end

		for i in g.size...n
			for j in g.size...n
				for k in 0...g.size
					c[i][j][k] = 0
				end
				for k in g.size...n
					c[i][j][k] = c2[i - g.size , j - g.size][k - g.size]
				end
			end
		end

		for i in 0...n
			for j in 0...n
				c[i][j] = Vector.elements(c[i][j],false)
			end
		end

		return Liealgebra.new(Matrix.rows(c,false))
	end

	def b1
		g = self.domain
		if g.size == 0
			return 0
		end
		m = self.to_a
		n = m[0].row_size
		m = (0...n).map{|i|
			a = m.map{|p| p.column(i)}
			Matrix.columns(a).to_vec
		}

		m = Subspace.new(m)

		m = m.map{|p| p.to_mat(g.size,n)}

		return m
	end

	def z1
		g = self.domain
		m = self.to_a

		c = g.str_const

		gs = g.size
		if gs == 0
			return 0
		end
		vs = m[0].row_size
		n = gs * vs
		b = []
		for i in 0...gs
			for j in 0...gs
				for k in 0...vs
					a = (0...n).map{0}
					for l in 0...vs
						a[gs * l + j] = a[gs * l + j] + m[i][k,l]
						a[gs * l + i] = a[gs * l + i] - m[j][k,l]
					end
					for l in 0...gs
						a[gs * k + l] = a[gs * k + l] - c[i,j][l]
					end
					b << a
				end
			end
		end

		a = Matrix.rows(b,false).gauss

		a = a.map{|p| p.to_mat(gs,vs)}

		return a

	end

	def h1
		return self.z1.size - self.b1.size
	end

end

class Liealgebra

	def adrep
		c = self.str_const
		a = (0...self.size).map{|i|
			m = (0...self.size).map{|j|
				(0...self.size).map{|k|
					c[i,k][j]
				}
			}
			Matrix.rows(m,false)
		}
		return Representation.new(self,a)
	end

	def center2
		if self.cent == nil
			self.cent = self.adrep.to_mat.gauss
		end
		return self.cent
	end

end


class LinearLiealgebra

	def to_rep
		return Representation.new(self , self.to_a)
	end

	def *(g)
		self.to_rep * g
	end

end

