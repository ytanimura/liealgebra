load "liealg.rb"

class LinearLiealgebra < Liealgebra

	attr_accessor( "m" , "mb" )

	def initialize(m)
		n = m[0].row_size
		a = m.map{|p| p.to_vec}
		self.mb = Subspace.new(a)
		self.m = self.mb.map{|p| p.to_mat(n,n)}

		super(nil)
	end

	def size
		return self.m.size
	end

	def to_a
		return (0...self.size).map{|i| self[i]}
	end

	def [](n)
		return self.m[n]
	end

	def str_const
		a = self.m.map{|p| p.to_vec}
		zero = Vector.zero(a.size)
		c = (0...a.size).map{(0...a.size).map{zero}}
		
		for i in 0...a.size
			for j in (i+1)...a.size
				bra = self[i].bracket(self[j])
				c[i][j] = bra.to_vec.wrt(a)
				c[j][i] = - c[i][j]
			end
		end

		self.c = Matrix.rows(c,false)
	end

	def mprint
		for i in 0...(self.size - 1)
			self[i].mprint
			print "\n"
		end
		self[self.size - 1].mprint
	end

	def sublinliealg(b)
		return LinearLiealgebra.new(b.map{|p| p.in_mat(self)})
	end

	def LinearLiealgebra.gl(n)
		e = []
		for i in 0...n
			for j in 0...n
				e << Matrix.elementary(j,i,n)
			end
		end

		return LinearLiealgebra.new(e)
	end
	
	def LinearLiealgebra.sl(n)
		m = (1...n).map{|i| Matrix.elementary(0,0,n) - Matrix.elementary(i,i,n)}

		for i in 0...n
			for j in 0...n
				if i != j
					m << Matrix.elementary(j,i,n)
				end
			end
		end

		m = LinearLiealgebra.new(m)
		return m
	end

	def LinearLiealgebra.o(p,q)
		n = p + q

		m = []

		for j in 0...p
			for i in 0...j
				m << Matrix.elementary(i,j,n) - Matrix.elementary(j,i,n)
			end
		end


		for j in p...n
			for i in 0...p
				m << Matrix.elementary(i,j,n) + Matrix.elementary(j,i,n)
			end
			for i in p...j
				m << Matrix.elementary(i,j,n) - Matrix.elementary(j,i,n)
			end
		end

		return LinearLiealgebra.new(m)
	end

	def LinearLiealgebra.heisenberg(n)
		m = (1...n+1).map{|i| Matrix.elementary(0,i,n+2)}
		m = m + (1...n+1).map{|i| Matrix.elementary(i,n+1,n+2)}
		m << Matrix.elementary(0,n+1,n+2)

		return LinearLiealgebra.new(m)
	end

end

class Liealgebra

	def derivation
		c = self.str_const
		n = self.size
		a = []
		for i in 0...n
			for j in i...n
				for l in 0...n
					s = Array.new(n*n,0)
					for k in 0...n
						s[l + n * k] = s[l + n * k] + c[i,j][k]
						s[k + n * i] = s[k + n * i] - c[k,j][l]
						s[k + n * j] = s[k + n * j] - c[i,k][l]
					end
					a << s
				end
			end
		end

		a = Matrix.rows(a,true).gauss.map{|p| p.to_mat(n,n)}
		return LinearLiealgebra.new(a)
	end

end

class Vector

	def in_mat(m)
		s = Matrix.zero(m[0].row_size)
		for i in 0...m.size
			s = s + m[i] * self[i]
		end
		return s
	end

end

class Matrix

	def wrt(m)
		return self.to_vec.wrt(m.mb)
	end

end



class Subspace

	def sublinliealg(m)
		return m.sublinliealg(self)
	end

end
