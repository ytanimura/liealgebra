load "vec_mat.rb"

class Subspace < Array

	def initialize(v)
		return super(schmidt(v))
	end

	def Subspace.std(n)
		Matrix.identity(n).to_subsp
	end

	def in(other)
		zero = Vector.zero(self[0].size)
		self.map{|p|
			if p.normal(other) != zero
				return false
			end
		}
		return true
	end

	def ==(other)
		if other.class != Subspace
			return false
		elsif self.size != other.size
			return false
		else
			return self.in(other)
		end
	end

	def cap(other)
		b = self.map{|p| p.proj(other)}
		return Subspace.new(b)
	end

	def +(other)
		return Subspace.new(super)
	end

	def to_mat
		Matrix.columns(self)
	end

	def normal
		b = Subspace.std(self[0].size)
		s = b.map{|p| p.normal(self)}
		return Subspace.new(s)
	end

	def bprint
		self.map{|p| p p}
		nil
	end

end

class Vector

	def wrt(b)
		v = b.map{|p| (self * p) / (p * p)}
		return Vector.elements(v,false)
	end

	def P(i,j)
		a = self.to_a
		s = a[i]
		a[i] = a[j]
		a[j] = s
		return Vector.elements(a,false)
	end

end

class Matrix

	def to_subsp
		a = self.to_a
		a = a.map{|p| Vector.elements(p,false)}
		return Subspace.new(a)
	end

	def P(i,j)
		a = self.to_a
		s = a[i]
		a[i] = a[j]
		a[j] = s
		return Matrix.rows(a,false)
	end

	def Q(i,r)
		a = self.to_a
		a[i] = a[i].map{|p| p * r}
		return Matrix.rows(a,false)
	end

	def R(i,j,r)
		a = self.to_a
		a[j] = (0...a[i].size).map{|k| a[j][k] + r * a[i][k]}
		return Matrix.rows(a,false)
	end

	def rP(i,j)
		a = self.transpose.to_a
		s = a[i]
		a[i] = a[j]
		a[j] = s
		return Matrix.columns(a)
	end

	def where_is_zero(i,k)
		for j in i...self.row_size
			if self[j,k] != 0
				return [ self[j,k] , self.P(i,j) ]
			end
		end
		return [0,self]
	end

	def sub_solution(i)
		if self.row_size < i
			a = (0...self.row_size).map{|j| - self[j,i]}
			a = a + Array.new(i - self.row_size,0)
		else
			a = (0...i).map{|j| - self[j,i]}
		end
		a << 1
		a = a + Array.new(self.column_size - i - 1,0)
		return Vector.elements(a,false)
	end

	def sub_gauss(i,s)
		m = self
		m = m.Q(i,1/s)
		for j in 0...m.row_size
			if j != i
				m = m.R( i , j , - m[j,i])
			end
		end
		return m
	end

	def gauss
		com = []
		m = self

		for i in 0...m.column_size

			for j in i...m.column_size
				a = m.where_is_zero(i,j)
				m = a[1]
				if a[0] != 0
					if j != i
						m = m.rP(i,j)
						com << [i,j]
					end
					break
				end
			end
			if a[0] == 0
				break
			end
			m = m.sub_gauss(i,a[0])
		end

		if a[0] != 0
			return Subspace.new([])
		end

		sol = (i...m.column_size).map{|j| m.sub_solution(j)}

		sol = sol.map{|p|
			for i in 0...com.size
				p = p.P(com[com.size - i - 1][0],com[com.size - i - 1][1])
			end
			p
		}

		return Subspace.new(sol)
	end

end



