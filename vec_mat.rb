require "mathn"
require "matrix"

class Vector
	def *(other)
		if (other.class == Vector)
			s = 0
			(0...self.size).map{|i| s = s + self[i] * other[i]}
			return s
		else
			a = self.map{|p| p * other}
			return Vector.elements(a,false)
		end
	end

	def -@
		return self * (-1)
	end

	def Vector.zero(n)
		a = (0...n).map{0}
		return Vector.elements(a,false)
	end

	def normal(b)
		w = self
		b.map{|p| w = w - p * ( ( self * p) / ( p * p ) ) }
		w
	end

	def proj(b)
		return self - self.normal(b)
	end

	def clean
		s = 1
		self.map{|p| s = s.lcm(p.denominator)}
		a = self.map{|p| p * s}

		s = 0
		a.map{|p| s = s.gcd(p.numerator)}
		if s != 0
			a = a.map{|p| p / s}
		end

		return Vector.elements(a,false)
	end

	def to_mat(n,m)
		a = (0...n).map{|i|
			(0...m).map{|j|
				self[i * m + j]
			}
		}
		return Matrix.columns(a)
	end

end

def schmidt(v)
	if v.size == 0
		return []
	end
	w = []
	i = 0
	zero = Vector.zero(v[0].size)
	for k in 0...v.size
		w[i] = v[k].normal(w).clean
		if w[i] == zero
			w[i] = nil
			w = w.compact
		else
			i = i + 1
		end
	end
	return w
end

class Matrix

	def Matrix.elementary(k,l,n)
		a = (0...n).map{(0...n).map{0}}
		a[k][l] = 1
		return Matrix.rows(a,false)
	end

	def anti_sym?
		n = self.row_size
		for i in 0...n
			for j in 0...n
				if self[i,j] != - self[j,i]
					return false
				end
			end
		end

		return true
	end


	def to_vec
		a = []
		m = self.transpose.to_a
		m.map{|p| a = a + p}
		return Vector.elements(a,false)
	end

	def bracket(other)
		return self * other - other * self
	end

	def mprint
		a = self.to_a
		a.map{|p|
			p.map{|q|
				print q
				print " "
			}
			print "\n"
		}
		nil
	end


end

