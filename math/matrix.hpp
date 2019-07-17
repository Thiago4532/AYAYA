#ifndef _matrix_hpp
#define _matrix_hpp

#include <stdexcept>
#include <ostream>
#include <functional>
#include <iostream>
namespace math {

template<typename t, int n, int m>
class matrix {
public:
	matrix(t i_val = t()) {
		for (int i = 0; i < n; i++)
			for(int j = 0; j < m; j++)
				v[i][j] = i_val;
	}

	matrix(std::initializer_list<std::initializer_list<t>> const& l) {
		if(l.size() != n)
			throw std::runtime_error("matrix(initializer_list): row size error");

		auto it = l.begin();
		for(int i = 0; i < n; i++){
			if(it->size() != m)
				throw std::runtime_error("matrix(initializer_list): column size error");

			auto it2 = it->begin();
			for(int j = 0; j < m; j++)
				v[i][j] = *(it2++);

			it++;
		}
	}

	t* operator[](int i) {
		if(i < 0 || i > n)
			throw std::runtime_error("t& operator[](int i): i value error");

		return v[i];
	}

	const t* operator[](int i) const {
		if(i < 0 || i > n)
			throw std::runtime_error("const t* operator[](int i) const: i value error");

		return v[i];
	}

	matrix<t, m, n> transpose() const {
		matrix<t, m, n> ans;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				ans[i][j] = v[j][i];
		return ans;
	}

	matrix<t, m, n> operator~() const {
		return transpose();
	}

	matrix<t, n, m>& apply_hadamard(matrix<t, n, m>& mat) {
		for (int i = 0; i < n; i++)
			for(int j = 0; j < m; j++)
				v[i][j] *= mat[i][j];
		return *this;
	}

	matrix<t, n, m> operator+(matrix<t, n, m> const &rhs) const {
		matrix<t, n, m> ans;
		for(int i = 0; i < n; i++)
			for(int j = 0; j < m; j++)
				ans[i][j] = v[i][j] + rhs[i][j];
	}

	matrix<t, n, m>& operator+=(matrix<t, n, m> const& rhs) {
		for(int i = 0; i < n; i++)
			for(int j = 0; j < m; j++)
				v[i][j] += rhs[i][j];
		return *this;
	}

	matrix<t, n, m> operator-(matrix<t, n, m> const &rhs) const {
		matrix<t, n, m> ans;
		for(int i = 0; i < n; i++)
			for(int j = 0; j < m; j++)
				ans[i][j] = v[i][j] - rhs[i][j];
	}

	matrix<t, n, m>& operator-=(matrix<t, n, m> const& rhs) {
		for(int i = 0; i < n; i++)
			for(int j = 0; j < m; j++)
				v[i][j] -= rhs[i][j];
		return *this;
	}

	template<int p>
	matrix<t, n, p> operator*(matrix<t, m, p> const& rhs) const {
		matrix<t, n, p> ans;
		for(int i = 0; i < n; i++)
			for(int j = 0; j < p; j++)
				for(int k = 0; k < m; k++)
					ans[i][j] += (v[i][k] * rhs[k][j]);

		return ans;
	}

	matrix<t, n, m>& operator*=(matrix<t, m, m> rhs) {
		(*this) = (*this)*rhs;
		return *this;
	}

private:
	t v[n][m];
};

// MULTIPLICATION BY A SCALAR

template<typename t, int n, int m>
matrix<t, n, m> operator*(matrix<t, n, m> const& mat, t const& x) {
	matrix<t, n, m> ans;
	for(int i = 0; i < n; i++)
		for(int j = 0; j < m; j++)
			ans[i][j] = mat[i][j]*x;
	return ans;
}

template<typename t, int n, int m>
matrix<t, n, m> operator*(t const& x, matrix<t, n, m> const& mat) {
	matrix<t, n, m> ans;
	for(int i = 0; i < n; i++)
		for(int j = 0; j < m; j++)
			ans[i][j] = x*mat[i][j];
	return ans;
}

template<typename t, int n, int m>
matrix<t, n, m>& operator*=(matrix<t, n, m>& mat, t const& x) {
	for(int i = 0; i < n; i++)
		for(int j = 0; j < m; j++)
			mat[i][j] *= x;
	return mat;
}

// APPLY VECTORIZED FUNCTIONS
template<typename t, typename lambda, int n, int m>
matrix<t, n, m> apply(matrix<t, n, m> const& mat, lambda const& func) {
	matrix<t, n, m> ans;
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < m; j++) {
			ans[i][j] = func(mat[i][j]);
		}
	}
	return ans;
}

template<typename t, int n, int m>
matrix<t, n, m>& apply_ref(matrix<t, n, m>& mat, std::function<t(t)> const& func) {
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < m; j++) {
			mat[i][j] = func(mat[i][j]);
		}
	}
	return mat;
}

// HADAMARD PRODUCT ( % )

template<typename t, int n, int m>
matrix<t, n, m> hadamard(matrix<t, n, m> const& lhs, matrix<t, n, m> const& rhs) {
	matrix<t, n, m> ans;
	for(int i = 0; i < n; i++) 
		for(int j = 0; j < m; j++)
			ans[i][j] = lhs[i][j] * rhs[i][j];
	return ans;
}

template<typename t, int n, int m>
matrix<t, n, m> operator%(matrix<t, n, m> const& lhs, matrix<t, n, m> const& rhs) {
	return hadamard(lhs, rhs);
}

template<typename t, int n, int m>
matrix<t, n, m> operator%=(matrix<t, n, m>& lhs, matrix<t, n, m> const& rhs) {
	return lhs.apply_hadamard(rhs);
}


// POWER OPERATOR

template<typename t, int n>
matrix<t, n, n> power(matrix<t, n, n> mat, int p) {
	if(p <= 0)
		throw std::runtime_error("matrix power: power error");

	bool used = false;
	matrix<t, n, n> ans;

	while(p) {
		if(p&1) {
			if(used) ans *= mat;
			else ans = mat, used = true;
		}

		mat *= mat;
		p /= 2;
	}

	return ans;
}

template<typename t, int n>
matrix<t, n, n> operator^(matrix<t, n, n> const& mat, int p) {
	return power(mat, p);
}

template<typename t, int n>
matrix<t, n, n>& operator^=(matrix<t, n, n>& mat, int p) {
	mat = mat^p;
	return mat;
}

// OUTPUT OPERATOR
template<typename t, int n, int m>
std::ostream& operator<<(std::ostream& out, matrix<t, n, m> const& mat) {
	for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++){
			out << mat[i][j];
			if(j != m-1) out << " ";
		}
		if(i != n-1) out << "\n";
	}
	return out;
}

template<typename t, int n>
using rvector = matrix<t, n, 1>;

template<int n, int m>
using imat = matrix<int, n, m>;

template<int n, int m>
using mat = matrix<double, n, m>;

template<int n>
using ivec = rvector<int, n>;

template<int n>
using vec = rvector<double, n>;

};


#endif
