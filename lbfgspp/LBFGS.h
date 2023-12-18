// Copyright (C) 2016-2022 Yixuan Qiu <yixuan.qiu@cos.name>
// Under MIT license

#ifndef LBFGSPP_LBFGS_H
#define LBFGSPP_LBFGS_H
#include "LBFGSpp/Param.h"
#include "LBFGSpp/BFGSMat.h"
#include "LBFGSpp/LineSearchNocedalWright.h"



inline int triangular_matrix_index(int n, int i, int j) {
    assert(j < n);
    assert(i <= j); 
    return i + j * (j + 1) / 2; 
}

template<typename T>
class triangular_matrix {
	std::vector<T> m_data;
	int m_dim;
public:
	int index(int i, int j) const { return triangular_matrix_index(m_dim, i, j); }
	int index_permissive(int i, int j) const { return (i < j) ? index(i, j) : index(j, i); }
	triangular_matrix() : m_dim(0) {}
	triangular_matrix(int n, const T& filler_val) : m_data(n*(n+1)/2, filler_val), m_dim(n) {} 
	int dim() const { return m_dim; }
        const T& operator()(int i) const { return m_data[i]; } 
        T& operator()(int i)       { return m_data[i]; } 
        const T& operator()(int i, int j) const { return m_data[index(i, j)]; } 
        T& operator()(int i, int j)       { return m_data[index(i, j)]; } 
};


typedef triangular_matrix<double> flmat;


namespace LBFGSpp {

///
/// L-BFGS solver for unconstrained numerical optimization
///
template <typename Scalar,
          template <class> class LineSearch = LineSearchNocedalWright>
//          template <class> class LineSearch = LineSearchBacktracking>
//          template <class> class LineSearch = LineSearchBracketing>
//          template <class> class LineSearch = LineSearchMoreThuente>
class LBFGSSolver
{
private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Map<Vector> MapVec;

    const LBFGSParam<Scalar>& m_param;  // Parameters to control the LBFGS algorithm
    BFGSMat<Scalar> m_bfgs;             // Approximation to the Hessian matrix
    Vector m_fx;                        // History of the objective function values
    Vector m_xp;                        // Old x
    Vector m_grad;                      // New gradient
    Scalar m_gnorm;                     // Norm of the gradient
    Vector m_gradp;                     // Old gradient
    Vector m_drt;                       // Moving direction

    // Reset internal variables
    // n: dimension of the vector to be optimized
    inline void reset(int n)
    {
        const int m = m_param.m;
        m_bfgs.reset(n, m);
        m_xp.resize(n);
        m_grad.resize(n);
        m_gradp.resize(n);
        m_drt.resize(n);
        if (m_param.past > 0)
            m_fx.resize(m_param.past);
    }

public:
    ///
    /// Constructor for the L-BFGS solver.
    ///
    /// \param param An object of \ref LBFGSParam to store parameters for the
    ///        algorithm
    ///
    LBFGSSolver(const LBFGSParam<Scalar>& param) :
        m_param(param)
    {
        m_param.check_param();
    }

    ///
    /// Minimizing a multivariate function using the L-BFGS algorithm.
    /// Exceptions will be thrown if error occurs.
    ///
    /// \param f  A function object such that `f(x, grad)` returns the
    ///           objective function value at `x`, and overwrites `grad` with
    ///           the gradient.
    /// \param x  In: An initial guess of the optimal point. Out: The best point
    ///           found.
    /// \param fx Out: The objective function value at `x`.
    ///
    /// \return Number of iterations used.
    ///
    template <typename Foo>
    inline int minimize(Foo& f, Vector& x, Scalar& fx, int border=1000)
    {
        using std::abs;
        // Dimension of the vector
        const int n = x.size();
        reset(n);
        // The length of lag for objective function value to test convergence
        const int fpast = m_param.past;
        fx = f(x, m_grad);
        m_gnorm = m_grad.norm();
        if (fpast > 0)
            m_fx[0] = fx;
        // Early exit if the initial x is already a minimizer
        if (m_gnorm <= m_param.epsilon || m_gnorm <= m_param.epsilon_rel * x.norm())
        {
            return 1;
        }
        // Initial direction
        m_drt.noalias() = -m_grad;
        // Initial step size
        Scalar step = Scalar(1) / m_drt.norm();
        // Number of iterations used
        int k = 1;
        for (;;)
        {
            k++;
            if(k == border) {
                return k;
            }
            // Save the curent x and gradient
            m_xp.noalias() = x;
            m_gradp.noalias() = m_grad;
            Scalar dg = m_grad.dot(m_drt);
            const Scalar step_max = m_param.max_step;

            // Line search to update x, fx and gradient
            if(dg<0)
                  LineSearch<Scalar>::LineSearch(f, m_param, m_xp, m_drt, step_max, step, fx, m_grad, dg, x);
            else
            {
                for(int i=0; i<m_drt.size(); i++)
                    m_drt[i]*=(-1);
                return k;
            }
            // New gradient norm
            m_gnorm = m_grad.norm();
            // Convergence test -- gradient
            if (m_gnorm <= m_param.epsilon || m_gnorm <= m_param.epsilon_rel * x.norm())
            {
                return k;
            }
            // Convergence test -- objective function value
            if (fpast > 0)
            {
                const Scalar fxd = m_fx[k % fpast];
                if (k >= fpast && abs(fxd - fx) <= m_param.delta * std::max(std::max(abs(fx), abs(fxd)), Scalar(1))){
                    return k;
                }

                m_fx[k % fpast] = fx;
            }
            // Maximum number of iterations
            if (m_param.max_iterations != 0 && k >= m_param.max_iterations)
            {
                return k;
            }
            if(m_bfgs.add_correction(x - m_xp, m_grad - m_gradp))
            {
                ;
            }
            else
            {
                return k;
            }
            // Recursive formula to compute d = -H * g
            if(m_bfgs.apply_Hv(m_grad, -Scalar(1), m_drt))
            {
            // Reset step = 1.0 as initial guess for the next line search
               step = Scalar(1);
               k++;
            }
            else
            {
                return k;
            }
        }
        return k;
    }

    ///
    /// Returning the gradient vector on the last iterate.
    /// Typically used to debug and test convergence.
    /// Should only be called after the `minimize()` function.
    ///
    /// \return A const reference to the gradient vector.
    ///
    const Vector& final_grad() const { return m_grad; }

    ///
    /// Returning the Euclidean norm of the final gradient.
    ///
    Scalar final_grad_norm() const { return m_gnorm; }
    
    
    inline double scalar_product(const Vector& a, const Vector& b, double n)
    {
            double tmp = 0;
//            VINA_FOR(i, n)
//                    tmp += a(i) * b(i);
            for (int i = 0; i < n; i++)
                tmp += a[i] * b[i];
            return tmp;
    }
    
    template <typename Foo>
    inline double fast_line_search(Foo& f, int n, const Vector& x, const Vector& g, const double f0,
                    const Vector& p, Vector& x_new, Vector& g_new, double& f1)
    { // returns alpha
            const double c0 = 0.0001;
            const unsigned max_trials = 10;
            const double multiplier = 0.5;
            double alpha = 1;
            const double pg = scalar_product(p, g, n);
            
            for (int trial = 0; trial < max_trials; trial++){
                x_new.noalias() = x + p * alpha;
                f1 = f(x_new, g_new);
                if (f1 - f0 < c0 * alpha * pg) // FIXME check - div by norm(p) ? no?
                    break;
                alpha *= multiplier; 
            }
            return alpha;
    }
    
    inline void set_diagonal(flmat& m, double x) {
        for(int i = 0; i < m.dim(); i++)
            m(i, i) = x;
    }
    
    void minus_mat_vec_product(const flmat& m, const Vector& in, Vector& out)
    {
        int n = m.dim();
        for(int i = 0; i < n; i++)
        {
            double sum = 0;
            for(int j = 0; j < n; j++)
                sum += m(m.index_permissive(i, j)) * in[j];
            out[i] = -sum;
        }
    }
    
    void subtract_change(Vector& b, const Vector& a, int n)
    { // b -= a
        for(int i = 0; i < n; i++)
            b[i] -= a[i];
    }
    
    inline bool bfgs_update(flmat& h, const Vector& p, const Vector& y, const double alpha)
    {
        const double yp = scalar_product(y, p, h.dim());
        if (alpha * yp < epsilon_fl)
            return false; 
        Vector minus_hy(y);
        minus_mat_vec_product(h, y, minus_hy);
        const double yhy = -scalar_product(y, minus_hy, h.dim());
        const double r = 1 / (alpha * yp); // 1 / (s^T * y) , where s = alpha * p // FIXME   ... < epsilon
        const int n = p.size();
        for(int i = 0; i < n; i++)
                for(int j = i; j < n; j++) // includes i
                        h(i, j) += alpha * r * (minus_hy(i) * p[j] + minus_hy(j) * p[i]) +
                            + alpha * alpha * (r * r * yhy + r) * p[i] * p[j]; // s * s == alpha * alpha * p * p
        return true;
    }
    
    template <typename Foo>
    inline double bfgs(Foo& f, Vector& x, int max_iters)
    { // x is I/O, final value is returned
	int n = x.size();
	flmat h(n, 0);
	set_diagonal(h, 1);
        Vector g(n);
        double f0 = f(x, g);
	Vector g_new(g);
	Vector x_new(x);
	double f_orig = f0;
	Vector g_orig(g);
	Vector x_orig(x);
	Vector p(g);
	for(int step = 0; step < max_iters; step++)
	{
		minus_mat_vec_product(h, g, p);
		double f1 = 0;
		double alpha;
                alpha = fast_line_search(f, n, x, g, f0, p, x_new, g_new, f1);
		if(alpha == 0) {
                    break; //line direction was wrong, give up
		}
		Vector y(g_new);
		subtract_change(y, g, n);
		double prevf0 = f0;
		f0 = f1;
		x.noalias() = x_new;
		g.noalias() = g_new; // dkoes - check the convergence of the new gradient
		double gradnormsq = scalar_product(g, g, n);
		if (!(gradnormsq >= 1e-4)) //slightly arbitrary cutoff - works with fp
		{
                    break; // breaks for nans too
		}
		if (step == 0)
		{
                    const double yy = scalar_product(y, y, n);
                    if (std::abs(yy) > epsilon_fl)
                            set_diagonal(h, alpha * scalar_product(y, p, n) / yy);
		}
		bool h_updated = bfgs_update(h, p, y, alpha);
	}
	if (!(f0 <= f_orig))
	{ // succeeds for nans too
            f0 = f_orig;
            x.noalias() = x_orig;
            g.noalias() = g_orig;
	}
	return f0;
    }
};

}  // namespace LBFGSpp

#endif  // LBFGSPP_LBFGS_H
