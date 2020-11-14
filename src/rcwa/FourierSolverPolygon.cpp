#include "FourierSolverPolygon.h"

#include "BlockToeplitzMatrixXcs.h"
#include <iostream>
#include <limits>
#include <cmath>
#include "MathFunction.h"

using namespace std;
using namespace Eigen;

constexpr scalar tolerance = 1e-14;

FourierSolverPolygon::FourierSolverPolygon(const Eigen::VectorXs& alphas,
	scalar dx, scalar dy):alphas_(alphas), dx_(dx), dy_(dy), offset_(0.)
{
	if (alphas_.size() < 3)
	{
		cerr << "[ERROR]: the number of vertices of a polygon should at least be 3!\n";
		exit(-1);
	}

	for (int v = 0; v < alphas_.size(); ++v)
	{
		if (alphas_(v) <= 0.)
		{
			cerr << "[ERROR]: only support positive alpha value!\n";
			exit(-1);
		}

		if (alphas_(v) > 0.5 * dx || alphas_(v) > 0.5 * dy)
		{
			cerr << "[WARNING]: dangerous alpha value, it is possibly larger than the region!\n";
		}
	}


}

void FourierSolverPolygon::solve(scalex f1, // filling
	scalex f0, // background
	int nx, int ny,
	Eigen::MatrixXcs& Ff,
	const std::vector<int>* layouts,
	std::vector<Eigen::MatrixXcs>* dFfs)
{
	bool evaluate_deriv = false;
	if (layouts != nullptr)
	{
		evaluate_deriv = true;
	}
	// clog << "f1: " << f1 << ", f0: " << f0 << endl;
	
	Ff.resize(nx*ny, nx*ny);

	int Nv = alphas_.size();
	Matrix2Xs gammas(2, Nv);

	for (int v = 0; v < Nv; ++v)
	{
		scalar alpha = alphas_(v);
		gammas(0, v) = alpha * cos(offset_-v*2.*Pi/Nv);
		gammas(1, v) = alpha * sin(offset_-v*2.*Pi/Nv);

		// clog << "gamma: " << gammas.col(v).transpose() << endl;
	}

	Matrix2Xs edges(2, Nv);
	for (int v = 0; v < Nv; ++v)
	{
		edges.col(v) = (v == Nv-1) ? (gammas.col(0) - gammas.col(v)) : (gammas.col(v+1) - gammas.col(v));
	}

	scalar area = 0.;
	for (int v = 0; v < Nv; ++v)
	{
		scalar a = alphas_(v);
		scalar b = (v == Nv - 1) ? alphas_(0) : alphas_(v+1);
		// clog << "a: " << a << ", b: " << b << endl;
		area += 0.5 * a * b * sin(2*Pi/Nv);
	}
	// clog << "area: " << area << endl;

	VectorXs darea;
	if (evaluate_deriv)
	{
		int Nd = layouts->size();
		darea.resize(Nd);
		darea.setZero();
		for (int i = 0; i < Nd; ++i)
		{
			int v = (*layouts)[i];
			if (v < 0 || v >= Nv)
			{
				cerr << "[ERROR]: invalid layout!\n";
				return;
			}
			scalar b = (v == Nv - 1) ? alphas_(0) : alphas_(v+1);
			darea(i) += 0.5 * b * sin(2*Pi/Nv);
			scalar a = (v == 0) ? alphas_(Nv-1) : alphas_(v-1);
			darea(i) += 0.5 * a * sin(2*Pi/Nv);
		}
	}

	int M = 2*nx-1;
	int N = 2*ny-1;

	MatrixXcs fourier2D(M, N);
	fourier2D.setZero();

	std::vector<MatrixXcs> dfourier2D;
	if (evaluate_deriv)
	{
		dfourier2D.resize(layouts->size());
		for (auto& df : dfourier2D)
		{
			df.resize(M, N);
			df.setZero();
		}
	}

	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			int m = i - (M - 1) / 2;
			int n = j - (N - 1) / 2;

			if (m == 0 && n == 0)
			{
				fourier2D(i, j) = area;
				if (evaluate_deriv)
				{
					int Nd = layouts->size();
					for (int k = 0; k < Nd; ++k)
					{
						dfourier2D[k](i, j) = darea[k];
					}
				}
			}
			else
			{
				Vector2s w;
				w << -2.*m*Pi/dx_, -2.*n*Pi/dy_;

				for (int v = 0; v < Nv; ++v)
				{
					const Vector2s& gamma = gammas.col(v);
					const Vector2s& a1 = edges.col(v);
					const Vector2s& a0 = v ? edges.col(v-1) : edges.col(Nv-1);
					
					scalex frac;
					if (fabs(w.dot(a1)) <= tolerance)
					{
						if (fabs(w.dot(a0)) <= tolerance)
						{
							cerr << "w perp to a0\n";
							cerr << "[ERROR]: we cannot handle more than two perpendicular cases!\n";
							exit(-1);
						}

						//clog << "w perpto a1\n";

						const Vector2s& gamman = gamma.normalized();
						const Vector2s& ap1 = (v == Nv-1) ? edges.col(0) : edges.col(v+1);

						if (fabs(w.dot(ap1)) <= tolerance)
						{
							cerr << "w perp to ap1\n";
							cerr << "[ERROR]: we cannot handle more than two perpendicular cases!\n";
							exit(-1);
						}

						const Vector2s& asum = a1 + a0;

						frac = 0.;
						frac += cross2(asum, gamman) * w.dot(ap1);
						frac += scalex(0., w.dot(gamman)) * cross2(a1, a0) * w.dot(ap1);
						frac += -cross2(ap1, gamman) * w.dot(a0);
						frac += cross2(ap1, a1) * w.dot(gamman);
						frac /= (-w.dot(gamman) * w.dot(a0) * w.dot(ap1));
						// cout << "frac prep: " << frac << endl;

					}
					else if (fabs(w.dot(a0)) <= tolerance)
					{
						//clog << "w perpto a0\n";
						//clog << "This case is handled in w perpto a1\n";
						frac = 0;
					}
					else
					{
						//clog << "otherwise\n";
						frac = cross2(a1, a0) / w.dot(a1) / w.dot(a0);	
					}

					//clog << "frac: " << frac << endl;
					//clog << "+=: " << frac * exp(scalex(0., w.dot(gamma))) << endl;
					fourier2D(i, j) += frac * exp(scalex(0., w.dot(gamma)));
				}

				//cout << "fourier2D(i, j): " << fourier2D(i, j) << endl;


				if (evaluate_deriv)
				{
					int Nd = layouts->size();
					for (int k = 0; k < Nd; ++k)
					{
						int v = (*layouts)[k];
						const Vector2s& gamma0 = gammas.col(v);

						// vi = v-2
						int vr = v-2;
						if (vr < 0) vr += Nv;
						const Vector2s& ar = edges.col(vr);
						if (fabs(w.dot(ar)) <= tolerance)
						{
							const Vector2s& gamma = gammas.col(vr);
							const Vector2s& gamman = gamma.normalized();

							scalex exponential_term = exp(scalex(0., w.dot(gamma)));

							const Vector2s& a0 = (vr < 0) ? edges.col(vr-1+Nv) : edges.col(vr-1);
							const Vector2s& a1 = edges.col(vr);
							const Vector2s& a2 = (v == 0) ? edges.col(Nv-1) : edges.col(v-1);

							Vector2s da2 = gamma0.normalized();

							scalar _u1 = cross2(a2, gamman);
							scalar _du1 = cross2(da2, gamman);

							scalar _v1 = w.dot(gamman) * w.dot(a2);
							scalar _dv1 = w.dot(gamman) * w.dot(da2);

							scalar _u2 = -cross2(a2, a1);
							scalar _du2 = -cross2(da2, a1);

							scalar _v2 = w.dot(a0) * w.dot(a2);
							scalar _dv2 = w.dot(a0) * w.dot(da2);

							scalar dfrac = (_du1 * _v1 - _u1 * _dv1) / _v1 / _v1
										 + (_du2 * _v2 - _u2 * _dv2) / _v2 / _v2;

							dfourier2D[k](i, j) += dfrac * exponential_term;
						}

						for (int vi = v-1; vi <= v+1; ++vi)
						{
							int v1 = vi;
							if (v1 > Nv - 1) v1 -= Nv;
							if (v1 < 0) v1 += Nv;

							const Vector2s& a1 = edges.col(v1);
							const Vector2s& gamma = gammas.col(v1);

							int v0 = vi - 1;
							if (v0 < 0) v0 += Nv;

							const Vector2s& a0 = edges.col(v0);

							Vector2s da1, da0;
							if (vi == v)
							{
								da1 = -gamma0.normalized();
								da0 = gamma0.normalized();
							}
							else if (vi == v-1)
							{
								da1 = gamma0.normalized();
								da0.setZero();
							}
							else
							{
								da1.setZero();
								da0 = -gamma0.normalized();
							}

							scalex exponential_term = exp(scalex(0., w.dot(gamma)));
							scalex dexponential_term = 0.;
							if (vi == v)
							{
								dexponential_term = exp(scalex(0., w.dot(gamma))) * scalex(0., w.dot(gamma.normalized()));
							}

							if (fabs(w.dot(a1)) <= tolerance)
							{
								// TODO
								int v2 = vi + 1;
								if (v2 > Nv - 1) v2 -= Nv;
								const Vector2s& a2 = edges.col(v2);

								Vector2s da2;
								if (vi == v - 1)
								{
									da2 = -gamma0.normalized();
								}
								else
								{
									da2.setZero();
								}

								if (fabs(w.dot(a0)) <= tolerance)
								{
									cerr << "w perp to a0\n";
									cerr << "[ERROR]: we cannot handle more than two perpendicular cases!\n";
									exit(-1);
								}

								if (fabs(w.dot(a2)) <= tolerance)
								{
									cerr << "w perp to a2\n";
									cerr << "[ERROR]: we cannot handle more than two perpendicular cases!\n";
									exit(-1);			

								}

								const Vector2s& gamman = gamma.normalized();
								scalex _u1 = -cross2(a1 + a0, gamman);
								scalex _du1 = -cross2(da1 + da0, gamman);

								scalex _v1 = w.dot(gamman) * w.dot(a0);
								scalex _dv1 = w.dot(gamman) * w.dot(da0);


								scalex _u2 = cross2(a2, gamman);
								scalex _du2 = cross2(da2, gamman);

								scalex _v2 = w.dot(gamman) * w.dot(a2);
								scalex _dv2 = w.dot(gamman) * w.dot(da2);

								scalex _u3 = scalex(0., cross2(a0, a1));
								scalex _du3 = scalex(0., cross2(da0, a1) + cross2(a0, da1));

								scalex _v3 = w.dot(a0);
								scalex _dv3 = w.dot(da0);

								scalex _u4 = -cross2(a2, a1);
								scalex _du4 = -cross2(a2, da1) - cross2(da2, a1);

								scalex _v4 = w.dot(a0) * w.dot(a2);
								scalex _dv4 = w.dot(da0) * w.dot(a2) + w.dot(a0) * w.dot(da2);

								auto GEN_DFRAC = [](scalex _u, scalex _v, scalex _du, scalex _dv)
								{
									return (_du * _v - _u * _dv) / (_v * _v);
								};

								scalex frac = _u1 / _v1 + _u2 / _v2 + _u3 / _v3 + _u4 / _v4;
								scalex dfrac = GEN_DFRAC(_u1, _v1, _du1, _dv1)
											 + GEN_DFRAC(_u2, _v2, _du2, _dv2)
											 + GEN_DFRAC(_u3, _v3, _du3, _dv3)
											 + GEN_DFRAC(_u4, _v4, _du4, _dv4);

								dfourier2D[k](i, j) += frac * dexponential_term + dfrac * exponential_term;
							}
							else if (fabs(w.dot(a0)) <= tolerance)
							{
								// do nothing
							}
							else
							{
								scalar numerator = cross2(a1, a0);
								scalar denominator = w.dot(a1) * w.dot(a0);

								scalar du = cross2(da1, a0) + cross2(a1, da0);
								scalar dv = w.dot(da1) * w.dot(a0) + w.dot(a1) * w.dot(da0);

								scalar frac = numerator / denominator;							
								scalar dfrac = (du * denominator - numerator * dv) / (denominator * denominator);

								dfourier2D[k](i, j) += frac * dexponential_term + dfrac * exponential_term;
								// cout << "+=: " << frac * dexponential_term + dfrac * exponential_term << endl;
							}

						}
					}
				}				
			}
		}
	}

	fourier2D *= (f1 - f0) / (dx_ * dy_);

	//cout << "multiply: " << (f1 - f0) / (dx_ * dy_) << endl;
	// the fourier transform from background is delta function
	fourier2D((M-1)/2, (N-1)/2) += f0;
	// cout << "======================\n";
	// cout << fourier2D << endl;
	// cout << "======================\n";
	Ff = BlockToeplitzMatrixXcs(fourier2D).toDense();

	if (layouts != nullptr)
	{
		int Nd = layouts->size();
		dFfs->resize(Nd);
		for (int k = 0; k < Nd; ++k)
		{
			dfourier2D[k] *= (f1 - f0) / (dx_ * dy_);
			(*dFfs)[k] = BlockToeplitzMatrixXcs(dfourier2D[k]).toDense();
		}
	}
}