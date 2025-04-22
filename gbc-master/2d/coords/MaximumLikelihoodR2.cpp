#include "../extra/math.h"
#include "MaximumLikelihoodR2.h"

using namespace gbc;

void gbc::MaximumLikelihoodR2::compute(const VertexR2 &p, std::vector<double> &b)
{
    b.clear();
    b.resize(n_v, 0.0);
    // Boundary.
    if (computeBoundaryCoordinates(p, b)) return;
    // Interior.
    Eigen::Vector2d pp;
    pp(0) = p[0];
    pp(1) = p[1];

    Eigen::VectorXd lambda = dmlCoordinates(pp, poly);
    for (int i = 0; i < n_v; ++i)
    {
        b[i] = lambda(i);
    }
    return;
}

Eigen::VectorXd MaximumLikelihoodR2::dmlCoordinates(Eigen::Vector2d v, Eigen::MatrixXd & polygon) {
  int n = polygon.cols();
  Eigen::VectorXd lambda(n);

  // projection
  Eigen::MatrixXd shifted_poly = polygon.colwise() - v;
  Eigen::VectorXd r = shifted_poly.colwise().norm();
  Eigen::MatrixXd proj_poly(2, n);
  for (size_t j = 0; j < n; ++j) {
    proj_poly.col(j) = shifted_poly.col(j) / r(j);
  }

  bool smooth = true;
  Eigen::VectorXd r2, r1;
  if (smooth) {
    shifted_poly.col(0) = proj_poly.col(n - 1) + proj_poly.col(0);
    for (size_t i = 1; i < n; ++i) {
      shifted_poly.col(i) = proj_poly.col(i - 1) + proj_poly.col(i);
    }
    r1 = shifted_poly.colwise().norm();
    for (size_t j = 0; j < n; ++j) {
      shifted_poly.col(j) /= r1(j);
    }

    for (size_t i = 0; i < n - 1; ++i) {
      proj_poly.col(i) = shifted_poly.col(i) + shifted_poly.col(i + 1);
    }
    proj_poly.col(n - 1) = shifted_poly.col(n - 1) + shifted_poly.col(0);
    r2 = proj_poly.colwise().norm();
    for (size_t j = 0; j < n; ++j) {
      proj_poly.col(j) /= r2(j);
    }
  }

  auto F = [&n, &proj_poly](Eigen::Vector2d Phi) -> double {
    double res = 0;
    for (size_t i = 0; i < n; ++i) {
      res -= log(n + proj_poly.col(i).transpose().row(0) * Phi);
    }
    return res;
  };

  auto F1 = [&n, &proj_poly](Eigen::Vector2d Phi) -> Eigen::Vector2d {
    Eigen::Vector2d res = Eigen::Vector2d::Zero(2);
    for (size_t i = 0; i < n; ++i) {
      res -= 1.0 / (n + proj_poly.col(i).transpose().row(0) * Phi) *
             proj_poly.col(i);
    }
    return res;
  };

  auto F2 = [&n, &proj_poly](Eigen::Vector2d Phi) -> Eigen::Matrix2d {
    Eigen::Matrix2d res = Eigen::Matrix2d::Zero();
    for (size_t i = 0; i < n; ++i) {
      res += proj_poly.col(i) * (proj_poly.col(i).transpose()) /
             (pow(n + proj_poly.col(i).transpose().row(0) * Phi, 2));
    }
    return res;
  };
  Eigen::Vector2d Phi;
  Eigen::VectorXd lambda1 = Eigen::VectorXd::Ones(n);
  Phi << 0.0, 0.0;
  if (!convex_opt(F, F1, F2, &Phi)) {
    std::cout << "error" << std::endl;
    return lambda1;
  }
  for (size_t i = 0; i < n; ++i) {
    lambda[i] = 1.0 / (n + proj_poly.col(i).transpose().row(0) * Phi);
  }
  if (smooth) {
    for (size_t i = 0; i < n; ++i) {
      lambda[i] = lambda[i] / r2[i];
    }
    lambda1[0] = (lambda[0] + lambda[n - 1]) / r1[0];
    for (size_t i = 1; i < n; ++i) {
      lambda1[i] = (lambda[i] + lambda[i - 1]) / r1[i];
    }

    for (size_t i = 0; i < n - 1; ++i) {
      lambda[i] = (lambda1[i] + lambda1[i + 1]);
    }
    lambda[n - 1] = (lambda1[n - 1] + lambda1[0]);
  }
  for (size_t i = 0; i < n; ++i) {
    lambda[i] = lambda[i] / r[i];
  }

  //  std::cout << (polygon * lambda / lambda.sum() - v).transpose() <<
  //  std::endl;

  return lambda / lambda.sum();
}



bool MaximumLikelihoodR2::convex_opt(std::function<double(Eigen::Vector2d)> f,
                std::function<Eigen::Vector2d(Eigen::Vector2d)> f1,
                std::function<Eigen::Matrix2d(Eigen::Vector2d)> f2,
                Eigen::Vector2d *Phi, double error, double rho, double sigma,
                int M) {
  Eigen::Matrix2Xd H = Eigen::Matrix2Xd::Zero(2, 2);
  Eigen::Vector2d g, Delta;
  int k = 0;
  std::vector<double> E;
  E.clear();
  bool flag = false;
  Eigen::Vector2d newPhi;
  while (k < 100) {
    g = f1(*Phi);
    // std::cout << "g:" << g.transpose() << "  |g|=" << g.norm() << std::endl;
    ++k;

    E.push_back(g.norm());

    if (E.back() <= error) {
      // std::cout << "k: " << k << std::endl;
      return true;
    }
    H = f2(*Phi);
    Delta = H.colPivHouseholderQr().solve(-g);
    // std::cout << "Delta=" << Delta.transpose()
    //           << "  H*d+g=" << (H * Delta + g).transpose() << std::endl
    //           << H << std::endl;
    //  Armijo linear search
    int mm = 0, mk = 0;

    while (mm < M) {
      if (f(*Phi + pow(rho, mm) * Delta) <
          f(*Phi) + sigma * pow(rho, mm) * g.dot(Delta)) {
        mk = mm;
        break;
      }
      mm++;
    }
    double alpha = pow(rho, mk);
    if (!flag && f1(*Phi + alpha * Delta).norm() > E.back()) {
      newPhi = *Phi;
      flag = true;
    }
    *Phi += alpha * Delta;

    // std::cout << "phi:" << Phi->transpose() << "f(phi)=" << f(*Phi)
    //           << std::endl;
  }
  // std::cout << "k: --" << k << std::endl;
  // for (size_t i = 0; i < E.size(); i++) {
  //   std::cout << "->" << E.at(i);
  // }
  *Phi = newPhi;
  return false;
}