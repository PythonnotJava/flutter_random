import 'dart:math' as math;

final class RandomGenerator {
  RandomGenerator._internal();
  static final RandomGenerator instance = RandomGenerator._internal();

  /// Standard && normal distribution.
  double StandardNormal(math.Random rd) {
    double u1 = rd.nextDouble();
    double u2 = rd.nextDouble();
    return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2);
  }

  double Normal(math.Random rd, double sigma, double mu) {
    double u1 = rd.nextDouble();
    double u2 = rd.nextDouble();
    double z = math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2);
    return z * sigma + mu;
  }

  /// Uniform distribution.
  double Uniformm(math.Random rd, {double lb = 0.0, double ub = 1.0}) {
    return lb + (ub - lb) * rd.nextDouble();
  }

  /// Binomial Distribution.
  int Binomial(math.Random rd, {required int n, required double p}) {
    int successCount = 0;
    for (int i = 0; i < n; i++) {
      if (rd.nextDouble() < p) {
        successCount++;
      }
    }
    return successCount;
  }

  /// Chi-square distribution.
  /// - Method 1:
  /// If k independent random variables Z1, Z2, ..., Zk obey the standard normal distribution (mean 0, variance 1),
  /// then the sum of the squares of these k random variables obeys the chi-square distribution with k degrees of freedom.
  /// - Method 2: Directly generate according to the gamma distribution
  double Chisquare(math.Random rd, {required int df}) {
    double sum = 0.0;
    for (int i = 0; i < df; i++) {
      double z = StandardNormal(rd);
      sum += z * z;
    }
    return sum;
  }

  /// Exponential distribution.
  /// When solving the inverse function based on the cumulative distribution function,
  /// we get a term containing `ln(1-u)`, but `u` and `1-u` are standard uniform distributions,
  /// so we choose to use `u` instead of `1-u`.
  double Exponential(math.Random rd, {required double lambda}) {
    double u = rd.nextDouble();
    return -math.log(u) / lambda;
  }

  /// F distribution.
  double F(math.Random rd, {required int d1, required int d2}) {
    return (Chisquare(rd, df: d1) / d1) / (Chisquare(rd, df: d2) / d2);
  }

  /// Gamma distribution.
  /// Method from : https://dl.acm.org/doi/10.1145/358407.358414
  double Gamma(math.Random rd, {required double k, required double theta}) {
    if (k < 1) {
      return Gamma(rd, k: k + 1, theta: 1.0) *
          math.pow(rd.nextDouble(), 1.0 / k) *
          theta;
    }
    double d = k - 1.0 / 3.0;
    double c = 1.0 / math.sqrt(9.0 * d);
    double x, v, u;
    while (true) {
      x = StandardNormal(rd);
      v = 1.0 + c * x;
      if (v <= 0.0) {
        continue;
      }
      v = v * v * v;
      u = rd.nextDouble();
      if (u <= 1.0 - 0.0331 * x * x * x * x) {
        return d * v * theta;
      }
      if (math.log(u) <= 0.5 * x * x + d * (1.0 - v + math.log(v))) {
        return d * v * theta;
      }
    }
  }

  /// Beta distribution.
  /// When any shape parameter is non-integer,
  /// we use Cheng BB/BC method, specific source: https://dl.acm.org/doi/10.1145/359460.359482,
  /// otherwise we use two [Gamma] distribution sampling methods
  double Beta(math.Random rd, {required double a, required double b}) {
    if (a == a.toInt() && b == b.toInt()) {
      return Beta_by_Gamma2(rd, a0: a, b0: b);
    } else {
      return math.min(a, b) > 1
          ? Beta_BB(rd, a0: a, b0: b)
          : Beta_BC(rd, a0: a, b0: b);
    }
  }

  double Beta_BB(math.Random rd, {required double a0, required double b0}) {
    final a = math.min(a0, b0);
    final b = math.max(a0, b0);
    final alpha = a + b;
    final beta = math.sqrt((alpha - 2) / (2 * a * b - alpha));
    final gamma = a + 1 / beta;
    while (true) {
      final u1 = rd.nextDouble();
      final u2 = rd.nextDouble();
      final v = beta * math.log(u1 / (1 - u1));
      final w = a * math.exp(v);
      final z = u1 * u1 * u2;
      final r = gamma * v - math.log(4);
      final s = a + r - w;
      if (s + 1 + math.log(5) >= 5 * z) {
        return (a == a0) ? w / (b + w) : b / (b + w);
      }
      final t = math.log(z);
      if (s >= t) {
        return (a == a0) ? w / (b + w) : b / (b + w);
      }
      if (r + alpha * math.log(alpha / (b + w)) >= t) {
        return (a == a0) ? w / (b + w) : b / (b + w);
      }
    }
  }

  double Beta_BC(math.Random rd, {required double a0, required double b0}) {
    final a = math.max(a0, b0);
    final b = math.min(a0, b0);
    final alpha = a + b;
    final beta = 1 / b;
    final delta = 1 + a - b;
    final k1 = 8 * (0.0138889 + 0.0416667 * b) / (a * beta - 0.777778);
    final k2 = 0.25 + (0.5 + 0.25 / delta) * b;

    while (true) {
      final u1 = rd.nextDouble();
      final u2 = rd.nextDouble();

      if (u1 >= 0.5) {
        final z = u1 * u1 * u2;
        if (z <= 0.25) {
          final v = beta * math.log(u1 / (1 - u1));
          final w = a * math.exp(v);
          return (a == a0) ? w / (b + w) : b / (b + w);
        }
        if (z >= k2) continue;
        final v = beta * math.log(u1 / (1 - u1));
        final w = a * math.exp(v);
        if (alpha * math.log(alpha / (b + w)) + v - math.log(4) >=
            math.log(z)) {
          return (a == a0) ? w / (b + w) : b / (b + w);
        }
      } else {
        final y = u1 * u2;
        final z = u1 * y;
        if (0.25 * u2 + z - y >= k1) continue;
        final v = beta * math.log(u1 / (1 - u1));
        final w = a * math.exp(v);
        if (alpha * math.log(alpha / (b + w)) + v - math.log(4) >=
            math.log(z)) {
          return (a == a0) ? w / (b + w) : b / (b + w);
        }
      }
    }
  }

  double Beta_by_Gamma2(math.Random rd,
      {required double a0, required double b0}) {
    double gamma1 = Gamma(rd, k: a0, theta: 1.0);
    double gamma2 = Gamma(rd, k: b0, theta: 1.0);
    return gamma1 / (gamma1 + gamma2);
  }

  /// Dirichlet distribution.
  List<double> Dirichlet(math.Random rd, {required List<double> alpha}) {
    int dim = alpha.length;
    List<double> samples = List.filled(dim, 0.0);
    double sum = 0.0;
    for (int i = 0; i < dim; i++) {
      samples[i] = Gamma(rd, k: alpha[i], theta: 1.0);
      sum += samples[i];
    }
    return samples.map((x) => x / sum).toList();
  }

  /// Geometric distribution.
  /// Using the inverse transform method.
  /// It is best not to use `ceil`, because when [p] is 1, the denominator is 0!
  int Geometric(math.Random rd, {required double p}) {
    return (math.log(1 - rd.nextDouble()) / math.log(1 - p)).floor() + 1;
  }

  /// Gumbel distribution.
  double Gumbel(math.Random rd, {required double loc, required double scale}) {
    double u = rd.nextDouble();
    return loc - scale * math.log(-math.log(u));
  }

  /// Hypergeometric distribution.
  /// Algorithm using accurate simulation sampling.
  int Hypergeometric(math.Random rd,
      {required int N, required int K, required int n}) {
    int successes = 0;
    int remainingSuccesses = K;
    int remainingTotal = N;
    for (int i = 0; i < n; i++) {
      double p = remainingSuccesses / remainingTotal;
      if (rd.nextDouble() < p) {
        successes++;
        remainingSuccesses--;
      }
      remainingTotal--;
    }
    return successes;
  }

  /// Laplace distribution.
  double Laplace(math.Random rd, {required double mu, required double b}) {
    double u = rd.nextDouble() - 0.5;
    return mu - b * u.sign * math.log(1 - 2 * u.abs());
  }

  /// Logistic distribution.
  double Logistic(math.Random rd, {required double mu, required double s}) {
    double f = rd.nextDouble();
    return mu + math.log(f / (1 - f)) * s;
  }

  /// Lognormal distribution.
  double Lognormal(math.Random rd,
      {required double mu, required double sigma}) {
    double sd = StandardNormal(rd);
    return math.exp(mu + sigma * sd);
  }

  /// Multinomial distribution.
  List<int> Multinomial(math.Random rd,
      {required int n, required List<double> p}) {
    final k = p.length;
    final counts = List<int>.filled(k, 0);
    final cumProbs = List<double>.filled(k, 0);
    cumProbs[0] = p[0];
    for (int i = 1; i < k; i++) {
      cumProbs[i] = cumProbs[i - 1] + p[i];
    }

    for (int i = 0; i < n; i++) {
      final u = rd.nextDouble();
      for (int j = 0; j < k; j++) {
        if (u < cumProbs[j]) {
          counts[j]++;
          break;
        }
      }
    }

    return counts;
  }

  /// Poisson distribution.
  /// The Knuth method performs well when [lambda] is less than 10,
  /// but when [lambda] is large, many random numbers need to be multiplied in each loop, and the efficiency becomes low.
  /// A good solution is to use the Ahrens-Dieter method.
  /// For details, please see https://dl.acm.org/doi/10.1145/355993.355997.
  /// Our algorithm still uses the Knuth method.
  int Poisson(math.Random rd, {required double lambda}) {
    final double L = math.exp(-lambda);
    int k = 0;
    double p = 1.0;
    do {
      k += 1;
      double u = rd.nextDouble();
      p *= u;
    } while (p > L);
    return k - 1;
  }

  /// Cauchy distribution.
  double Cauchy(math.Random rd, {required double x0, required double gamma}) {
    return x0 + gamma * math.tan(math.pi * (rd.nextDouble() - 0.5));
  }

  /// Pareto distribution.
  /// The parameter [ia] is the reciprocal of alpha.
  double Pareto(math.Random rd, {required double xm, required double ia}) {
    return xm / math.pow(1 - rd.nextDouble(), ia);
  }

  /// Rayleigh distribution.
  double Rayleigh(math.Random rd, {required double sigma}) {
    return sigma * math.sqrt(-2 * math.log(1 - rd.nextDouble()));
  }

  /// Triangular distribution.
  double Triangular(math.Random rd,
      {required double a, required b, required c}) {
    var fc = (c - a) / (b - a);
    var u = rd.nextDouble();
    if (u < fc) {
      return a + math.sqrt(u * (b - a) * (c - a));
    } else {
      return b - math.sqrt((1 - u) * (b - a) * (b - c));
    }
  }

  /// Gaussian Error Function.
  /// By: https://wikimedia.org/api/rest_v1/media/math/render/svg/22c92344e736efd3a7f58953eeb9e9cf2d74fa37
  /// from: https://en.wikipedia.org/wiki/Error_function
  double erf(double x) {
    if (x == 0)
      return 0;
    else if (x < 0) return -erf(-x);
    double t = 1 / (1 + 0.3275911 * x);
    return 1.0 -
        math.exp(-x * x) *
            (0.254829592 * t -
                0.284496736 * t * t +
                1.421413741 * t * t * t -
                1.453152027 * math.pow(t, 4) +
                1.061405429 * math.pow(t, 5));
  }

  /// Wald distribution.
  /// Reference: https://www.sciencedirect.com/science/article/pii/S0167947309001078
  double Wald(math.Random rd, {required double mu, required double lambda}) {
    double phiCdf(double x) => 0.5 * (1 + erf(x / math.sqrt(2)));

    double inverseGaussianCdf(double x) {
      if (x <= 0) return 0.0;
      double term1 = math.sqrt(lambda / x) * (x / mu - 1);
      double term2 = math.sqrt(lambda / x) * (x / mu + 1);
      return phiCdf(term1) + math.exp(2 * lambda / mu) * phiCdf(-term2);
    }

    double inverseGaussianInverseCdf(double p, {double tolerance = 1e-8}) {
      double left = 1e-8;
      double right = mu * 10;
      double mid;
      while ((right - left) > tolerance) {
        mid = (left + right) / 2;
        double fMid = inverseGaussianCdf(mid);
        if (fMid < p) {
          left = mid;
        } else {
          right = mid;
        }
      }
      return (left + right) / 2;
    }

    double p = rd.nextDouble();
    return inverseGaussianInverseCdf(p);
  }

  /// Weibull distribution.
  double Weibull(math.Random rd, {required double lambda, required double k}) {
    double p = rd.nextDouble();
    return lambda * math.pow(-math.log(1 - p), 1 / k);
  }

  /// Von Mises distribution.
  /// Reference: https://academic.oup.com/jrsssc/article/28/2/152/6953743
  double Vonmises(math.Random rd, {required double k, required double mu}) {
    final double tau = 1 + math.sqrt(1 + 4 * k * k);
    final double rho = (tau - math.sqrt(2 * tau)) / (2 * k);
    final double r = (1 + rho * rho) / (2 * rho);

    while (true) {
      final double u1 = rd.nextDouble();
      final double z = math.cos(math.pi * u1);
      final double f = (1 + r * z) / (r + z);
      final double c = k * (r - f);

      final double u2 = rd.nextDouble();
      if (c * (2 - c) - u2 > 0) {
        final double u3 = rd.nextDouble();
        final double theta = (u3 - 0.5).sign * math.acos(f);
        return theta + mu;
      }
      if (math.log(c / u2) + 1 - c < 0) {
        continue;
      }
    }
  }

  /// Student_t distribution.
  /// By: https://wikimedia.org/api/rest_v1/media/math/render/svg/2fafe190857ecb2fdb4f99faf0fcc2ef1a1b2cb1
  /// from: https://en.wikipedia.org/wiki/Noncentral_t-distribution
  double Student_t(math.Random rd, {required int v, required double mu}) {
    double Z = StandardNormal(rd);
    double V = Chisquare(rd, df: v);
    return (Z + mu) * math.sqrt(v / V);
  }

  /// Frechet distribution, positive alpha indicates the shape, positive s indicates the proportion, and m indicates the position of the minimum value.
  double Frechet(math.Random rd,
      {required double alpha, required double s, required double m}) {
    return m + s / math.pow(-math.log(rd.nextDouble()), 1 / alpha);
  }
}
