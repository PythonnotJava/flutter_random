part of 'core.dart';

/// # en
/// - Gamma distribution.The shape parameter [alpha] and scale parameter [theta] are both positive
/// - Reference:
/// - - https://dl.acm.org/doi/10.1145/358407.358414
/// - - https://en.wikipedia.org/wiki/Gamma_distribution
/// # zh
/// - Gamma分布。形状参数alpha和尺度参数theta均为正数
/// - 参考：
/// - - https://dl.acm.org/doi/10.1145/358407.358414
/// - - https://en.wikipedia.org/wiki/Gamma_distribution
class Gamma extends ContinuousDistribution {
  final double alpha;
  final double theta;
  
  late final Unifrom _standardUniform;
  Gamma({required this.alpha, required this.theta, math.Random? random})
      : assert(alpha > 0 && theta > 0), super(random: random) {
    _standardUniform = Unifrom.standard(random: this._random);
  }
  
  double get lamda => 1.0 / theta;
  double get fixCoef => GammaFunction.lanczos(alpha);

  @override
  double cdf(num x) {
    return GammaFunction.adaptiveSimpsonGammaLower(alpha: alpha, x: lamda * x) / fixCoef;
  }

  @override
  double pdf(num x) {
    return math.pow(lamda, alpha) * math.pow(x, alpha - 1) * math.exp(-lamda * x) / fixCoef;
  }

  @override
  double ppf(double p) {
    // TODO: implement ppf
    throw UnimplementedError();
  }

  @override
  double rng(){
    double _rng({required double alpha, required double theta}) {
      if (alpha < 1) {
        return _rng(alpha: alpha + 1, theta: 1.0) *
            math.pow(_random.nextDouble(), 1.0 / alpha) *
            theta;
      }
      double d = alpha - 1.0 / 3.0;
      double c = 1.0 / math.sqrt(9.0 * d);
      double x, v, u;
      while (true) {
        x = _standardUniform.rng();
        v = 1.0 + c * x;
        if (v <= 0.0) {
          continue;
        }
        v = v * v * v;
        u = _random.nextDouble();
        if (u <= 1.0 - 0.0331 * math.pow(x, 4)) {
          return d * v * theta;
        }
        if (math.log(u) <= 0.5 * x * x + d * (1.0 - v + math.log(v))) {
          return d * v * theta;
        }
      }
    }
    return _rng(alpha: alpha, theta: theta);
  }
}

