part of 'core.dart';

/// # en
/// - Normal distribution, where [miu] is the mean and [sigma] is a non-negative standard deviation
/// - Reference:
/// - - https://en.wikipedia.org/wiki/Normal_distribution
/// - - https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
/// # zh
/// - 正态分布，其中miu是均值，标准差sigma是非负数
/// - 参考：
/// - - https://en.wikipedia.org/wiki/Normal_distribution
/// - - https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
class Normal extends ContinuousDistribution with ErrorFunction {
  final double sigma;
  final double miu;
  double? _nextRng;

  Normal({required this.miu, required this.sigma, math.Random? random})
      : assert(sigma >= 0),
        super(random: random);

  Normal.standard({math.Random? random})
      : sigma = 1,
        miu = 0,
        super(random: random);

  double get variance => sigma * sigma;

  @override
  double cdf(num x) {
    return 0.5 * (1 + erf((x - miu) / (sigma * math.sqrt2)));
  }

  @override
  double pdf(num x) {
    return math.exp(-(x - miu) * (x - miu) / (2 * variance)) / math.sqrt(2 * math.pi * variance);
  }

  @override
  double ppf(double p) {
    assert(p >= 0 && p <= 1);
    return miu + sigma * math.sqrt2 * erfInv(2 * p - 1);
  }

  @override
  double rng() {
    if (_nextRng != null){
      double value = _nextRng!;
      _nextRng = null;
      return value;
    }
    double u1 = _random.nextDouble();
    double u2 = _random.nextDouble();
    double r = math.sqrt(-2.0 * math.log(u1));
    double theta = 2.0 * math.pi * u2;
    _nextRng = r * math.sin(theta) * sigma + miu;
    return r * math.cos(theta) * sigma + miu;
  }
}
