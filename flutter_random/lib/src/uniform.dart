part of 'core.dart';

/// # en
/// - Here we only discuss the continuous Uniform Distribution.
/// - Reference: https://en.wikipedia.org/wiki/Continuous_uniform_distribution
/// # zh
/// - 这里仅考虑连续均匀分布
/// - 参考：https://en.wikipedia.org/wiki/Continuous_uniform_distribution
class Unifrom extends ContinuousDistribution {
  final double lb;
  final double ub;
  late final double gap;

  Unifrom({required this.lb, required this.ub, math.Random? random})
      : assert(lb <= ub),
        super(random: random) {
    gap = ub - lb;
  }

  /// # en
  /// - Standard uniform distribution
  /// # zh
  /// - 标准均匀分布
  Unifrom.standard({math.Random? random})
      : lb = 0,
        ub = 1,
        gap = 1,
        super(random: random);

  @override
  double cdf(num x) {
    if (x < lb)
      return 0;
    else if (x >= ub)
      return 1;
    else
      return (x - lb) / gap;
  }

  @override
  double ppf(double p) {
    assert(p >= 0 && p <= 1);
    return lb + p * (ub - lb);
  }

  @override
  double rng() {
    return lb + _random.nextDouble() * gap;
  }

  @override
  double pdf(num x) {
    return lb <= x && x <= ub ? 1.0 / gap : 0.0;
  }
}
