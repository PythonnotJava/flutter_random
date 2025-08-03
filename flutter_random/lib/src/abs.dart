part of 'core.dart';

/// # en
/// - Without passing in the [Random] class, a thread-safe private Random class is generated and passed in the global seed control.
/// The default value is null, indicating that the seed changes over time.
/// # zh
/// - 在不传入Random类的情况下，生成线程安全的私有Random类并传入全局种子控制，默认是null值，表示种子随时间变化
int? _changeableSeed;

/// # en
/// | Concept | Full Name                        | Discrete (e.g., Dice) | Continuous (e.g., Normal) | Description                                            |
/// |---------|----------------------------------|------------------------|----------------------------|--------------------------------------------------------|
/// | **PMF** | Probability Mass Function        | ✅ Yes                 | ❌ No                      | Probability of exact values, e.g., P(X = 3) = 1/6      |
/// | **PDF** | Probability Density Function     | ❌ No                  | ✅ Yes                     | Describes density, not direct probability              |
/// | **CDF** | Cumulative Distribution Function | ✅ Yes                 | ✅ Yes                     | Probability that X ≤ x                                 |
/// | **PPF** | Percent Point Function (CDF⁻¹)   | ✅ Yes                 | ✅ Yes                     | Inverse of CDF: given p, return x                      |
/// | **RNG** | Random Number Generator (sample) | ✅ Yes                 | ✅ Yes                     | Generate random samples following the distribution     |
/// # zh
/// | 概念      | 中文名                  | 离散分布（如骰子） | 连续分布（如正态） | 说明                                       |
/// |-----------|-------------------------|--------------------|--------------------|--------------------------------------------|
/// | **PMF**   | 概率质量函数            | ✅ 有               | ❌ 没有             | 离散分布中每个具体取值的概率，如 P(X = 3) = 1/6 |
/// | **PDF**   | 概率密度函数            | ❌ 没有             | ✅ 有               | 连续分布中表示密度的函数，不能直接代表概率     |
/// | **CDF**   | 累积分布函数            | ✅ 有               | ✅ 有               | 小于等于某个值的概率，如 P(X ≤ x)                |
/// | **PPF**   | 分位点函数（CDF 的反函数） | ✅ 有               | ✅ 有               | 给定概率 p，返回使 CDF(x) = p 的 x             |
/// | **RNG**   | 随机采样器              | ✅ 有               | ✅ 有               | 从该分布中生成符合规律的随机样本                |
abstract interface class StatsDistribution {
  late final math.Random _random;
  StatsDistribution({math.Random? random}){
    this._random = random ?? math.Random(_changeableSeed);
  }

  double cdf(num x);
  num ppf(double p);
  num rng();
}

mixin ErrorFunction on StatsDistribution {
  /// Gaussian Error Function.
  /// By: https://wikimedia.org/api/rest_v1/media/math/render/svg/22c92344e736efd3a7f58953eeb9e9cf2d74fa37
  /// from: https://en.wikipedia.org/wiki/Error_function
  double erf(double x){
    if (x == 0) return 0;
    else if (x < 0) return -erf(-x);
    double t = 1 / (1 + 0.3275911 * x);
    return 1.0 - math.exp(-x * x) * (
        0.254829592 * t
            - 0.284496736 * t * t
            + 1.421413741 * t * t * t
            - 1.453152027 * math.pow(t, 4) + 1.061405429 * math.pow(t, 5));
  }

  /// # en
  /// By: https://wikimedia.org/api/rest_v1/media/math/render/svg/f35ab4db6e3d906cb76de918bad2cf0ab9ea1c85
  /// # zh
  /// 来源：https://wikimedia.org/api/rest_v1/media/math/render/svg/f35ab4db6e3d906cb76de918bad2cf0ab9ea1c85
  double erfInv(double z) {
    if (z <= -1 || z >= 1) {
      throw ArgumentError("Input value must be in the open interval (-1, 1).");
    }
    final double sqrtPiOver2 = math.sqrt(math.pi) / 2;
    final coeffs = <double>[
      1,
      0.261799388,
      0.143931731,
      0.0976636195,
      0.073299079,
      0.058372501,
    ];
    double sum = 0;
    for (int k = 0; k < coeffs.length; k++) {
      int exponent = 2 * k + 1;
      sum += coeffs[k] * math.pow(z, exponent);
    }

    return sqrtPiOver2 * sum;
  }
}

abstract interface class ContinuousDistribution extends StatsDistribution {
  ContinuousDistribution({math.Random? random}) : super(random: random);
  double pdf(num x);
}
abstract interface class DiscreteDistribution extends StatsDistribution  {
  DiscreteDistribution({math.Random? random}) : super(random: random);
  double pmf(num x);
}