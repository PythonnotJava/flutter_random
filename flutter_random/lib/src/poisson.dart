part of 'core.dart';

class Poisson extends DiscreteDistribution {
/// Gamma 采样（近似法，使用 Marsaglia-Tsang）
//   double gammaSample(Random rng, double shape) {
//     if (shape < 1) {
//       return gammaSample(rng, 1 + shape) * pow(rng.nextDouble(), 1 / shape);
//     }
//
//     final d = shape - 1 / 3.0;
//     final c = 1 / sqrt(9 * d);
//
//     while (true) {
//       double x, v;
//       do {
//         x = normalSample(rng);
//         v = 1 + c * x;
//       } while (v <= 0);
//
//       v = v * v * v;
//       final u = rng.nextDouble();
//
//       if (u < 1 - 0.0331 * (x * x * x * x)) return d * v;
//       if (log(u) < 0.5 * x * x + d * (1 - v + log(v))) return d * v;
//     }
//   }
//
//   /// 二项采样（伯努利法，小 n）
//   int binomialSample(Random rng, int n, double p) {
//     var k = 0;
//     for (var i = 0; i < n; i++) {
//       if (rng.nextDouble() < p) k++;
//     }
//     return k;
//   }
//
//   /// 泊松采样：Gamma 方法（PG）
//   int poissonPG(Random rng, double lambda) {
//     const cutoff = 24.0;
//     const d = 7 / 8.0;
//     var k = 0;
//     var w = lambda;
//
//     while (true) {
//       if (w < cutoff) {
//         // 使用乘法法
//         final limit = exp(-w);
//         var prod = 1.0;
//         var count = 0;
//         while (prod > limit) {
//           prod *= rng.nextDouble();
//           count++;
//         }
//         return k + count - 1;
//       }
//
//       final n = (d * w).floor();
//       final x = gammaSample(rng, n.toDouble());
//       if (x > w) {
//         // Case C: Binomial
//         final p = w / x;
//         final bin = binomialSample(rng, n - 1, p);
//         return k + bin;
//       } else {
//         // Case D: Poisson
//         k += n;
//         w -= x;
//       }
//     }
//   }
  @override
  double cdf(num x) {
    // TODO: implement cdf
    throw UnimplementedError();
  }

  @override
  double pmf(num x) {
    // TODO: implement pmf
    throw UnimplementedError();
  }

  @override
  num ppf(double p) {
    // TODO: implement ppf
    throw UnimplementedError();
  }

  @override
  num rng() {
    // TODO: implement rng
    throw UnimplementedError();
  }

}