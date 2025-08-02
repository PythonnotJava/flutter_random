import 'dart:math' as math;

part 'abs.dart';
part 'uniform.dart';

mixin ErrorFunction {
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
}