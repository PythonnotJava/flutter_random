part of 'core.dart';

typedef SingleVarFunc = double Function(double);

/// # en:
/// - Adaptive Simpson algorithm is used to solve the integral.
/// Adaptive Simpson algorithm is used to solve the integral.
/// The parameters are as follows:
///  - [func]: integrand.
///  - [b], [a]: current upper and lower limits of the integral.
///  - [eps]: allowed error (the smaller the error, the more accurate it is).
///  - [depth]: recursive depth control.
///  The smaller [eps] and the larger the [depth], the more accurate the solution will be,
///  but at the same time it will be more computationally expensive.
///  Reference: https://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
/// # zh
/// - 自适应辛普森积分法
/// 采用自适应辛普森算法求解积分。
/// 参数如下：
/// - [func]：被积函数。
/// - [b], [a]：当前积分的上限和下限。
/// - [eps]：允许误差（误差越小，精度越高）。
/// - [depth]：递归深度控制。
/// [eps] 越小，[depth] 越大，解的精度越高，
/// 但同时也会增加计算开销。
/// 参考：https://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
final class AdaptiveSimpson {
  final SingleVarFunc func;
  final double a;
  final double b;
  const AdaptiveSimpson({required this.func, required this.a, required this.b});

  static double _simpsonWithCache(
          double fa, double fm, double fb, double a, double b) =>
      (fa + 4 * fm + fb) * (b - a) / 6;

  static double _adaptiveSimpsonCache(
    SingleVarFunc f,
    double a,
    double b,
    double eps,
    double whole,
    int depth,
    double fa,
    double fb,
    double fm,
  ) {
    double c = (a + b) / 2;
    double leftMid = (a + c) / 2;
    double rightMid = (c + b) / 2;

    double flm = f(leftMid);
    double frm = f(rightMid);

    double left = _simpsonWithCache(fa, flm, fm, a, c);
    double right = _simpsonWithCache(fm, frm, fb, c, b);

    double delta = left + right - whole;

    if (depth <= 0 || delta.abs() <= 15 * eps) {
      return left + right + delta / 15;
    }

    return _adaptiveSimpsonCache(
            f, a, c, eps / 2, left, depth - 1, fa, fm, flm) +
        _adaptiveSimpsonCache(f, c, b, eps / 2, right, depth - 1, fm, fb, frm);
  }

  double adaptiveSimpson({double eps = 1e-8, int maxDepth = 20}) {
    double fa = func(a);
    double fb = func(b);
    double m = (a + b) / 2;
    double fm = func(m);
    double whole = _simpsonWithCache(fa, fm, fb, a, b);
    return _adaptiveSimpsonCache(func, a, b, eps, whole, maxDepth, fa, fb, fm);
  }
  
  double call() => adaptiveSimpson();
}

final class GammaFunction {
  /// # en
  /// - Approximate implementation of the gamma function based on the Lanczos approximation algorithm
  /// - Reference: https://en.wikipedia.org/wiki/Lanczos_approximation
  /// # zh
  /// - 基于Lanczos approximation算法的近似实现gamma函数
  /// - 参考：https://en.wikipedia.org/wiki/Lanczos_approximation
  static double lanczos(double z) {
    const double g = 7.0;
    const List<double> p = [
      0.99999999999980993,
      676.5203681218851,
      -1259.1392167224028,
      771.32342877765313,
      -176.61502916214059,
      12.507343278686905,
      -0.13857109526572012,
      9.9843695780195716e-6,
      1.5056327351493116e-7
    ];
    if (z < 0.5) {
      return math.pi / (math.sin(math.pi * z) * lanczos(1 - z));
    } else {
      z -= 1.0;
      double x = p[0];
      for (int i = 1; i < p.length; i++) {
        x += p[i] / (z + i);
      }
      double t = z + g + 0.5;
      return math.sqrt(2 * math.pi) * math.pow(t, z + 0.5) * math.exp(-t) * x;
    }
  }

  /// # en
  /// - Based on the adaptive Simpson integration method, where a and b are the upper and lower limits of the integral, respectively.
  /// # zh
  /// - 基于自适应辛普森积分法实现，其中a和b分别是积分的上下限
  static double adaptiveSimpson({required double z, required double a, required double b}) { 
    SingleVarFunc func = (double t) => math.pow(t, z - 1) * math.exp(-t);
    return AdaptiveSimpson(func: func, a: a, b: b)();
  }

  /// # en
  /// - Incomplete gamma function, only considers the incomplete lower gamma function
  /// - Reference: https://en.wikipedia.org/wiki/Incomplete_gamma_function
  /// # zh
  /// - 不完全gamma函数，仅考虑不完全下gamma函数
  /// - 参考：https://en.wikipedia.org/wiki/Incomplete_gamma_function
  static double adaptiveSimpsonGammaLower({required double alpha, required double x}){
    assert(x >= 0 && alpha > 0);
    if (x == 0)
      return 0;
    return adaptiveSimpson(z: alpha, a: 0, b: x);
  }
}

