part of 'core.dart';

/// | 概念      | 全称                             | 离散分布（如骰子） | 连续分布（如正态） | 说明                              |
/// | ------- | ------------------------------ | --------- | --------- | ------------------------------- |
/// | **PMF** | Probability Mass Function      | ✅ 有       | ❌ 没有      | 离散型的“每个点的概率”，如 $P(X = 3) = 1/6$ |
/// | **PDF** | ProbabilityDii znsity Function   | ❌ 没有      | ✅ 有       | 连续型的“概率密度”，不能直接代表概率             |
/// | **CDF** | Cumulative Distribution Func.  | ✅ 有       | ✅ 有       | 累积分布函数： $P(X \le x)$            |
/// | **PPF** | Percent Point Function (CDF⁻¹) | ✅ 有       | ✅ 有       | 反函数，给定概率 $p$，返回 $x$             |
/// | **SF**  | Survival Function              | ✅ 有       | ✅ 有       | 尾部概率： $P(X > x) = 1 - CDF(x)$   |
/// | **ISF** | Inverse Survival Function      | ✅ 有       | ✅ 有       | SF 的反函数： $x = \text{ISF}(p)$    |
/// | **RNG** | Random Number Generator (采样)   | ✅ 有       | ✅ 有       | 从该分布中随机采样值                      |

