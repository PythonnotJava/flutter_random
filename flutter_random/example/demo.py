import numpy as np
from scipy.integrate import quad

def f(x):
    return np.exp(-x**2) * np.sin(10*x)

result, error = quad(f, 0, 2)
print(f"积分结果：{result}, 误差估计：{error}")
print(result - 0.10088272642690274)